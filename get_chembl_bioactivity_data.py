"""
get_chembl_bioactivity_data.py

Reusable function based on TeachOpenCADD T001 (Compound data acquisition - ChEMBL):
https://projects.volkamerlab.org/teachopencadd/talktorials/T001_query_chembl.html

Given a UniProt ID, fetches the corresponding ChEMBL target, pulls IC50 bioactivity
data (human, exact measurements, binding assays), links it to compound SMILES, and
returns a tidy pd.DataFrame with pIC50 values.
"""

from __future__ import annotations

import math

import pandas as pd
from chembl_webresource_client.new_client import new_client
from tqdm.auto import tqdm


def convert_ic50_to_pic50(ic50_value: float) -> float:
    """Convert an IC50 value in nM to a pIC50 value."""
    return 9 - math.log10(ic50_value)


def fetch_chembl_targets(uniprot_id: str) -> pd.DataFrame:
    """
    Query ChEMBL for all target(s) matching a UniProt accession and display them.

    A UniProt accession can map to several ChEMBL target types (single protein,
    protein family, chimeric construct, protein-protein interaction, etc.), so
    inspect this table and pick the row index of the target you actually want
    before calling get_chembl_bioactivity_data(..., target_index=<your choice>).

    Parameters
    ----------
    uniprot_id : str
        UniProt accession of the target of interest, e.g. "P00533" (EGFR).

    Returns
    -------
    pd.DataFrame
        Columns: organism, pref_name, target_chembl_id, target_type.
        The DataFrame index is what you pass as `target_index` downstream.

    Raises
    ------
    ValueError
        If no ChEMBL target is found for the given UniProt ID.
    """
    targets_api = new_client.target
    targets = targets_api.get(target_components__accession=uniprot_id).only(
        "target_chembl_id", "organism", "pref_name", "target_type"
    )
    targets_df = pd.DataFrame.from_records(targets)

    if targets_df.empty:
        raise ValueError(f"No ChEMBL target found for UniProt ID '{uniprot_id}'.")

    print(f"Found {targets_df.shape[0]} ChEMBL target(s) matching '{uniprot_id}':\n")
    print(targets_df.to_string())
    print(
        "\nInspect the table above, then call get_chembl_bioactivity_data("
        f"'{uniprot_id}', target_index=<row>) with your chosen row index."
    )

    return targets_df


def get_chembl_bioactivity_data(
    uniprot_id: str,
    target_index: int | None = None,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Fetch bioactivity + compound data from ChEMBL for a given target (by UniProt ID).

    Follows the same query/filter/merge logic as TeachOpenCADD talktorial T001:
    - Look up ChEMBL target(s) matching the UniProt accession
    - Fetch IC50 bioactivities (relation '=', assay_type 'B') for the selected target
    - Filter to nM units, drop NaNs/duplicates
    - Fetch compound structures (canonical SMILES) for the resulting molecules
    - Merge bioactivity + compound data and compute pIC50

    Recommended workflow (no blocking input() prompts, safe for any environment
    including agent-driven IDEs where stdin isn't interactive):

        targets_df = fetch_chembl_targets("P00533")   # inspect the printed table
        df = get_chembl_bioactivity_data("P00533", target_index=0)  # pick a row

    Parameters
    ----------
    uniprot_id : str
        UniProt accession of the target of interest, e.g. "P00533" (EGFR).
    target_index : int or None, default None
        Row index into the queried targets DataFrame to use.
        - If None: falls back to the first row where target_type == "SINGLE PROTEIN"
          and organism == "Homo sapiens" (or row 0 if no such match exists), and
          prints which target was auto-selected so you can verify it.
        - If an int is given, that index is used directly. Run fetch_chembl_targets()
          first to see the available indices for this UniProt ID.
    show_progress : bool, default True
        Whether to show a tqdm progress bar while downloading compound records.

    Returns
    -------
    pd.DataFrame
        Columns: molecule_chembl_id, IC50, units, smiles, pIC50
        (IC50 in nM; pIC50 computed as 9 - log10(IC50)).

    Raises
    ------
    ValueError
        If no ChEMBL target is found for the given UniProt ID.
    """
    # --- API resource objects ---
    targets_api = new_client.target
    compounds_api = new_client.molecule
    bioactivities_api = new_client.activity

    # --- Get target data ---
    targets = targets_api.get(target_components__accession=uniprot_id).only(
        "target_chembl_id", "organism", "pref_name", "target_type"
    )
    targets_df = pd.DataFrame.from_records(targets)

    if targets_df.empty:
        raise ValueError(f"No ChEMBL target found for UniProt ID '{uniprot_id}'.")

    # --- Resolve which target row to use ---
    if target_index is not None:
        chosen_index = target_index
    else:
        single_protein = targets_df[
            (targets_df["target_type"] == "SINGLE PROTEIN")
            & (targets_df["organism"] == "Homo sapiens")
        ]
        if not single_protein.empty:
            chosen_index = single_protein.index[0]
        else:
            chosen_index = 0
            print(
                f"Warning: no 'SINGLE PROTEIN' human target found for '{uniprot_id}'. "
                f"Defaulting to row 0:\n{targets_df.iloc[0]}"
            )

    target = targets_df.iloc[chosen_index]
    chembl_id = target.target_chembl_id
    print(
        f"Using target: {target.pref_name} ({chembl_id}, "
        f"{target.target_type}, {target.organism})\n"
    )

    # --- Fetch bioactivity data ---
    bioactivities = bioactivities_api.filter(
        target_chembl_id=chembl_id, type="IC50", relation="=", assay_type="B"
    ).only(
        "activity_id",
        "assay_chembl_id",
        "assay_description",
        "assay_type",
        "molecule_chembl_id",
        "type",
        "standard_units",
        "relation",
        "standard_value",
        "target_chembl_id",
        "target_organism",
    )
    bioactivities_df = pd.DataFrame.from_dict(bioactivities)

    if bioactivities_df.empty:
        raise ValueError(f"No IC50 bioactivity data found for target '{chembl_id}'.")

    # Drop non-standardized columns (present when using .filter without .only trimming units/value)
    bioactivities_df = bioactivities_df.drop(columns=["units", "value"], errors="ignore")

    # --- Preprocess and filter bioactivity data ---
    bioactivities_df = bioactivities_df.astype({"standard_value": "float64"})
    bioactivities_df.dropna(axis=0, how="any", inplace=True)
    bioactivities_df = bioactivities_df[bioactivities_df["standard_units"] == "nM"]
    bioactivities_df.drop_duplicates("molecule_chembl_id", keep="first", inplace=True)
    bioactivities_df.reset_index(drop=True, inplace=True)
    bioactivities_df.rename(
        columns={"standard_value": "IC50", "standard_units": "units"}, inplace=True
    )

    # --- Fetch compound data ---
    compounds_provider = compounds_api.filter(
        molecule_chembl_id__in=list(bioactivities_df["molecule_chembl_id"])
    ).only("molecule_chembl_id", "molecule_structures")

    compounds = list(tqdm(compounds_provider, disable=not show_progress))
    compounds_df = pd.DataFrame.from_records(compounds)

    # --- Preprocess and filter compound data ---
    compounds_df.dropna(axis=0, how="any", inplace=True)
    compounds_df.drop_duplicates("molecule_chembl_id", keep="first", inplace=True)

    canonical_smiles = []
    for _, compound in compounds_df.iterrows():
        try:
            canonical_smiles.append(compound["molecule_structures"]["canonical_smiles"])
        except KeyError:
            canonical_smiles.append(None)
    compounds_df["smiles"] = canonical_smiles
    compounds_df.drop("molecule_structures", axis=1, inplace=True)
    compounds_df.dropna(axis=0, how="any", inplace=True)

    # --- Merge bioactivity and compound data ---
    output_df = pd.merge(
        bioactivities_df[["molecule_chembl_id", "IC50", "units"]],
        compounds_df,
        on="molecule_chembl_id",
    )
    output_df.reset_index(drop=True, inplace=True)

    # --- Add pIC50 values ---
    output_df["pIC50"] = output_df["IC50"].apply(convert_ic50_to_pic50)

    # --- Sort by potency (highest pIC50 first) ---
    output_df.sort_values(by="pIC50", ascending=False, inplace=True)
    output_df.reset_index(drop=True, inplace=True)

    return output_df[["molecule_chembl_id", "IC50", "units", "smiles", "pIC50"]]


if __name__ == "__main__":
    # Step 1: fetch and inspect candidate targets for the UniProt ID
    targets_df = fetch_chembl_targets("P00533")  # EGFR

    # Step 2: run extraction using the row index you want (inspect printed
    # table above first; here we use 0 as in the original T001 talktorial)
    df = get_chembl_bioactivity_data("P00533", target_index=0)
    print(f"Retrieved {df.shape[0]} compounds with bioactivity data.")
    print(df.head())

    # If you skip target_index, it auto-selects the first SINGLE PROTEIN /
    # Homo sapiens match (with a printed warning if none exists):
    # df = get_chembl_bioactivity_data("P00533")
