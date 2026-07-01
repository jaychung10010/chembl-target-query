# get-chembl-bioactivity-data

Small, reusable Python utility for pulling target-specific bioactivity and compound
data from [ChEMBL](https://www.ebi.ac.uk/chembl/) given a UniProt ID — adapted from
[TeachOpenCADD T001: Compound data acquisition (ChEMBL)](https://projects.volkamerlab.org/teachopencadd/talktorials/T001_query_chembl.html).

Given a UniProt accession (e.g. `P00533` for EGFR), it:

1. Looks up the matching ChEMBL target(s)
2. Fetches IC50 bioactivity data (human, exact measurements, binding assays)
3. Filters to nM units and deduplicates by compound
4. Fetches canonical SMILES for the resulting compounds
5. Merges bioactivity + compound data and computes pIC50

... and returns a tidy `pandas.DataFrame` with columns:

| molecule_chembl_id | IC50 | units | smiles | pIC50 |
|---|---|---|---|---|

## Installation

```bash
git clone https://github.com/<your-username>/get-chembl-bioactivity-data.git
cd get-chembl-bioactivity-data
pip install -r requirements.txt
```

## Usage

A UniProt accession can map to multiple ChEMBL target entries (single protein,
protein family, chimeric construct, protein-protein interaction, etc.), so the
workflow is split into two steps — inspect, then extract. This avoids blocking
`input()` prompts, so it works in any environment (plain scripts, Jupyter,
agent-driven IDEs like Antigravity/Cursor where stdin isn't interactive).

```python
from get_chembl_bioactivity_data import fetch_chembl_targets, get_chembl_bioactivity_data

# Step 1: see what ChEMBL targets match this UniProt ID
targets_df = fetch_chembl_targets("P00533")  # EGFR

# Step 2: extract bioactivity + compound data for the target you want
# (inspect targets_df above, then pick its row index)
df = get_chembl_bioactivity_data("P00533", target_index=0)

df.head()
```

**Skipping `target_index`** falls back to auto-selecting the first
`SINGLE PROTEIN` + `Homo sapiens` match (printing a warning if none exists) —
useful for unattended/batch runs over many targets:

```python
targets = ["P00533", "Q00534", "P07900"]  # example UniProt IDs
results = {uid: get_chembl_bioactivity_data(uid) for uid in targets}
```

## API

### `fetch_chembl_targets(uniprot_id: str) -> pd.DataFrame`
Queries and prints all ChEMBL targets matching a UniProt accession. Returns the
DataFrame so you can inspect `target_type`, `organism`, and `pref_name` before
choosing which row to extract.

### `get_chembl_bioactivity_data(uniprot_id, target_index=None, show_progress=True) -> pd.DataFrame`
Runs the full extraction pipeline for the selected target and returns the merged,
filtered bioactivity + compound DataFrame with pIC50 values.

### `convert_ic50_to_pic50(ic50_value: float) -> float`
Converts an IC50 value in nM to pIC50 (`9 - log10(IC50)`).

## Notes

- Query speed depends heavily on how much bioactivity data exists for the target
  (a well-studied kinase like EGFR can take from several seconds up to ~2 minutes)
  and on EBI server load — there's no documented SLA for the public ChEMBL API.
- This follows the exact filtering logic from TeachOpenCADD T001: IC50 measurements
  only, exact relation (`=`), binding assays (`B`), nM units, first-seen compound
  kept on duplicates.

## Acknowledgments

Built on the [`chembl_webresource_client`](https://github.com/chembl/chembl_webresource_client)
and adapted from the [TeachOpenCADD](https://github.com/volkamerlab/teachopencadd)
platform (Volkamer Lab, Charité/FU Berlin).

## License

MIT — see [LICENSE](LICENSE).
