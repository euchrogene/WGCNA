# WGCNA Co-expression Network Analysis

*EuchroGene WGCNA v2.0 — for EuchroGene members.*

Publication-grade weighted gene co-expression network analysis built on [WGCNA](https://doi.org/10.1186/1471-2105-9-559). Takes a gene expression matrix — the count matrices from the EuchroGene Salmon and STAR/RSEM pipelines drop straight in — and runs the full workflow: filtering and variance-stabilizing normalization, batch-effect diagnostics, network construction, module detection, hub genes with explicit statistical criteria, FDR-controlled module–trait association, bootstrap module robustness, GO enrichment for any species, and Cytoscape-ready export.

Every parameter is recorded at run time and the methods text is generated from those records, so the manuscript description always matches what was executed.

---

## Installation

### 0. Install EG_tools &nbsp; *(skip if already installed)*
```
wget https://github.com/euchrogene/EG_tools/raw/refs/heads/main/EG_tools
sudo chmod 777 EG_tools
sudo mv EG_tools /usr/bin
```

### 1. Install
```
sudo EG_tools install -r https://github.com/euchrogene/WGCNA.git -d WGCNA -e WGCNA_v.2.0 -m "End-to-end weighted gene co-expression network analysis"
```

### 2. Display installed software
```
EG_tools
```

### 3. Show help contents
```
WGCNA_v.2.0
```

### 4. Uninstall
```
sudo EG_tools uninstall -t WGCNA_v.2.0 -i managene7/wgcna_package:v.2.0
```

Docker image: `managene7/wgcna_package:v.2.0`

---

## Quick Start

```bash
# Standard run (the mode is optional - 'all' is the default)
WGCNA_v.2.0 -i gene_counts.csv -traits phenotypes.csv

# Publication run: module robustness + GO enrichment for any plant
WGCNA_v.2.0 all -i gene_counts.csv -traits traits.csv \
    -stability true -emapper species.emapper.annotations

# Tissue atlas: keep tissue-specific genes, export every module
WGCNA_v.2.0 all -i gene_counts.csv -traits tissue_indicators.csv \
    -minfrac 0.08 -n all

# Not sure how to write the trait file? Generate templates from your own samples
WGCNA_v.2.0 template -i gene_counts.csv -o trait_templates
```

Run `WGCNA_v.2.0 -help` for the full option list.
Modes are optional — with none given the complete pipeline (`all`) runs. Other modes: `run`, `traits`, `hubs`, `export`, `enrich`, `stability`, `template`.

---

## Inputs

| Input | Required | Notes |
|---|---|---|
| `-i` | yes | Expression matrix: genes as rows (first column = gene IDs), one column per sample |
| `-traits` | recommended | One row per sample, first column = sample ID matching the matrix header; **all other columns numeric** (encode categories as 0/1) |
| `-covariates` | optional | Batch, run, library prep… for batch-effect diagnostics (defaults to `-traits`) |

> `WGCNA_v.2.0 template -i gene_counts.csv` writes ready-to-edit trait and covariate templates pre-filled with your sample names, plus a one-page format guide.

> Supply **counts** with `-type counts`. VST is the correct normalization for count data; a TPM matrix declared as counts is rejected rather than silently transformed.

> WGCNA needs enough samples for stable correlations — at least 15, preferably 20+.

---

## Outputs

```
<out>/
├── 0_WGCNA_Report.html              main report: methods, figures, tables
├── Methods_Publication.txt          methods section, ready to paste
├── <p>_modules.csv                  gene → module assignment
├── <p>_hub_genes.csv                hubs: kME, p, FDR, kWithin, gene significance
├── <p>_module_trait_stats.csv       module × trait: correlation, p, FDR
├── <p>_module_stability.csv         bootstrap Jaccard            [-stability]
├── <p>_PC_covariate_association.csv batch diagnostics
├── <p>_GO.csv                       GO enrichment per module
├── <p>_analysis_params.csv          every parameter and version used
├── Figures                          PDF (submission) + 300 dpi PNG
└── per module: Cytoscape edges/nodes, visual style, GraphML, network figure
```

---

## Key Options

| Option | Default | Purpose |
|---|---|---|
| `-o` / `-type` | `wgcna_results` / `counts` | Output folder / data type (`counts`, `tpm`, `fpkm`) |
| `-p` / `-max_mem` | `30` / `150` | CPU threads / memory budget in GB |
| `-minfrac` | `0.2` | Fraction of samples a gene must pass; lower it for tissue atlases |
| `-network` | `signed hybrid` | `unsigned` \| `signed` \| `signed hybrid` |
| `-power` | auto | Manual soft-thresholding power |
| `-kmethreshold` / `-hubtopn` | `0.8` / `10` | Hub criteria: minimum \|kME\| and max hubs per module |
| `-gsthreshold` | `0.2` | Combined MM + GS criterion (needs `-traits`) |
| `-stability` / `-nboot` | `false` / `50` | Bootstrap module robustness |
| `-testexpr` | — | Independent matrix → `modulePreservation` (Zsummary) |
| `-n` / `-preselect` | `20` / `500` | Modules exported / genes per module before the TOM |
| `-seed` | `12345` | Random seed (block assignment is randomised) |

**GO enrichment** works for any species — supply one of `-emapper` (eggNOG-mapper), `-interpro` (InterProScan), `-annotation` (`GeneID,GO_ID`), or `-species` (`arabidopsis`|`human`|`mouse`).

---

## Notes

- **Memory scales with the square of the block size.** The default 30,000-gene block needs roughly 33 GB; `-max_mem` caps the block size accordingly and the estimated peak is printed before construction.
- **Hub genes need a threshold, not just a rank.** Genes must clear \|kME\| ≥ 0.8 at FDR < 0.05 before the strongest per module are reported, and intramodular connectivity (kWithin) is reported alongside kME. The same set drives the hub table, the network figure and the Cytoscape export.
- **Module–trait p-values are FDR-corrected** across the full module × trait matrix — the FDR, not the raw p-value, drives every downstream selection.
- **Batch effects are quantified, not assumed away.** Each leading principal component is tested against every covariate, so hidden structure shows up before it reaches the modules.
- `-stability true` rebuilds the network on resampled samples (50× by default) and is slow — use it for the final publication run.

---

## Citation

> Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 9: 559.

> Zhang B, Horvath S (2005) A general framework for weighted gene co-expression network analysis. *Stat Appl Genet Mol Biol* 4: Article 17.

> Langfelder P, Zhang B, Horvath S (2008) Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. *Bioinformatics* 24(5): 719–720.

> Love MI, Huber W, Anders S (2014) DESeq2. *Genome Biology* 15: 550. — *when counts are normalized with VST*

> Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate. *J R Stat Soc B* 57(1): 289–300.

> Hennig C (2007) Cluster-wise assessment of cluster stability. *Comput Stat Data Anal* 52(1): 258–271. — *when bootstrap stability is reported*

> Langfelder P, Luo R, Oldham MC, Horvath S (2011) Is my network module preserved and reproducible? *PLoS Comput Biol* 7(1): e1001057. — *when module preservation is reported*

> Wu T, Hu E, Xu S, et al. (2021) clusterProfiler 4.0. *The Innovation* 2(3): 100141. — *when GO enrichment is used*

> Shannon P, Markiel A, Ozier O, et al. (2003) Cytoscape. *Genome Research* 13(11): 2498–2504. — *when networks are visualized in Cytoscape*

> EuchroGene WGCNA Pipeline v2.0 (2026). EuchroGene, LLC.

`Methods_Publication.txt` and the methods inside `0_WGCNA_Report.html` are parameterized on the actual run settings and list only the references they cite.

**Support:** bioinformatics@euchrogene.com
