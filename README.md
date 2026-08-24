# WGCNA Co-expression Network Pipeline

*EuchroGene WGCNA v2.1 — for EuchroGene members.*

A containerised weighted gene co-expression network analysis that takes an expression matrix and returns modules, phenotype associations, hub genes, functional enrichment and submission-ready figures. Networks are built with [WGCNA](https://doi.org/10.1186/1471-2105-9-559), counts are normalised with [DESeq2](https://doi.org/10.1186/s13059-014-0550-8), batch effects are removed with [limma](https://doi.org/10.1093/nar/gkv007) or [ComBat](https://doi.org/10.1093/biostatistics/kxj037), Gene Ontology enrichment runs through [clusterProfiler](https://doi.org/10.1016/j.xinn.2021.100141), and network figures are drawn with [igraph](https://igraph.org) alongside [Cytoscape](https://doi.org/10.1101/gr.1239303) tables.

Figures are written at the final printed size with embedded fonts, no text below 5 pt and collision-avoided gene labels, so they can go into a manuscript without being redrawn. The HTML report carries a Methods section generated from the settings the run actually used, with the tool versions queried at run time.

---

## Installation

**0. Install EG_tools**

```bash
wget https://github.com/euchrogene/EG_tools/raw/main/EG_tools
chmod 777 EG_tools
sudo mv EG_tools /usr/bin/
```

**1. Remove a previous version**

Check what is installed first, so the executable name and image tag are the ones
your machine actually has:

```bash
EG_tools
```

If v.2.0 is listed, remove it before installing v2.1. Note the dots: the old
names carry them, the new ones do not.

```bash
sudo EG_tools uninstall -t WGCNA_v.2.0 -i managene7/wgcna_package:v.2.0
```

Confirm the old executable is gone and the old image has been released:

```bash
which WGCNA_v.2.0                            # should print nothing
docker images | grep wgcna_package
```

If the old image is still listed, remove it by hand and reclaim the space:

```bash
docker rmi managene7/wgcna_package:v.2.0
```

Results from earlier runs are not touched. Reports written by v.2.0 remain
readable, but their Methods text and figures were produced by the older code,
so do not mix them with v2.1 output in one manuscript.

**2. Install**

```bash
sudo EG_tools install -r https://github.com/euchrogene/WGCNA.git \
                      -d WGCNA \
                      -e WGCNA_v2.1 \
                      -m "Weighted gene co-expression network analysis"
```

**3. Display installed software**

```bash
EG_tools
```

**4. Show help contents**

```bash
WGCNA_v2.1
```

**5. Uninstall**

```bash
sudo EG_tools uninstall -t WGCNA_v2.1 -i managene7/wgcna_package:v2.1
```

Docker image: `managene7/wgcna_package:v2.1`

Both halves have to move together on an upgrade. The seven analysis scripts live
inside the image and all seven changed in v2.1, one of them newly added, so a
v2.1 wrapper against a v.2.0 image is refused before the first stage:

```
[Error] The image does not contain the expected workers: export_multi_network.R.
[Error] Rebuild and reinstall the image; do not edit the R scripts on the host.
```

---

## Quick Start

**1. Minimum run, counts only**

```bash
WGCNA_v2.1 -i counts.csv -o wgcna_results
```

**2. With phenotypes, the usual case**

```bash
WGCNA_v2.1 -i counts.csv \
           -traits phenotypes.csv \
           -species "Vaccinium virgatum" \
           -o wgcna_results
```

**3. Publication run: enrichment, robustness, provenance**

```bash
WGCNA_v2.1 -i counts.csv \
           -traits phenotypes.csv \
           -emapper eggnog_annotations.tsv \
           -species "Vaccinium virgatum" \
           -assembly VacVir_v2 \
           -accession PRJNA1024567 \
           -stability true \
           -o wgcna_results
```

**4. Batch effect removed, biology preserved**

```bash
WGCNA_v2.1 -i counts.csv \
           -traits phenotypes.csv \
           -covariates sample_metadata.csv \
           -batch sequencing_run \
           -batchprotect tissue \
           -o wgcna_results
```

**5. Figures sized for a specific journal**

```bash
WGCNA_v2.1 -i counts.csv -traits phenotypes.csv \
           -figwidthmm 89 -basept 8 -fillby gs \
           -o wgcna_results
```

**6. Phenotype templates, and a check before a long run**

```bash
WGCNA_v2.1 template -i counts.csv -o phenotype_templates
WGCNA_v2.1 template -i counts.csv -traits my_phenotypes.csv
```

The first writes a file per phenotype form plus a format guide; the second
reports how every column will be read and whether the sample IDs match.

**7. Re-run one stage against existing results**

```bash
WGCNA_v2.1 export -o wgcna_results -traits phenotypes.csv -topk 50
```

Modes: `all` (default, may be omitted), `run`, `traits`, `hubs`, `export`, `enrich`, `stability`, `template`.

---

## Inputs

| Input | Required | Notes |
|---|---|---|
| `-i` expression matrix | yes | CSV, genes in rows, samples in columns, first column the gene ID. Raw counts by default; use `-type tpm` or `-type fpkm` for normalised input. |
| `-traits` phenotype table | no | CSV, first column the sample ID matching the expression header. Continuous, ordinal, binary and multi-level categorical columns are all accepted, as numbers or as text; the encoding applied to each is recorded. Generate a template and check your file with the `template` mode. |
| `-traitlog2` | `none` | log2-transform trait columns before analysis: `all` (every measured continuous column) or a comma-separated list. Missing cells (blank, `NA`, `ND`, `N/A`, `NULL`, `#N/A`) stay missing; zeros or negatives stop the run (replace below-LOQ zeros with LOD/2 or blank them). |
| `-traittypes` declarations | no | CSV of `trait,type` with type in `continuous`, `ordinal`, `binary`, `categorical`, `ignore`. Only needed to mark ordinal columns and to drop identifiers or free text. |
| `-covariates` metadata | no | CSV in the same shape as `-traits`. Defaults to the trait file. Supplies the batch column for `-batch`. |
| `-emapper` / `-interpro` / `-annotation` | no | Functional annotation for GO enrichment of non-model species. `-orgdb` covers Arabidopsis, human and mouse instead. |
| `-testexpr` independent matrix | no | Enables `modulePreservation` and a Zsummary per module. |
| `-geneannot` gene labels | no | Two columns, gene ID and display label, used for figure labels in place of the raw IDs. |

Absolute paths, relative paths and directories outside the working directory all resolve; the wrapper mounts what each stage needs.

---

## Outputs

```
wgcna_results/
├── 0_WGCNA_Report.html              # main report: Methods, tables, figures, references
├── RUN_SPEC.json                    # parameters, versions, worker checksums, input checksums
├── wgcna_modules.csv                # gene to module assignment  <- use this downstream
├── wgcna_normalized_expression.csv  # matrix the network was built from
├── wgcna_analysis_params.csv        # every network parameter as applied
├── wgcna_hub_genes.csv              # hub genes with kME, kWithin, GS and FDR
├── wgcna_hub_params.csv             # the hub criterion and the counts it produced
├── wgcna_module_trait_stats.csv     # module x trait: correlation, effect size, p, FDR
├── wgcna_trait_encoding.csv         # how each trait column was interpreted and tested
├── wgcna_correlations.csv           # module x trait correlation matrix
├── wgcna_GO.csv                     # enriched GO terms per module
├── wgcna_enrichment_params.csv      # cutoffs, ontology, redundancy collapse
├── wgcna_module_stability.csv       # bootstrap Jaccard per module
├── wgcna_module_preservation.csv    # Zsummary per module, with -testexpr
├── wgcna_stability_params.csv       # the resampling design
├── wgcna_multi_network.pdf/.png     # one network across all associated modules
├── wgcna_<module>_network_plot.pdf  # per-module network figure
├── wgcna_<module>_network_params.csv
├── wgcna_<module>_edges.txt         # Cytoscape edge table
├── wgcna_<module>_nodes.txt         # Cytoscape node table
├── wgcna_<module>_network.graphml
├── wgcna_sft_plot.pdf               # soft-threshold selection
├── wgcna_dendrogram.pdf             # gene dendrogram with module colours
├── wgcna_TOM_heatmap.pdf
├── wgcna_sample_PCA.pdf             # hidden structure and covariate testing
├── wgcna_trait_heatmap.pdf          # module-trait relationships
├── wgcna_eigengene_vs_trait.pdf     # the data behind each significant association
└── wgcna_MM_vs_GS_scatterplots.pdf
```

Every figure is written as a vector PDF for submission and a PNG for the report.

---

## Key Options

| Option | Default | Purpose |
|---|---|---|
| `-i` | (required) | Expression matrix |
| `-o` | `wgcna_results` | Output directory |
| `-type` | `counts` | `counts` (VST) \| `tpm` \| `fpkm` (log2) |
| `-traits` | none | Phenotype table; enables module-trait analysis |
| `-traittypes` | none | Column type declarations; needed only for ordinal columns and to drop identifiers |
| `-species` | none | Organism named in the generated Methods |
| `-assembly` | none | Reference assembly named in the Methods |
| `-accession` | none | Data accession stated in the Methods |
| `-p` | `16` | Threads; a restrained default on a 32-core host |
| `-max_mem` | `150` | Memory budget in GB; also caps the block size |
| `-network` | `signed hybrid` | `unsigned` \| `signed` \| `signed hybrid` |
| `-tomtype` | `signed` | Topological overlap type |
| `-power` | auto | Soft-thresholding power; auto uses the scale-free criterion |
| `-minmodule` | `30` | Minimum module size |
| `-merge` | `0.25` | Eigengene dissimilarity cut height for merging |
| `-batch` | none | Covariate column whose effect is removed before the network is built |
| `-batchprotect` | none | Columns preserved during that removal; name your biological variable here |
| `-batchmethod` | `removebatcheffect` | `removebatcheffect` (limma) \| `combat` (sva) |
| `-hubcriterion` | `auto` | `mm_gs` when traits are given, else `mm`; the hub count follows from the criterion |
| `-kmethreshold` | `0.8` | Minimum \|kME\| for a hub |
| `-gsthreshold` | `0.2` | Minimum \|GS\| for the combined criterion |
| `-hubfdr` | `0.05` | FDR applied to kME and GS |
| `-exportby` | `auto` | `trait` exports every phenotype-associated module regardless of size |
| `-traitfdr` | `0.05` | FDR that defines "associated" |
| `-topk` | `30` | Genes drawn per module network; a floor, hubs are added on top |
| `-nhub` | `0` | `0` marks every gene meeting the hub criterion; a positive value is an arbitrary cut |
| `-labeltop` | `-1` | How many hubs carry a text label; markers stay on the rest |
| `-fillby` | `auto` | `gs` (gene significance) for phenotype-selected modules, `kme` otherwise; either can be forced |
| `-figwidthmm` | `180` | Final printed width; `89` single column, `174` for Cell |
| `-basept` | `7` | Body text size in points; use `8` for the PNAS 6 pt floor |
| `-palette` | `okabe_ito` | Colour-blind safe module colours in the combined figure |
| `-combined` | `true` | Draw one network across all exported modules |
| `-stability` | `false` | Bootstrap module robustness; costs roughly `-nboot` times stage 1 |
| `-nboot` | `50` | Bootstrap replicates |
| `-emapper` | none | eggNOG-mapper annotations for GO enrichment |
| `-orgdb` | none | Built-in database: `arabidopsis` \| `human` \| `mouse` |
| `-gosimplify` | `true` | Collapse redundant GO terms |
| `-seed` | `12345` | Random seed, recorded in the Methods |

Run `WGCNA_v2.1` with no arguments for the full option list.

---

## Notes

**Sample size.** WGCNA is unreliable below roughly 15 samples: correlation estimates are noisy and small modules may not reproduce. The pipeline warns below that threshold and the generated Methods says so, which is better raised by the authors than by a reviewer. With few samples, run `-stability true`.

**Batch effects.** `-covariates` alone only tests the principal components against each covariate and reports what it finds. Removal happens only when `-batch` names a column. Give `-batchprotect` as well: without it the batch is removed blind to the biology, and where the two are correlated the signal goes with it. If they are fully confounded the run stops rather than producing a silently gutted network.

**Phenotype forms.** Continuous, ordinal, binary and multi-level categorical columns are all accepted, written as numbers or as text. A column with three or more text levels is expanded into one indicator per level, each tested as a point-biserial correlation, and is additionally tested for an overall effect on the eigengene by Kruskal-Wallis; quote the omnibus test for association with the trait as a whole. Every correlation and omnibus test shares one Benjamini-Hochberg family.

**Ordinal is the one form that cannot be detected.** A 0-5 severity score and a measurement in grams both look like "numeric with several values". Left alone, an ordinal column is tested with Pearson, which assumes the steps are equally spaced. Declare it in `-traittypes` and it is tested with Spearman instead, which uses only the order. Declare `ignore` for plot numbers, dates and free-text notes, which would otherwise be read as categorical variables with as many levels as there are samples.

**Check the file first.** `WGCNA_v2.1 template -i counts.csv -traits your_file.csv` reports how each column will be read and whether the sample IDs match. A pipeline run stops before any stage starts when the IDs do not match, listing every unmatched name on both sides together with the likely fix (a case variant on the other side, or a summary column such as `avg` sitting in the expression matrix). Fix the names and rerun; nothing is computed until they agree.

**Hub genes.** The count is a result of the criterion, not a setting, and the figures show the set the criterion defines rather than a fixed number. If a module returns more hubs than a panel can carry, tighten the criterion with `-kmethreshold` (0.85 or 0.9), which is a rule a reader can evaluate, rather than capping the count with `-nhub`, which is not. `-labeltop` thins the text labels without changing which genes are marked as hubs.

**Blocking.** When the gene set exceeds the block size, module detection runs in several blocks and genes in different blocks cannot join the same module. The report states the block count and the limitation when it applies.

**Figures.** Written at the width given by `-figwidthmm` so the journal does not reduce them, with fonts embedded and a 5 pt floor. Check `-figwidthmm 174` for Cell and `-basept 8` for PNAS.

**Stale images.** The workers live inside the image. The wrapper verifies all seven before the first stage and records their SHA-256 checksums in `RUN_SPEC.json` and the report, so a mismatch between wrapper and image fails immediately instead of surfacing later as an unexplained argument error. If it reports a missing worker, rebuild and reinstall the image.

---

## Citation

> Langfelder, P. and Horvath, S. (2008) WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 9: 559.

> Zhang, B. and Horvath, S. (2005) A general framework for weighted gene co-expression network analysis. *Statistical Applications in Genetics and Molecular Biology* 4: Article 17.

> Langfelder, P., Zhang, B. and Horvath, S. (2008) Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. *Bioinformatics* 24: 719-720.

> Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B* 57: 289-300.

> Csardi, G. and Nepusz, T. (2006) The igraph software package for complex network research. *InterJournal, Complex Systems* 1695.

> Shannon, P., Markiel, A., Ozier, O., et al. (2003) Cytoscape: a software environment for integrated models of biomolecular interaction networks. *Genome Research* 13: 2498-2504.

> Love, M.I., Huber, W. and Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15: 550. — *when -type counts is used*

> Ritchie, M.E., Phipson, B., Wu, D., et al. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 43: e47. — *when -batch is used*

> Johnson, W.E., Li, C. and Rabinovic, A. (2007) Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics* 8: 118-127. — *when -batchmethod combat is used*

> Ashburner, M., Ball, C.A., Blake, J.A., et al. (2000) Gene Ontology: tool for the unification of biology. *Nature Genetics* 25: 25-29. — *when enrichment is run*

> Wu, T., Hu, E., Xu, S., et al. (2021) clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *The Innovation* 2: 100141. — *when enrichment is run*

> Cantalapiedra, C.P., Hernandez-Plaza, A., Letunic, I., et al. (2021) eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Molecular Biology and Evolution* 38: 5825-5829. — *when -emapper is used*

> Jones, P., Binns, D., Chang, H.-Y., et al. (2014) InterProScan 5: genome-scale protein function classification. *Bioinformatics* 30: 1236-1240. — *when -interpro is used*

> Hennig, C. (2007) Cluster-wise assessment of cluster stability. *Computational Statistics and Data Analysis* 52: 258-271. — *when -stability true is used*

> Langfelder, P., Luo, R., Oldham, M.C. and Horvath, S. (2011) Is my network module preserved and reproducible? *PLoS Computational Biology* 7: e1001057. — *when -testexpr is used*

> EuchroGene WGCNA Pipeline v2.1 (2026). EuchroGene, LLC.

The Methods section inside the HTML report is generated from the actual run settings and the tool versions queried at run time, ready to copy into a manuscript.

**Support:** bioinformatics@euchrogene.com
