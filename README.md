# Human kidney slice culture: figure and table reproduction

Code accompanying the manuscript **Human Kidney Precision Cut Slice Cultures Are a Flexible Platform to Study Mechanisms of Glomerular Disease** by Shreeram Akilesh and colleagues.

This repository reproduces the GeoMx spatial transcriptomic figures and supplemental tables for human precision cut kidney slice cultures. Analyses examine culture medium and oxygenation, transcriptional changes over 72 hours, IFNγ responses, and modulation by the JAK1/2 inhibitor baricitinib in glomerular and tubulointerstitial compartments.

The output numbering follows the October 9, 2026 manuscript, with the subsequent corrected supplemental table order.

## Repository contents

- `reproduce_figures_tables.R`: figure and table generation script.
- `GSE345813_GeoMx_Analysis_Complete_Workspace.RData.gz`: required saved analysis workspace, hosted on GEO and downloaded separately (not committed to this repository).
- `LICENSE`: MIT license for the code.

The script loads a previously prepared workspace. It is a reproduction/export pipeline, not a complete workflow starting from raw sequencing reads. Some analyses are recomputed from workspace expression matrices and metadata; other outputs use saved analysis results.

## Requirements

- R 4.5.1 or a compatible version. The pipeline has been run successfully in the author's local R environment; compatibility with other versions has not been systematically tested.
- R packages: `ggplot2`, `pheatmap`, `patchwork`, `reshape2`, `limma`, `fgsea`, `ggrepel`, and `openxlsx`.
- The saved workspace, approximately 234 MB.

Install packages from an R console:

```r
install.packages(c("ggplot2", "pheatmap", "patchwork", "reshape2",
                   "ggrepel", "openxlsx", "BiocManager"))
BiocManager::install(c("limma", "fgsea"))
```

The saved workspace contains precomputed ssGSEA scores. Running this export script does not require recomputing them with GSVA.

## Obtain the workspace

Download **`GSE345813_GeoMx_Analysis_Complete_Workspace.RData.gz`** from
[GEO accession GSE345813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345813).

[Direct workspace download](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE345813&format=file&file=GSE345813%5FGeoMx%5FAnalysis%5FComplete%5FWorkspace%2ERData%2Egz)

The GEO record and file are currently private and are expected to become public when the paper is published. Until release, access requires the appropriate GEO authorization. No reviewer credentials are included in this repository.

Place the downloaded file beside `reproduce_figures_tables.R`, retaining the exact filename. **Do not decompress it:** R's `load()` reads the gzip-compressed workspace directly. The contents are the same saved analysis objects used by the original pipeline.

After the GEO file is publicly available, it can also be downloaded from bash:

```bash
curl --fail --location \
  'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE345813&format=file&file=GSE345813%5FGeoMx%5FAnalysis%5FComplete%5FWorkspace%2ERData%2Egz' \
  --output GSE345813_GeoMx_Analysis_Complete_Workspace.RData.gz
```

The workspace is hosted on GEO rather than GitHub; Git LFS is not required.

## Run the pipeline

From the repository directory:

```bash
Rscript reproduce_figures_tables.R > reproduce_figures_tables.log 2>&1
```

Outputs are written to `manuscript_10092026_outputs/`, created automatically. To watch progress from another terminal:

```bash
tail -f reproduce_figures_tables.log
```

Optional paths can be supplied through environment variables:

```bash
SLICE_WORKSPACE="/path/to/GSE345813_GeoMx_Analysis_Complete_Workspace.RData.gz" \
SLICE_OUTPUT_DIR="/path/to/results" \
Rscript reproduce_figures_tables.R > reproduce_figures_tables.log 2>&1
```

Use an empty output directory when switching script versions so obsolete files from previous numbering are not mistaken for current outputs. The script overwrites matching output filenames. It sets a random seed for plot jitter; identical rendering across different R/package versions is not guaranteed.

## Generated figures

| Manuscript item | Content |
| --- | --- |
| Figure 3 | Media-associated differential expression volcano plots, including comparisons with T0 and between media |
| Figure 4A/B | Hallmark ssGSEA score changes in glomerular and tubulointerstitial ROIs |
| Figure 5A/B | Temporal gene-expression heatmaps over 72 hours |
| Figure 6A | APOL1 forest plot across five donors |
| Figure 6B | Hallmark interferon gamma response scores across treatment conditions in MTS008 and MTS029 |
| Figure 7A | Hallmark fGSEA: IFNγ versus control |
| Figure 7B | Hallmark fGSEA: baricitinib versus DMSO |
| Supplemental Figure 2 | Expression distributions before and after quantile normalization |
| Supplemental Figure 3 | PCA before and after quantile normalization |
| Supplemental Figure 5A | Glomerular APOL1 expression across treatments |
| Supplemental Figure 6 | Hallmark fGSEA: IFNγ versus IFNγ plus baricitinib |
| Supplemental Figure 7 | IFNγ-associated pathway-score changes in glomeruli and tubulointerstitium |

Figures are exported as PDFs. Some manuscript figures are exported as separate panels or compartment-specific files and require assembly for publication.

## Generated supplemental tables

Excel workbooks use the following numbering:

| Table | Output filename | Content |
| --- | --- | --- |
| 3 | `Suppl_Table_03_GeoMx_ROI_Summary.xlsx` | ROI metadata summaries and totals by donor |
| 4 | `Suppl_Table_04_Media_Effects_DEG_Results.xlsx` | Media-associated differential expression results |
| 5 | `Suppl_Table_05_Media_Effects_ssGSEA_Pathway_Stats.xlsx` | Media-associated Hallmark pathway statistics |
| 6 | `Suppl_Table_06_Temporal_Gene_Clusters_72hr.xlsx` | Temporal gene patterns and expression summaries |
| 7 | `Suppl_Table_07_IFN_vs_Control_MetaAnalysis_DEGs.xlsx` | IFNγ meta-analysis and per-donor effect estimates |
| 8 | `Suppl_Table_08_Baricitinib_vs_DMSO_MetaAnalysis_DEGs.xlsx` | Baricitinib meta-analysis, concordance, and per-donor estimates |
| 9 | `Suppl_Table_09_fGSEA_Results_All_Comparisons.xlsx` | fGSEA results for IFNγ/control, baricitinib/DMSO, and IFNγ/IFNγ plus baricitinib |
| 10 | `Suppl_Table_10_Glom_vs_TI_IFN_Pathway_Comparison.xlsx` | Compartment-specific pathway changes and treatment-by-compartment interaction FDR |

The large quantile-normalized expression workbook is no longer a numbered supplemental table and is not exported by this script.

## Analyses

Glomerular and tubulointerstitial expression data are quantile-normalized separately. Media-associated differential expression uses limma with donor adjustment. The IFNγ response combines within-donor estimates across five donors; baricitinib comparisons involve two donors. Hallmark pathway analyses use precomputed ssGSEA scores and saved preranked fGSEA results. The compartment comparison uses a treatment-by-compartment interaction model for donor MTS029, followed by Benjamini–Hochberg correction.

The temporal analyses use a single donor, MTS028. ROIs are multicellular regions, and multiple ROIs from one donor are not independent donor replicates. See the manuscript for experimental design, thresholds, and interpretation.

## Separate source material

The following items require source images or data outside the GeoMx workspace:

- Figures 1 and 2: workflow, lightsheet microscopy, LDH assays, and histology.
- Figure 6C: APOL1 RNA in situ hybridization.
- Figure 7C: mechanism artwork, `Baricitinib JAK-STAT-OxPhos.pdf`.
- Supplemental Figures 1 and 4: histology and NPHS1/SLC6A19 immunohistochemistry.
- Supplemental Figure 5B/C: APOL1 qPCR and immunohistochemistry.
- Supplemental Tables 1 and 2: donor information and primer sequences.
- Supplemental Movies 1 and 2: lightsheet microscopy movies.

The script writes `Manuscript_Output_Checklist.csv` identifying these separate source requirements. It does not generate these assets or claim to reproduce complete publication composites.

## Data availability

Raw and processed spatial transcriptomic data are associated with GEO accession **[GSE345813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345813)**. Access depends on the dataset's release status.

An interactive data explorer is available at:
[Human kidney slice GeoMx Explorer](https://faculty.washington.edu/shreeram/GeoMx_Explorer/).

The saved R workspace is a distinct input required by this script; use the direct workspace link above after GEO release.

## Citation and contact

Please cite the associated manuscript when using these analyses. A journal citation and DOI can be added after publication.

Contact: Shreeram Akilesh, [shreeram@uw.edu](mailto:shreeram@uw.edu).

## License

Code is distributed under the MIT license in `LICENSE`. Dataset and third-party artwork reuse remains subject to the applicable source terms. BioRender artwork is not relicensed by the repository's MIT code license.
