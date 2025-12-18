[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/k4IGslji)
[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-2e0aaae1b6195c2367325f4f02e2d04e9abb55f0b24a779b69b11b9e10269abc.svg)](https://classroom.github.com/online_ide?assignment_repo_id=21939418&assignment_repo_type=AssignmentRepo)

# Reproducing Shiao et al. 2015: Expression Divergence of Chemosensory Genes in Drosophila sechellia

## Paper Reference
**Expression Divergence of Chemosensory Genes between *Drosophila sechellia* and Its Sibling Species and Its Implications for Host Shift**  
Shiao et al., *Genome Biology and Evolution*, 2015

### Members
* 劉愷崴, 111752011
* 周彥廷, 110601043
* 曹柏泓, 113753108
* 柏慶銘, 110703028

---

## Demo

### Reproduce Table 1 (Differential Expression of Chemosensory Genes)
```bash
# Usage: Rscript code/reproduce_table1.R [base_dir] [output_dir]

# With default paths:
Rscript code/reproduce_table1.R

# With custom input and output directories:
Rscript code/reproduce_table1.R "/path/to/data" "/path/to/output"
```

---

## Folder Organization

Idea by Noble WS (2009) [A Quick Guide to Organizing Computational Biology Projects.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424) PLoS Comput Biol 5(7): e1000424.

### docs
* Presentation slides: `1131_bioinformatics_FP_group2.pptx`

### data

| Dataset | Source | Format | Description |
|---------|--------|--------|-------------|
| GSE67587_RAW | GEO | FPKM tracking files | *D. sechellia* Taiwan samples (6 files) |
| GSE67861_RAW | GEO | FPKM tracking files | *D. sechellia* Japan samples (6 files) |
| GSE67862_RAW | GEO | FPKM tracking files | *D. simulans* Japan samples (6 files) |
| GSE99545_RAW | GEO | FPKM tracking files | *D. melanogaster* samples (2 files) |
| all_3spe-exist+noDup_paper1.csv | Ortholog mapping | CSV | Orthologous genes across 3 Drosophila species (~10,035 genes) |

### code

#### Packages Used
**Original packages from the paper:**
* `NOISeq` - Differential expression analysis with upper quartile normalization

**Additional packages:**
* R: `dplyr`, `tidyr`

#### Analysis Script
| Script | Description |
|--------|-------------|
| `reproduce_table1.R` | Reproduces Table 1 using NOISeq differential expression. Accepts optional `base_dir` (input) and `output_dir` arguments. |

#### Analysis Steps
1. Load ortholog table and map gene IDs across species
2. Build expression matrix from FPKM tracking files
3. Apply Upper Quartile Normalization (NOISeq `uqua()`)
4. Run differential expression analysis:
   - `noiseq(replicates="no")` for D. melanogaster comparisons (1 replicate)
   - `noiseq()` for D. sechellia vs D. simulans (3 replicates each)
5. Filter DEGs at q ≥ 0.95 threshold
6. Generate Table 1

### results

| Output | Description |
|--------|-------------|
| `expression_matrix_raw.csv` | Raw FPKM expression matrix |
| `expression_matrix_normalized.csv` | Upper quartile normalized expression matrix |
| `DEG_table_Dsec_M_TW_vs_Dmel_M_TW.csv` | DEG analysis: *D. sechellia* vs *D. melanogaster* (Male, Taiwan) |
| `DEG_table_Dsec_F_TW_vs_Dmel_F_TW.csv` | DEG analysis: *D. sechellia* vs *D. melanogaster* (Female, Taiwan) |
| `DEG_table_Dsec_M_JP_vs_Dsim_M_JP.csv` | DEG analysis: *D. sechellia* vs *D. simulans* (Male, Japan) |
| `DEG_table_Dsec_F_JP_vs_Dsim_F_JP.csv` | DEG analysis: *D. sechellia* vs *D. simulans* (Female, Japan) |
| `Table1_reproduced_R.csv` | Reproduced Table 1: FPKM and log2-ratios of 18 chemosensory DEGs |

---

## What We Reproduced

### Table 1
Expression levels (FPKMs) and Log2-ratios of 18 chemosensory genes differentially expressed in *D. sechellia* compared with *D. melanogaster* or *D. simulans*:

**Upregulated (13 genes):** Obp19a, Obp50a, Obp56d, CheA75a, CheA87a, Or23a, Or35a, Or56a, Or67b, Or85b, Or85c, Gr64f, Ir84a

**Downregulated (5 genes):** Obp83a, Obp99c, Obp99d, Or9a, Or42b

### Fig. 1(a) 
Venn diagram showing differentially expressed genes (DEG) in D. sechellia and other two species. `plot_venn_diagram.R` takes the 4 output files starting with `DEG` (`DEG_table_Dsec_M_TW_vs_Dmel_M_TW`, `DEG_table_Dsec_M_JP_vs_Dsim_M_JP`, `DEG_table_Dsec_F_TW_vs_Dmel_F_TW` and `DEG_table_Dsec_F_JP_vs_Dsim_F_JP`) and exports two Vwnn diagrams: `up.jpeg` (up-regulated DEGs) and `down.jpeg` (down-regulated DEGs).

### Fig. 2(a)
Heat maps of expression profiles of Or genes.
`heat_map.Rmd` takes the file `fpkm.csv` and recreates Fig.2(a) from the original article.
---

## References

### Paper
* Shiao MS, et al. (2015) Expression Divergence of Chemosensory Genes between *Drosophila sechellia* and Its Sibling Species and Its Implications for Host Shift. *Genome Biology and Evolution*, 7(6):1563-1576.

### Packages
* Tarazona S, et al. (2015) NOISeq: Exploratory analysis and differential expression for RNA-seq data. *Bioinformatics*, 31(5):1873-1875.

### Data
* GEO Accession: GSE67587, GSE67861, GSE67862, GSE99545
