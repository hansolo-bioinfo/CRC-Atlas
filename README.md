# CRCAtlas

CRC Atlas is an interactive web server that is built on R Shiny for researchers to investigate landscape of gene-coexpressing modules (GEM) in curated single cell RNAseq (scRNAseq) data of colorectal cancer. In CRC Atlas, we have collected 279 samples and 626858 cells from 8 scRNAseq datasets.

Citation: (biorxiv link if applicable)

<img width="1897" alt="Screenshot 2024-01-21 at 16 34 30" src="https://github.com/hansolo-bioinfo/CRCAtlas/assets/65295899/b4483fc2-d893-4bdf-81e4-899683a1aef0">

## Installation

```R
devtools::install_github("hansolo-bioinfo/CRCAtlas")
library(CRCAtlas)
CRCAtlasApp()
```

## How to Use

CRC Atlas offers four functions for researchers to look into the distribution of GEMs in the colorectal cancer

### 1) GEM Information

Look into the distribution of each gene-coexpressing modules (GEM) and associate top weighted genes by selecting GEM in below dropdown menu. Left panel shows the UMAP, top 20 genes and distribution of the inquired GEM over cell subtypes. Right panel displays the dot plots of all useful GEMs of this major cell type.

<img width="1880" alt="Screenshot 2024-01-21 at 16 34 48" src="https://github.com/hansolo-bioinfo/CRCAtlas/assets/65295899/5e225c41-0c4e-4e08-b248-f392fdf6c075">

### 2) GEM-LR Correlation

Look into the nonlinear Sigmoid relationship between gene-coexpressing modules (GEM) and ligand -> receptor (LR) pair in our curated CRC Atlas (n=151). On the right panel, you can also find the Top 10 significant GEM-LR pairs for your search, ranked by RSE (residual standard error).

For example, the top correlated LRs with TNK26:

<img width="1903" alt="Screenshot 2024-01-21 at 16 35 04" src="https://github.com/hansolo-bioinfo/CRCAtlas/assets/65295899/14f743da-921a-4dcc-b57a-a205dfc6a5fb">

Or the top correlated GEMs with the pair CD274 -> PDCD1:

<img width="1901" alt="Screenshot 2024-01-21 at 16 35 21" src="https://github.com/hansolo-bioinfo/CRCAtlas/assets/65295899/e71d9fe1-adaa-42c9-bfd8-b650cddaa9ef">


### 3) GEM-LR Causality

Conditional independence is an important concept in causal analysis, and in our setting, if a pair of highly correlated GEMs becomes independent when conditioning on a ligand -> receptor (LR) pair, the LR is likely involved in transmitting signal between cells expressing the GEMs. Finding a GEM X-LR-GEM Y triplet leads to the hypothesis depicted in below figure. Herein, we provide the triplet results that pass the conditional independence test. (choose at least one)

<img width="1888" alt="Screenshot 2024-01-21 at 16 36 13" src="https://github.com/hansolo-bioinfo/CRCAtlas/assets/65295899/66dd045c-3942-459f-b80e-9c70adf90697">

### 4) Deconvolution

To deconvolute the bulk-RNAseq expression data using our curated GEM gene list, you will need to upload the gene expression in <.csv> format of which the rowname is gene symbol and colname is sample ID. The data is better deconvoluted if the input matrix is log2(TPM+1) or log2(FKPM+1). 

<img width="1909" alt="Screenshot 2024-01-21 at 16 36 31" src="https://github.com/hansolo-bioinfo/CRCAtlas/assets/65295899/ceacc4b0-e944-4eb0-984f-7e2d5920d60d">
