# BCRi

Tools for exploring the interplay between B cell receptor (BCR) sequence diversity and cellular phenotypes at single-cell resolution.

## Overview

The BCRi toolset incorporates adaptive immune receptor repertoire (AIRR) data at single-cell resolution that has been annotated for cell state or any other phenotypic covariate of interest. BCRi then constructs the repertoire-phenotype joint distribution in which each row corresponds to a distinct BCR sequence—by default, the sequence encoding the VDJ heavy chain—and each column represents a cell phenotype. A symmetric similarity matrix `k` defines the relationship between all pairs of BCR sequences; in the simplest case, `k` is the identity matrix spanning the `n` unique sequences.

BCRi leverages an information theoretic framework applied to these data structures to quantify relationships between cell phenotype and immune repertoire. The package formalizes repertoire diversity metrics based on Shannon entropy that capture both varying abundances of BCR clonotypes and clonal families, as well as the varying sequence similarity at single-cell resolution.

## Key Features

- Construct repertoire–phenotype joint distributions directly from annotated single-cell AIRR datasets.
- Quantify diversity with Shannon entropy-based metrics that incorporate sequence similarity.
- Compute similarity-sensitive diversity profiles, meta-diversity summaries, and bootstrap-supported uncertainty estimates.
- Integrate seamlessly with Seurat metadata through helper utilities like `functional_diversity_from_seurat()`.
- Extend analyses with customizable affinity matrices and phenotype groupings.

## Installation

Install the development version from GitHub with `remotes` (or `devtools`):

```r
install.packages("remotes")
remotes::install_github("DoaneAS/BCRi")
```

Alternatively, clone the repository and install from source:

```r
git clone https://github.com/DoaneAS/BCRi.git
setwd("BCRi")
devtools::install()
```

## Getting Started

```r
library(BCRi)

# Example: run diversity analysis on a prepared AIRR data.frame
diversity_metrics <- functional_diversity(
  db = airr_table,
  groupID = "clone_001",
  phenotype_var = "subset",
  group = "subject_id"
)

# Example: compute diversity directly from a Seurat object
diversity_from_seurat <- functional_diversity_from_seurat(
  seurat_obj = seurat_object,
  groupID = "subject_A",
  phenotype_var = "celltype",
  group = "subject_id"
)
```

## Documentation

- Function reference: see the `man/` directory or run `help(package = "BCRi")` after installation.
- Tutorials and vignettes: forthcoming. Contributions and suggestions are welcome—feel free to open an issue or submit a pull request.

## Contributing

1. Fork the repository and create a feature branch.
2. Make your changes, adding tests and documentation when appropriate.
3. Submit a pull request describing the motivation and approach.

Please follow the existing code style and include reproducible examples when reporting issues.

## License

Specify your license here (e.g., MIT, GPL-3). Update this section once a license has been chosen.


