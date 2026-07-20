# Introduction

This workflow organizes the supplied annotated Seurat objects and downstream R analyses into reproducible Snakemake rules. This is part of the data analysis for the manuscript "UTX preserves remote memory by orchestrating engram persistence through chromatin remodeling."

## Inputs

- `data/annotation_data.rds`: annotated all-cell Seurat object.
- `data/neurons_subclass.rds`: neuronal Seurat object with subclass labels.

The data files are hard-linked to the supplied files to avoid duplicating large RDS files. Replace the files in `data/` with copies when the workflow is moved to a different volume or computer.

## Analyses

- Cell-class and neuronal-subclass composition with t-SNE visualizations.
- Differential expression, functional enrichment, volcano plots, and heatmap for all neuronal cells.
- STRING PPI networks, betweenness visualizations, and hub-gene tables.
- NeuronChat comparison across the four experimental groups.
- Differential expression and functional enrichment restricted to subclass C0.

## Run

Create the Snakemake environment, then run from this directory:

```
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
snakemake --use-conda --cores 1
```

`neuronchat` uses `CellChat` and `NeuronChat` in addition to the Conda environment. Install both R packages in the environment before enabling the `neuronchat` rule, following their package installation instructions. All other rules can be run independently, for example `snakemake --use-conda --cores 1 results/neurons`.

## Results

Each rule writes to its own directory below `results/`. The neuronal differential-expression rule also writes tab-separated tables that are used as the explicit input to the PPI rule.
