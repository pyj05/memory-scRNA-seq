# Script consolidation map

| Workflow script | Consolidated source scripts | Scope |
| --- | --- | --- |
| `scripts/01_composition_tsne.R` | `tsne_proportion_plot.R`, `tsne_proportion_subclass.R` | Cell-class and neuronal-subclass composition and t-SNE plots |
| `scripts/02_neuronal_deg_enrichment.R` | `go_kegg_reactome.R`, `filter_enrichment.R`, `volcano_DEG.R`, `retrieval_vs_acquisition_enrichment.R`, `retrieval_vs_acquisition_volcano.R`, `retrieval_vs_acquisition_enrichment_plot.R`, `enrichment_plot.R`, `enrichment_plot_nolabel.R`, `heatmap.R` | All-neuronal differential expression, enrichment, selected enrichment tables, and visualizations |
| `scripts/03_ppi_analysis.R` | `PPI.R`, `PPI_betweenness.R`, `PPI_hub.R` | STRING PPI networks, betweenness plots, and hub-gene tables |
| `scripts/04_neuronchat_analysis.R` | `neuronchat.R` | NeuronChat interaction analysis |
| `scripts/05_C0_deg_enrichment.R` | `C0_go_kegg_reactome.R`, `C0_volcano_DEG.R`, `C0_retrieval_vs_acquisition_volcano.R`, `C0_enrichment_plot.R` | C0-only differential expression, enrichment, and visualizations |

The merged scripts retain the original analysis order and parameters. They add only explicit Snakemake inputs, workflow-local output paths, differential-expression table export for the PPI dependency, and the C0 object reference required by the C0 retrieval-versus-acquisition comparison.
