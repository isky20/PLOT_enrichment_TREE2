# Dendrogram Plotting Function

## Overview

The `create_dendrogram_plots` function processes biological data from a CSV file, normalizes subsets based on categories, generates dendrogram plots, and saves them as SVG files. This function is useful for visualizing relationships between biological samples based on gene/protein presence.

## Function Definition

```r
create_dendrogram_plots("ATTR_AL_FULL_filter_adj.csv", "ATTR_AL_FULL_filter",
                        width = 200, offset= 6, xlim_max= 25, height = 300)
