# hhmR
## Hierarchical Heatmaps

This package allows users to create high-quality heatmaps from labelled hierarchical data. Specifically, it requires input data in the form of a two-level hierarchical structure. It will produce a heatmap where each row and column represent a *category* at the lower level. Rows and columns are then grouped into a higher-level *groupings*. Names for each higher-level *category* and *groupings* are shown in the margins. While other packages (e.g. `dendextend`) allow heatmap rows and columns to be arranged by groups only, `hhmR` also allows the labelling of the data at both the category and group level.

The package has two main functions: `hhm` and `tshhm`. These are explained in more detail below.

<br>

### `hhm` (Hierarchical Heatmap)

`hhm` creates a labelled heatmap from hierarchical data. This function is useful to create a heatmap where the categories shown on both the x and y axis are grouped in some way. This heatmap will order the categories by their assigned groupings, and present both the categories and grouping labels along the axes. An example using geographic data might be a series of smaller geographies (lower categories) which aggregate into larger geographical regions (upper groups).

This function requires a data.frame containing columns which specify the lower categories (`ylower`, `xlower`) and upper groupings (`yupper`, `xupper`). These categories and groupings are used to arrange and label rows and columns on a heatmap. The data.frame must contain a `values` variable containing values to populate the heatmap. Note that the groupings will, by default, be arranged alphabetically (top to bottom / left to right). The ordering of the groups can be manually specified by converting `yupper` and/or `xupper` to factors. In this case, the groupings will be ordered based on the ordering of the factor levels provided.

Below is an example of the `hhm` function's application. For a more in-depth description of it's usage, see the [package vignette](https://sgmmahon.github.io/hhmR/articles/hhmR_overview.html).

```
# Import package
library(hhmR)

# Import toy demonstration dataset (see `?example_migration` for details)
data(example_migration)

# Intial heatmap
hierarchical_heatmap = hhm(df = example_migration,
                           ylower = "Origin County",
                           xlower = "Destination County",
                           yupper = "Origin Region",
                           xupper = "Destination Region",
                           values = "Migration",
                           yttl_width = 0.22,
                           xttl_height = 0.4)

# View result
hierarchical_heatmap
```
![ ](https://github.com/sgmmahon/hhmR/blob/main/docs/images/hhm_initial_output.png?raw=true)

<br>

### `tshhm` (Time-Series Hierarchical Heatmap)

`tshhm` creates a labelled time-series heatmap from heirarchical data. This function is useful if you wish to create a time-series heatmap where the categories shown on the y-axis can be grouped in some way. This heatmap will order the categories by their assigned group and present both the categories and group labels along the y-axis. An example might be a series of smaller geographies (lower categories) which aggregate into larger geographical regions (upper groups).

This function requires a data.frame containing columns that specify the lower categories (`lower`) and upper groups (`upper`) that each value corresponds to. These categories and groups will be used to arrange and label the rows of the heatmap. The data.frame must also contain a `values` variable, containing the values used to populate the heatmap, and a `times` variable, containing the time period during which each value was observed. Note that the groups in `upper` will by default be arranged alphabetically (top to bottom). The ordering of the groups can be manually specified by converting `upper` to a factor. In this case, the groups will be ordered based on the ordering of the factor levels. The ordering of rows within each group can also be specified using the `sort_lower` variable.

Below is an of how the `tshhm` function can be used. For more details, see the [package vignette](https://sgmmahon.github.io/hhmR/articles/hhmR_overview.html).
```
# Import packages
library(dplyr)

# Import toy demonstration dataset (see `?example_time_series` for details)
data(example_time_series)

# Intial heatmap
time_series_heatmap = tshhm(df = example_time_series,
                            lower  = "County",
                            upper  = "Region",
                            times  = "Year",
                            values = "Immigration",
                            yttl_width  = 0.25)

# View result
time_series_heatmap
```
![ ](https://github.com/sgmmahon/hhmR/blob/main/docs/images/tshhm_initial_output.png?raw=true)
