# hhmR
## Hierarchical Heatmaps

Allows users to create high-quality heatmaps from labelled hierarchical data. Specifically, it requires input data in the form of a two-level hierarchical structure. It will produce a heatmap where each row and column represent a *category* at the lower level. Rows and columns are then grouped into a higher-level *groupings*. Names for each higher-level *category* and *groupings* are shown in the margins. While other packages (e.g. `dendextend`) allow heatmap rows and columns to be arranged by groups only, `hhmR` also allows the labelling of the data at both levels categories and groupings.

The package has two main functions: `hhm` and `tshhm`. These are explained in more detail below.

### `hhm` (Hierarchical Heatmap)

`hhm` creates a labelled heatmap from hierarchical data. This function is useful to create a heatmap where the categories shown on both the x and y axis are grouped in some way. This heatmap will order the categories by their assigned groupings, and present both the categories and grouping labels along the axes. An example using geographic data might be a series of smaller geographies (lower categories) which aggregate into larger geographical regions (upper groups).

This function requires a data.frame containing columns which specify the lower categories (`ylower`, `xlower`) and upper groupings (`yupper`, `xupper`). These categories and groupings are used to arrange and label rows and columns on a heatmap. The data.frame must contain a `values` variable containing values to populate the heatmap. Note that the groupings will, by default, be arranged alphabetically (top to bottom / left to right). The ordering of the groups can be manually specified by converting `yupper` and/or `xupper` to factors. In this case, the groupings will be ordered based on the ordering of the factor levels provided.

Below are some examples of how the `hhm` function can be used.
```
# Import package
library(hhm)

# Import toy demonstration dataset (see `?example_migration` for see details)
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

# Remove diagonal from heatmap (i.e. hide static populations)
removed_diag         = hhm(df = example_migration,
                           ylower = "Origin County",
                           xlower = "Destination County",
                           yupper = "Origin Region",
                           xupper = "Destination Region",
                           values = "Migration",
                           yttl_width = 0.22,
                           xttl_height = 0.4,
                           rm_diag = TRUE)

# View result
removed_diag

# Nomalise the legend
normalised_lgd       = hhm(df = example_migration,
                           ylower = "Origin County",
                           xlower = "Destination County",
                           yupper = "Origin Region",
                           xupper = "Destination Region",
                           values = "Migration",
                           yttl_width = 0.22,
                           xttl_height = 0.4,
                           rm_diag = TRUE,
                           norm_lgd = TRUE)

# View result
normalised_lgd

# Manually define colour scheme for heatmap (uses viridis colour scheme)
viridis_12 = c("#440154FF","#482173FF","#433E85FF","#38598CFF","#2D708EFF","#25858EFF",
               "#1E9B8AFF","#2BB07FFF","#51C56AFF","#85D54AFF","#C2DF23FF","#FDE725FF")

# Assign continuous colour scheme
cont_clrs            = hhm(df = example_migration,
                           ylower = "Origin County",
                           xlower = "Destination County",
                           yupper = "Origin Region",
                           xupper = "Destination Region",
                           values = "Migration",
                           yttl_width = 0.22,
                           xttl_height = 0.4,
                           rm_diag = TRUE,
                           norm_lgd = TRUE,
                           cclrs = viridis_12)

# View result
cont_clrs

# Break legends into a specified number of bins
# (of equal intervals between 0 and the maximum value in `values`)
bins_15              = hhm(df = example_migration,
                           ylower = "Origin County",
                           xlower = "Destination County",
                           yupper = "Origin Region",
                           xupper = "Destination Region",
                           values = "Migration",
                           yttl_width = 0.22,
                           xttl_height = 0.4,
                           rm_diag = TRUE,
                           bins = 15)

# View result
bins_15

# Manually break data into categories using user-specified intervals.
# In this instance, the `hhmR` function `log_seq` has been used to create a
# vector of logarithmicly increasing values between 1 and the maximum value
# in the dataset not on the diagonal.
cbrks = log_seq(example_migration[example_migration[["Origin County"     ]] !=
                                  example_migration[["Destination County"]],] %>%
                .$Migration %>% max(), 12, rmv_extremes = TRUE)

# Manually assign legend categories
legend_cats          = hhm(df = example_migration,
                           ylower = "Origin County",
                           xlower = "Destination County",
                           yupper = "Origin Region",
                           xupper = "Destination Region",
                           values = "Migration",
                           yttl_width = 0.22,
                           xttl_height = 0.4,
                           rm_diag = TRUE,
                           cbrks = cbrks)

# View result
legend_cats

# Manually assign colours to legend categories
cat_clrs             = hhm(df = example_migration,
                           ylower = "Origin County",
                           xlower = "Destination County",
                           yupper = "Origin Region",
                           xupper = "Destination Region",
                           values = "Migration",
                           yttl_width = 0.22,
                           xttl_height = 0.4,
                           rm_diag = TRUE,
                           cbrks = cbrks,
                           cclrs = viridis_12)

# View result
cat_clrs
```

### `tshhm` (Time-Series Hierarchical Heatmap)

`tshhm` creates a labelled time-series heatmap from heirarchical data. This function is useful if you wish to create a time-series heatmap where the categories shown on the y-axis can be grouped in some way. This heatmap will order the categories by their assigned group and present both the categories and group labels along the y-axis. An example might be a series of smaller geographies (lower categories) which aggregate into larger geographical regions (upper groups).

This function requires a data.frame containing columns that specify the lower categories (`lower`) and upper groups (`upper`) that each value corresponds to. These categories and groups will be used to arrange and label the rows of the heatmap. The data.frame must also contain a `values` variable, containing the values used to populate the heatmap, and a `times` variable, containing the time period during which each value was observed. Note that the groups in `upper` will by default be arranged alphabetically (top to bottom). The ordering of the groups can be manually specified by converting `upper` to a factor. In this case, the groups will be ordered based on the ordering of the factor levels. The ordering of rows within each group can also be specified using the `sort_lower` variable.

Below are some examples of how the `tshhm` function can be used.
```
# Import packages
library(dplyr)

# Import toy demonstration dataset (see `?example_time_series` for see details)
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

# Arrange counties within each region by total number of immigrants
# across all five years (ascending from top to bottom)
sort_ascending      = tshhm(df = example_time_series,
                            lower  = "County",
                            upper  = "Region",
                            times  = "Year",
                            values = "Immigration",
                            sort_lower = "sum_ascend",
                            yttl_width  = 0.25)

# View result
sort_ascending

# Nomalise the legend
normalised_lgd      = tshhm(df = example_time_series,
                            lower  = "County",
                            upper  = "Region",
                            times  = "Year",
                            values = "Immigration",
                            sort_lower = "sum_ascend",
                            norm_lgd = TRUE,
                            yttl_width  = 0.25)

# View result
normalised_lgd

# Manually define colour scheme for heatmap (uses viridis colour scheme)
viridis_12 = c("#440154FF","#482173FF","#433E85FF","#38598CFF","#2D708EFF","#25858EFF",
               "#1E9B8AFF","#2BB07FFF","#51C56AFF","#85D54AFF","#C2DF23FF","#FDE725FF")

# Assign continuous colour scheme
cont_clrs           = tshhm(df = example_time_series,
                            lower  = "County",
                            upper  = "Region",
                            times  = "Year",
                            values = "Immigration",
                            sort_lower = "sum_ascend",
                            norm_lgd = TRUE,
                            cclrs = viridis_12,
                            yttl_width  = 0.25)

# View result
cont_clrs

# Assign colour for NA values
na_clrs             = tshhm(df = example_time_series,
                            lower  = "County",
                            upper  = "Region",
                            times  = "Year",
                            values = "Immigration",
                            sort_lower = "sum_ascend",
                            norm_lgd = TRUE,
                            cclrs = viridis_12,
                            na_colour = "grey80",
                            yttl_width  = 0.25)

# View result
na_clrs

# Break legends into a specified number of bins
# (of equal intervals between 0 and the maximum value in `values`)
bins_15             = tshhm(df = example_time_series,
                            lower  = "County",
                            upper  = "Region",
                            times  = "Year",
                            values = "Immigration",
                            sort_lower = "sum_ascend",
                            bins = 15,
                            yttl_width  = 0.25)

# View result
bins_15

# Manually break data into categories using user-specified intervals.
# In this instance, the `hhmR` function `log_seq` has been used to create a
# vector of logarithmicly increasing values between 1 and the maximum value
# in the dataset.
cbrks = log_seq(example_time_series %>% .$Immigration %>% max(na.rm = TRUE),
                12, rmv_extremes = TRUE)

# Manually assign legend categories
legend_cats         = tshhm(df = example_time_series,
                            lower  = "County",
                            upper  = "Region",
                            times  = "Year",
                            values = "Immigration",
                            sort_lower = "sum_ascend",
                            cbrks = cbrks,
                            yttl_width  = 0.25)

# View result
legend_cats

# Manually assign colours to legend categories
cat_clrs            = tshhm(df = example_time_series,
                            lower  = "County",
                            upper  = "Region",
                            times  = "Year",
                            values = "Immigration",
                            sort_lower = "sum_ascend",
                            cbrks = cbrks,
                            cclrs = viridis_12,
                            na_colour = "grey80",
                            yttl_width  = 0.25)

# View result
cat_clrs

# Manually define order of x-axis and groups using factor levels
new_time_series = example_time_series %>%
                  mutate(Year   = factor(Year,
                                         levels = c(2012,2011,2014,
                                                    2013,2015)),
                         Region = factor(Region,
                                         levels = c("North","Midlands",
                                                    "South West",
                                                    "South East")))

# Manually define order of x-axis and groups
rearrange_axes      = tshhm(df = new_time_series,
                            lower  = "County",
                            upper  = "Region",
                            times  = "Year",
                            values = "Immigration",
                            sort_lower = "sum_ascend",
                            cbrks = cbrks,
                            cclrs = viridis_12,
                            na_colour = "grey80",
                            yttl_width  = 0.25)

# View result
rearrange_axes
```
