utils::globalVariables(".")
#' @title Time-series Hierarchical Heatmap
#'
#' @description Creates a labelled time-series heatmap from heirarchical data. This
#' function is useful if you wish to create a time-series heatmap where the
#' categories shown on the y axis can be grouped in some way. This heatmap
#' will order the categories by their assigned group and present both the categories
#' and group labels along the y-axis. An example might be series of smaller
#' geographies (lower categories) which aggregate into larger geographical regions
#' (upper groups).
#'
#' @param df A data.frame with containing values with which to populate the heatmap.
#' The data.frame must include columns specifying the lower categories (`lower`) and
#' upper groups (`upper`) that each value corresponds to. These categories and
#' groups will be used to arrange and label the rows of the heatmap. `df` must also
#' contain a `values` variable, containing the values used to populate the heatmap,
#' and a `times` variable, containing the time period during which each value was
#' observed. Note that the groups in `upper` will by default be arranged
#' alphabetically (top to bottom). The ordering of the groups can be manually
#' specified by converting `upper` to a factor. In this case, the groups
#' will be ordered based on the ordering of the factor levels. The ordering of rows
#' within each group can also be specified using the `sort_lower` variable.
#' @param lower A column in `df` containing the categories that will be presented
#' as rows along the y-axis of the heatmap.
#' @param upper A column in `df` containing the groupings that will be used to
#' arrange the heatmap rows.
#' @param times A column in `df` containing the time-period during which each
#' each value in `values` was observed.
#' @param values A column in `df` containing the values used to populate the
#' heatmap.
#' @param sort_lower Option to define how rows (lower) within each group (upper)
#' are ordered. The default options is `alphabetical`, which orders rows in
#' alphabetical order from top to bottom. Other options include `sum_ascend` and
#' `mean_ascend`, which order rows in ascending order (top to bottom) based on
#' the row totals and row means respectively. This order can be reversed with the
#' options `sum_descend` and `mean_descend`.
#' @param lgttl Option to manually define legend title.
#' @param bins Option to break the data into a specified number of groups
#' (defaults to `NULL`). The thresholds between these groups will be equally
#' spaced between the minimum and maximum values observed in `values`.
#' @param cbrks Vector of custom breaks, if users wish to use a discrete legend
#' colour scheme (defaults to `NULL`). For example, a supplied vector of `c(5,10,
#' 20)` would break he values up into 5 ordered groups of ranges 0, 0-5, 5-10,
#' 10-20 and 20+.
#' @param cclrs Vector of hexcodes, which to create a custom legend colour scheme
#' (defaults to `NULL`). If `cbrks` is supplied, `cclrs` must have a length
#' two longer than `cbrks`. If `bins` is supplied, `cclrs` must have a length
#' equal to the values provided to `bins`.
#' @param norm_lgd Normalised to between 0 and 1 in legend (defaults to `FALSE`).
#' Allows for consistency when comparing heatmaps across different datasets. At
#' present, this only works all heatmap values are positive.
#' @param lgdps If using custom breaks, define the number of decimal points to
#' round the legend scale to (defaults to 0). If `norm_lgd` is `TRUE`, it will
#' default to 3.
#' @param na_colour Option to define the colour of NA values in the legend (defaults
#' to `NULL`, meaning NA values will be assigned no colour).
#' @param xttl_height The space allocated to the title on the x-axis as a
#' proportion of the heatmap's height (defaults to 0.05).
#' @param yttl_width The space allocated to the group titles on the y-axis as a
#' proportion of the heatmap's width (defaults to 0.15).
#'
#' @return A ggplot object containing the final heatmap.
#'
#' @examples
#' library(dplyr)
#'
#' # Import toy demonstration dataset (see `?example_time_series` for see details)
#' data(example_time_series)
#'
#' # Intial heatmap
#' time_series_heatmap = tshm(df = example_time_series,
#'                            lower  = "County",
#'                            upper  = "Region",
#'                            times  = "Year",
#'                            values = "Immigration")
#'
#' # View result
#' time_series_heatmap
#'
#' # Arrange counties within each region by total number of immigrants
#' # across all five years (ascending from top to bottom)
#' sort_ascending      = tshm(df = example_time_series,
#'                            lower  = "County",
#'                            upper  = "Region",
#'                            times  = "Year",
#'                            values = "Immigration",
#'                            sort_lower = "sum_ascend")
#'
#' # View result
#' sort_ascending
#'
#' # Increase spacing between plots
#' increase_spaces     = tshm(df = example_time_series,
#'                            lower  = "County",
#'                            upper  = "Region",
#'                            times  = "Year",
#'                            values = "Immigration",
#'                            sort_lower = "sum_ascend",
#'                            plot_spacers = 1)
#'
#' # View result
#' increase_spaces
#'
#' # Nomalise the legend
#' normalised_lgd      = tshm(df = example_time_series,
#'                            lower  = "County",
#'                            upper  = "Region",
#'                            times  = "Year",
#'                            values = "Immigration",
#'                            sort_lower = "sum_ascend",
#'                            norm_lgd = TRUE)
#'
#' # View result
#' normalised_lgd
#'
#' # Manually define colour scheme for heatmap (uses viridis colour scheme)
#' viridis_12 = c("#440154FF","#482173FF","#433E85FF","#38598CFF","#2D708EFF","#25858EFF",
#'                "#1E9B8AFF","#2BB07FFF","#51C56AFF","#85D54AFF","#C2DF23FF","#FDE725FF")
#'
#' # Assign continuous colour scheme
#' cont_clrs           = tshm(df = example_time_series,
#'                            lower  = "County",
#'                            upper  = "Region",
#'                            times  = "Year",
#'                            values = "Immigration",
#'                            sort_lower = "sum_ascend",
#'                            norm_lgd = TRUE,
#'                            cclrs = viridis_12)
#'
#' # View result
#' cont_clrs
#'
#' # Assign colour for NA values
#' na_clrs             = tshm(df = example_time_series,
#'                            lower  = "County",
#'                            upper  = "Region",
#'                            times  = "Year",
#'                            values = "Immigration",
#'                            sort_lower = "sum_ascend",
#'                            norm_lgd = TRUE,
#'                            cclrs = viridis_12,
#'                            na_colour = "grey80")
#'
#' # View result
#' na_clrs
#'
#' # Break legends into a specified number of bins
#' # (of equal intervals between 0 and the maximum value in `values`)
#' bins_15             = tshm(df = example_time_series,
#'                            lower  = "County",
#'                            upper  = "Region",
#'                            times  = "Year",
#'                            values = "Immigration",
#'                            sort_lower = "sum_ascend",
#'                            bins = 15)
#'
#' # View result
#' bins_15
#'
#' # Manually break data into categories using user-specified intervals.
#' # In this instance, the `hhmR` function `log_seq` has been used to create a
#' # vector of logarithmicly increasing values between 1 and the maximum value
#' # in the dataset.
#' cbrks = log_seq(example_time_series %>% .$Immigration %>% max(na.rm = TRUE),
#'                 12, rmv_extremes = TRUE)
#'
#' # Manually assign legend categories
#' legend_cats         = tshm(df = example_time_series,
#'                            lower  = "County",
#'                            upper  = "Region",
#'                            times  = "Year",
#'                            values = "Immigration",
#'                            sort_lower = "sum_ascend",
#'                            cbrks = cbrks)
#'
#' # View result
#' legend_cats
#'
#' # Manually assign colours to legend categories
#' cat_clrs            = tshm(df = example_time_series,
#'                            lower  = "County",
#'                            upper  = "Region",
#'                            times  = "Year",
#'                            values = "Immigration",
#'                            sort_lower = "sum_ascend",
#'                            cbrks = cbrks,
#'                            cclrs = viridis_12,
#'                            na_colour = "grey80")
#'
#' # View result
#' cat_clrs
#'
#' # Manually define order of x-axis and groups using factor levels
#' new_time_series = example_time_series %>%
#'                   mutate(Year   = factor(Year,
#'                                          levels = c(2012,2011,2014,
#'                                                     2013,2015)),
#'                          Region = factor(Region,
#'                                          levels = c("North","Midlands",
#'                                                     "South West",
#'                                                     "South East")))
#'
#' # Manually define order of x-axis and groups
#' rearrange_axes      = tshm(df = new_time_series,
#'                            lower  = "County",
#'                            upper  = "Region",
#'                            times  = "Year",
#'                            values = "Immigration",
#'                            sort_lower = "sum_ascend",
#'                            cbrks = cbrks,
#'                            cclrs = viridis_12,
#'                            na_colour = "grey80")
#'
#' # View result
#' rearrange_axes
#' @export
#' @importFrom dplyr group_split
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom purrr map
#' @importFrom purrr list_flatten
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_minimal
#' @importFrom grid unit
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 coord_cartesian
#' @importFrom patchwork plot_spacer
#' @importFrom rlang sym
#' @importFrom rlang .data
#' @importFrom utils globalVariables
tshm = function(df,lower,upper,times,values,sort_lower="alphabetical",lgttl=NULL,bins=NULL,cbrks=NULL,cclrs=NULL,norm_lgd=F,lgdps=0,na_colour=NULL,xttl_height=0.05,yttl_width=0.15) {

  # Define max value supplied to `values`
  max_value = max(df[[values]], na.rm = TRUE)

  # Check that supplied model inputs are compatible and won't cause errors
  if (!is.null(bins) && !is.null(cbrks)) { stop("The inputs bins and cbrks should not be supplied at the same time.
bins is used to break the data into a specific number of groups with equal intervals between the min and max values.
cbrks is used to manually break the data into groups based on the supplied thresholds.
Please provide either one or the other.") }
  if (!is.null(bins) && !is.null(cclrs)) {
    if (bins != length(cclrs)) { stop("If both bins and cclrs are provideds, bins and cclrs must both be vectors with cclrs having a length equal to the value of bins.") }
  }
  if (!is.null(cbrks) && !is.null(cclrs)) {
    if ( length(cbrks) != (length(cclrs)-2) ) { stop("If both cbrks and cclrs are provided, cbrks and cclrs must both be vectors with cclrs having a length two longer than cbrks.") }
  }
  if (!is.null(cbrks) && norm_lgd) {
    if (min(range(cbrks)) <  0 || max(range(cbrks)) > 1) { stop("If normalising the values (norm_lgd == TRUE), all breaks provided to cbrks must be between 0 and 1.") }
  }
  if (!is.null(cbrks) && (cbrks %>% diff() %>% {. <= 0} %>% sum() %>% {. > 0})) { stop("Please ensure the values in cbrks are provided in ascending order.") }

  # Remove unwanted rows and format origin so geographies appear in alphabetical order
  df = df[,c(lower,upper,times,values)]

  # If no colour has been assigned for NA values, then remove them from the dataset
  if (is.null(na_colour)) {
    df = df %>% filter(!is.na(!!rlang::sym(values)))
  }

  # Define the groups to be shown along the y-axis
  # If ordering of groups already defined via factor ordering, take this as the order
  # the groups should appear (top to bottom). Otherwise, order groups alphabetically
  if (!is.null(df[[upper]] %>% levels())) {
    ygrps = df[[upper]] %>% levels()
  } else {
    ygrps = df[[upper]] %>% unique() %>% sort()
  }

  # Option to normalise values between 0-1 (only works if all values are positive)
  if (norm_lgd) {

    # If any values are negative, return error message
    if ((df[[values]] < 0) %>% sum(na.rm = TRUE) %>% {. > 0}) {stop("norm_lgd is only designed to be used if all values used to populate the heatmap are positive.")}

    # Otherwise normalise values
    df[[values]] = df[[values]] / max(df[[values]], na.rm = TRUE)

    # Unless a value other than zero is supplied (i.e. the user has manually specified a non-default value), set the number of decimal points shown in the legend to 3
    if (lgdps == 0) {
      lgdps = 3
    }
  }

  # Option to split legend into custom categories
  if (!is.null(cbrks)) { # If cbrks provided

    # Add the smallest value possible in R as a lower threshold to cbrks.
    # This ensures anything that is equal to, or less than, zero is included in the first group.
    cbrks = c(.Machine$double.xmin,cbrks)

    # Define names of custom breaks
    if (lgdps == 0) {
      brk_nms = c(0, paste(cbrks %>% .[2:length(.)] %>% c(0    ,.    ),
                           cbrks %>% .[2:length(.)] %>% c(.,max_value),
                           sep = "-"))
    } else if (norm_lgd) {
      brk_nms = c(0, paste(cbrks %>% .[2:length(.)] %>% c(0    ,.    ) %>% sprintf(fmt = paste0('%#.',lgdps,'f')),
                           cbrks %>% .[2:length(.)] %>% c(.,1        ) %>% sprintf(fmt = paste0('%#.',lgdps,'f')),
                           sep = "-"))
    } else {
      brk_nms = c(0, paste(cbrks %>% .[2:length(.)] %>% c(0    ,.    ) %>% sprintf(fmt = paste0('%#.',lgdps,'f')),
                           cbrks %>% .[2:length(.)] %>% c(.,max_value) %>% sprintf(fmt = paste0('%#.',lgdps,'f')),
                           sep = "-"))
    }

    # Create backup version of origin variable so rows can still be ordered by row total or row means
    if (sort_lower != "alphabetical") {
      df[[paste0(values,"_old")]] = df[[values]]
    }

    # Create discrete scale based on these custom breaks
    df[[values]] = df[[values]] %>% findInterval(cbrks) %>% {. + 1} %>% addNA() %>% factor(levels = 1:length(brk_nms), labels = brk_nms)

  } else if (!is.null(bins)) { # If bins provided

    # Assign breaks to be equidistant thresholds between zero and maximum observed values
    # Also add the minimum possible value above zero as the first break in the sequence
    cbrks = seq(0, max(df[[values]], na.rm = TRUE), length.out = bins - 1) %>% .[2:(length(.)-1)] %>% c(.Machine$double.xmin,.)

    # Define names of custom breaks
    if (lgdps == 0) { # If set to show whole numbers
      # Round all values other than the first one (which is the minimum possible value above zero), up to the nearest whole number
      cbrks = c(.Machine$double.xmin, cbrks %>% .[2:length(.)] %>% ceiling())
      # Assign break names between zero and the maximum observed value in the data
      brk_nms = c(0, paste(cbrks %>% .[2:length(.)] %>% c(0    ,.    ),
                           cbrks %>% .[2:length(.)] %>% c(.,max_value),
                           sep = "-"))
    } else if (norm_lgd) { # If data has been normalised
      # Assign break names between 0 and 1 to the specified number of decimal points (lgdps)
      brk_nms = c(0, paste(cbrks %>% .[2:length(.)] %>% c(0    ,.    ) %>% sprintf(fmt = paste0('%#.',lgdps,'f')),
                           cbrks %>% .[2:length(.)] %>% c(.,1        ) %>% sprintf(fmt = paste0('%#.',lgdps,'f')),
                           sep = "-"))
    } else { # If using non-rounded, non-normalised data
      # Assign break names between zero and the maximum observed value in the data to the specified number of decimal points (lgdps)
      brk_nms = c(0, paste(cbrks %>% .[2:length(.)] %>% c(0    ,.    ) %>% sprintf(fmt = paste0('%#.',lgdps,'f')),
                           cbrks %>% .[2:length(.)] %>% c(.,max_value) %>% sprintf(fmt = paste0('%#.',lgdps,'f')),
                           sep = "-"))
    }

    # Create backup version of origin variable so rows can still be ordered by row total or row means
    if (sort_lower != "alphabetical") {
      df[[paste0(values,"_old")]] = df[[values]]
    }

    # Create discrete scale based on these custom breaks
    df[[values]] = df[[values]] %>% findInterval(cbrks) %>%
      replace(. == length(cbrks), (length(cbrks)-.Machine$double.xmin)) %>%
      {. + 1} %>% addNA() %>%
      factor(levels = 1:length(brk_nms), labels = brk_nms)

  } else { # Otherwise define consistent legend scale range
    lg_lims = df[[values]] %>% range(na.rm = TRUE)
  }

  # Define legend title (if not defined by user)
  if (is.null(lgttl) && norm_lgd) {
    lgttl = "Normalised\nValues"
  } else if (is.null(lgttl) && norm_lgd == F) {
    lgttl = "Values"
  }

  # Create empty list to populate with ggplot heatmaps
  pl = list()

  # Define vectors capturing the number of lower categories in each upper group
  yglns = df %>% group_split(!!rlang::sym(upper)) %>% purrr::map(~ .[[lower]] %>% unique() %>% length()) %>% unlist()

  # For each y-axis group
  for (ygrp in 1:length(ygrps)) {

    # Filter group-level migration data to only include origin and destination regions of interest
    sdf = df[df[[upper]] == ygrps[ygrp],]

    # If breaking data into categories, and sorting values by row values, then use continuous version of values to sort rows
    if (is.null(cbrks) && sort_lower != "alphabetical") {
      sort_var = values
    } else if (!is.null(cbrks) && sort_lower != "alphabetical") {
      sort_var = paste0(values,"_old")
    }

    # Define how rows of the lower categories should be arranged
    if (sort_lower == "alphabetical") {
      sdf[[lower]] = factor(sdf[[lower]], levels = sdf[[lower]] %>% unique() %>% sort() %>% rev() )
    } else if (sort_lower == "sum_ascend") {
      sdf[[lower]] = factor(sdf[[lower]], levels = sdf %>% group_by(!!rlang::sym(lower)) %>%
                              summarise(sum = sum(!!rlang::sym(sort_var), na.rm = TRUE)) %>%
                              arrange(sum) %>% .[[lower]] %>% rev() )
    } else if (sort_lower == "sum_descend") {
      sdf[[lower]] = factor(sdf[[lower]], levels = sdf %>% group_by(!!rlang::sym(lower)) %>%
                              summarise(sum = sum(!!rlang::sym(sort_var), na.rm = TRUE)) %>%
                              arrange(sum) %>% .[[lower]] )
    } else if (sort_lower == "mean_ascend") {
      sdf[[lower]] = factor(sdf[[lower]], levels = sdf %>% group_by(!!rlang::sym(lower)) %>%
                              summarise(mean = mean(!!rlang::sym(sort_var), na.rm = TRUE)) %>%
                              arrange(mean) %>% .[[lower]] %>% rev() )
    } else if (sort_lower == "mean_descend") {
      sdf[[lower]] = factor(sdf[[lower]], levels = sdf %>% group_by(!!rlang::sym(lower)) %>%
                              summarise(mean = mean(!!rlang::sym(sort_var), na.rm = TRUE)) %>%
                              arrange(mean) %>% .[[lower]] )
    } else {
      stop("The variable sort_lower should be defined as one of the following: alphabetical, sum_ascend, sum_descend, mean_ascend, mean_descend.
See `?tshm` for details.")
    }

    # Define main plot
    p = ggplot(sdf, aes(.data[[times]], .data[[lower]])) +
      geom_tile(aes(fill = .data[[values]]), show.legend = TRUE) +
      theme_minimal() +
      theme(plot.margin = unit(rep(0,4), "cm"),
            axis.title.y = element_text(angle =  0, hjust = 20.0, vjust = 0.5),
            axis.ticks   = element_blank()) +
      #labs(y = ygrps[ygrp])
      labs(x = NULL, y = NULL)

    # Define colour scale
    if (is.null(na_colour)) {
      if (!is.null(cbrks) && !is.null(cclrs)) {
        p = p + scale_fill_manual(name = lgttl, values = cclrs                                   , drop = F, na.translate = FALSE)
      } else if (!is.null(cbrks) && is.null(cclrs)) {
        p = p + scale_fill_manual(name = lgttl, values = cg("grey90","#08306B",(length(cbrks)+1)), drop = F, na.translate = FALSE)
      } else if (is.null(cbrks) && !is.null(cclrs) && norm_lgd) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cclrs                 , limits = c(0,1) )
      } else if (is.null(cbrks) && !is.null(cclrs) && norm_lgd == F) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cclrs                 , limits = lg_lims)
      } else if (is.null(cbrks) && is.null(cclrs) && norm_lgd) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cg("grey90","#08306B"), limits = c(0,1) )
      } else if (is.null(cbrks) && is.null(cclrs) && norm_lgd == F) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cg("grey90","#08306B"), limits = lg_lims)
      }
    } else {
      if (!is.null(cbrks) && !is.null(cclrs)) {
        p = p + scale_fill_manual(name = lgttl, values = cclrs                                   , drop = F, na.value = na_colour)
      } else if (!is.null(cbrks) && is.null(cclrs)) {
        p = p + scale_fill_manual(name = lgttl, values = cg("grey90","#08306B",(length(cbrks)+1)), drop = F, na.value = na_colour)
      } else if (is.null(cbrks) && !is.null(cclrs) && norm_lgd) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cclrs                 , limits = c(0,1) , na.value = na_colour)
      } else if (is.null(cbrks) && !is.null(cclrs) && norm_lgd == F) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cclrs                 , limits = lg_lims, na.value = na_colour)
      } else if (is.null(cbrks) && is.null(cclrs) && norm_lgd) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cg("grey90","#08306B"), limits = c(0,1) , na.value = na_colour)
      } else if (is.null(cbrks) && is.null(cclrs) && norm_lgd == F) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cg("grey90","#08306B"), limits = lg_lims, na.value = na_colour)
      }
    }

    # To avoid creating multiple legends, if the user has specified a colour for NA values then do not show any legends
    # for plots where no NA values occured. This assumed that at least one plot the multiplot contains NA values, which seems
    # reasonable, as why else would the user bother to specify the colour of NA values)
    if (!is.null(na_colour) && (sdf[[values]] %>% is.na() %>% sum() %>% {. == 0}) ) {
      p = p + theme(legend.position = "none")
    }

    # If not bottom plot, remove x-axis text and title
    if (ygrp != length(ygrps)) {
      p = p + theme(axis.title.x = element_blank(),
                    axis.text.x  = element_blank())
    }

    # Add ggplot to plot list
    pl[[ygrp]] = p

  }

  # Define plot heights and widths (including group titles)
  wds = c(yttl_width,1)
  hts = c(yglns,(sum(yglns)*xttl_height))

  # Define plot spacer
  ps = plot_spacer()

  # Create empty list to be populated with plot titles
  yttls = list()

  # Define plot titles
  for (ygrp in 1:length(ygrps)) {
    yttls[[ygrp]] = plt_ttl(ygrps[ygrp],axs="y")
  }

  # Create empty list to populate with both plot title and heatmap tiles in the correct order
  plts = list()

  # Define counters for subsetting plot title and heatmap lists, to ensure they are ordered correctly
  i = 1
  j = 1

  # For each group (row), assign each plot title, then the heatmap tiles within that row to plts list
  for (ygrp in 1:length(yttls)) {

    # Add plot title to list
    plts[[i]] = yttls[[ygrp]]

    # Add heatmap plots to list
    plts[i+1] = pl[j]

    # Adjust counters
    i = i + 2
    j = j + 1
  }

  # Add x-axis plots
  plts[[length(plts)+1]] = ps
  plts[[length(plts)+1]] = plt_ttl(times, rotate_title = FALSE)

  # Define final plot
  plt = patchwork::wrap_plots(plts, widths = wds, heights = hts, guides = "collect")

  # Return final plot
  return(plt)
}
