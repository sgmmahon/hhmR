utils::globalVariables(".")
#' @title Hierarchical Heatmap
#'
#' @description Creates a labelled heatmap from heirarchical data. This function is
#' useful if you wish to create a heatmap where the categories shown on both the x
#' and y axis can be grouped in some way. This heatmap will order the categories by
#' their assigned group and present both the categories and group labels along the
#' axes. An example might be series of smaller geographies (lower categories) which
#' aggregate into larger geographical regions (upper groups).
#'
#' @param df A data.frame with containing values with which to populate the heatmap.
#' The data.frame must include columns specifying the lower categories (`ylower`,
#' `xlower`) and upper groups (`yupper`, `xupper`) that each value corresponds to.
#' These categories and groups will be used to arrange and label the rows and
#' columns of the heatmap. It must also contain a `values` variable containing the
#' values used to populate the heatmap. Note that the groups will by default be
#' arranged alphabetically (top to bottom / left to right). The ordering of the
#' groups can be manually specified by converting yupper and/or xupper to factors.
#' In this case, the groups will be ordered based on the ordering of the factor
#' levels.
#' @param ylower A column in `df` containing the categories that will be presented
#' as rows along the y-axis of the heatmap.
#' @param xlower A column in `df` containing the categories that will be presented
#' as columns along the x-axis of the heatmap.
#' @param yupper A column in `df` containing the groupings that will be used to
#' arrange the heatmap rows.
#' @param xupper A column in `df` containing the groupings that will be used to
#' arrange the heatmap columns.
#' @param values A column in `df` containing the values used to populate the
#' heatmap.
#' @param rm_diag Do not show values for categories along the x and y axes that
#' are identical (defaults to `FALSE`). This is particularly useful for
#' origin-destination heatmaps, where the user may want to hide the diagonal
#' values.
#' @param lgttl Option to manually define legend title.
#' @param bins Option to break the data into a specified number of groups
#' (defaults to `NULL`). The thresholds between these groups will be equally
#' spaced between zero and the maximum value observed in `values`.
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
#' @param xttl_height The space allocated to the group titles on the x-axis as a
#' proportion of the heatmap's height (defaults to 0.15).
#' @param yttl_width The space allocated to the group titles on the y-axis as a
#' proportion of the heatmap's width (defaults to 0.15).
#'
#' @return A ggplot object containing the final heatmap.
#'
#' @examples
#' # Import toy demonstration dataset (see `?example_migration` for see details)
#' data(example_migration)
#'
#' # Intial heatmap
#' hierarchical_heatmap = hhm(df = example_migration,
#'                            ylower = "Origin County",
#'                            xlower = "Destination County",
#'                            yupper = "Origin Region",
#'                            xupper = "Destination Region",
#'                            values = "Migration",
#'                            yttl_width = 0.22,
#'                            xttl_height = 0.4)
#'
#' # View result
#' hierarchical_heatmap
#'
#' # Remove diagonal from heatmap (i.e. hide static populations)
#' removed_diag         = hhm(df = example_migration,
#'                            ylower = "Origin County",
#'                            xlower = "Destination County",
#'                            yupper = "Origin Region",
#'                            xupper = "Destination Region",
#'                            values = "Migration",
#'                            yttl_width = 0.22,
#'                            xttl_height = 0.4,
#'                            rm_diag = TRUE)
#'
#' # View result
#' removed_diag
#'
#' # Nomalise the legend
#' normalised_lgd       = hhm(df = example_migration,
#'                            ylower = "Origin County",
#'                            xlower = "Destination County",
#'                            yupper = "Origin Region",
#'                            xupper = "Destination Region",
#'                            values = "Migration",
#'                            yttl_width = 0.22,
#'                            xttl_height = 0.4,
#'                            rm_diag = TRUE,
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
#' cont_clrs            = hhm(df = example_migration,
#'                            ylower = "Origin County",
#'                            xlower = "Destination County",
#'                            yupper = "Origin Region",
#'                            xupper = "Destination Region",
#'                            values = "Migration",
#'                            yttl_width = 0.22,
#'                            xttl_height = 0.4,
#'                            rm_diag = TRUE,
#'                            norm_lgd = TRUE,
#'                            cclrs = viridis_12)
#'
#' # View result
#' cont_clrs
#'
#' # Break legends into a specified number of bins
#' # (of equal intervals between 0 and the maximum value in `values`)
#' bins_15              = hhm(df = example_migration,
#'                            ylower = "Origin County",
#'                            xlower = "Destination County",
#'                            yupper = "Origin Region",
#'                            xupper = "Destination Region",
#'                            values = "Migration",
#'                            yttl_width = 0.22,
#'                            xttl_height = 0.4,
#'                            rm_diag = TRUE,
#'                            bins = 15)
#'
#' # View result
#' bins_15
#'
#' # Manually break data into categories using user-specified intervals.
#' # In this instance, the `hhmR` function `log_seq` has been used to create a
#' # vector of logarithmicly increasing values between 1 and the maximum value
#' # in the dataset not on the diagonal.
#' cbrks = log_seq(example_migration[example_migration[["Origin County"     ]] !=
#'                                   example_migration[["Destination County"]],] %>%
#'                 .$Migration %>% max(), 12, rmv_extremes = TRUE)
#'
#' # Manually assign legend categories
#' legend_cats          = hhm(df = example_migration,
#'                            ylower = "Origin County",
#'                            xlower = "Destination County",
#'                            yupper = "Origin Region",
#'                            xupper = "Destination Region",
#'                            values = "Migration",
#'                            yttl_width = 0.22,
#'                            xttl_height = 0.4,
#'                            rm_diag = TRUE,
#'                            cbrks = cbrks)
#'
#' # View result
#' legend_cats
#'
#' # Manually assign colours to legend categories
#' cat_clrs             = hhm(df = example_migration,
#'                            ylower = "Origin County",
#'                            xlower = "Destination County",
#'                            yupper = "Origin Region",
#'                            xupper = "Destination Region",
#'                            values = "Migration",
#'                            yttl_width = 0.22,
#'                            xttl_height = 0.4,
#'                            rm_diag = TRUE,
#'                            cbrks = cbrks,
#'                            cclrs = viridis_12)
#'
#' # View result
#' cat_clrs
#' @export
#' @importFrom dplyr group_split
#' @importFrom purrr map
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_void
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
hhm = function(df,ylower,yupper,xlower,xupper,values,rm_diag=F,lgttl=NULL,bins=NULL,cbrks=NULL,cclrs=NULL,norm_lgd=F,lgdps=0,xttl_height=0.15,yttl_width=0.15) {

  # Define max value supplied to `values`
  if (rm_diag) {
    max_value = max(df[df[[xlower]] != df[[ylower]],values], na.rm = TRUE)
  } else {
    max_value = max(df[[values]], na.rm = TRUE)
  }

  # Check that supplied model inputs are compatible and won't cause errors
  if (!is.null(bins) && !is.null(cbrks)) { stop("The inputs bins and cbrks should not be supplied at the same time.
bins is used to break the data into a specific number of groups with equal intervals between the min and max values.
cbrks is used to manually break the data into groups based on the supplied thresholds.
Please provide either one or the other.") }
  if (!is.null(bins) && !is.null(cclrs)) {
    if (bins != length(cclrs)) { stop("If both bins and cclrs are provideds, bins and cclrs must both be vectors with cclrs having a length equal to the value of bins.") }
  }
  if (!is.null(cbrks) && rm_diag) {
    if ( (min(cbrks) <= 0) || (max(cbrks) >= max_value) ) { stop(paste0("All values in cbrks must be between 0 and the largest value provided to `values`.
In this instance rm_diag == TRUE, so only values not on the diagonal are considered.
All values provided to cbrks should therefore be between greater than 0 and less than ",max_value,".")) }
  }
  if (!is.null(cbrks) && rm_diag == F) {
    if ( !is.null(cbrks) && (min(cbrks) <= 0) || (max(cbrks) >= max_value) ) { stop(paste0("All values in cbrks must be between 0 and the largest value provided to `values`.
In this instance all values provided to cbrks should therefore be between greater than 0 and less than ",max_value,".")) }
  }
  if (!is.null(cbrks) && !is.null(cclrs)) {
    if ( length(cbrks) != (length(cclrs)-2) ) { stop("If both cbrks and cclrs are provided, cbrks and cclrs must both be vectors with cclrs having a length two longer than cbrks.") }
  }
  if (!is.null(cbrks) && norm_lgd) {
    if (min(range(cbrks)) <  0 || max(range(cbrks)) > 1) { stop("If normalising the values (norm_lgd == TRUE), all breaks provided to cbrks must be between 0 and 1.") }
  }
  if (!is.null(cbrks) && (cbrks %>% diff() %>% {. <= 0} %>% sum() %>% {. > 0})) { stop("Please ensure the values in cbrks are provided in ascending order.") }

  # Remove unwanted rows and format origin so geographies appear in alphabetical order
  df = df[,c(ylower,xlower,yupper,xupper,values)]

  # Define the groups to be shown along the x and y axes
  # If ordering of groups already defined via factor ordering, take this as the order
  # the groups should appear (top to bottom / left to right). # Otherwise, order groups alphabetically
  if (!is.null(df[[xupper]] %>% levels())) {
    xgrps = df[[xupper]] %>% levels()
  } else {
    xgrps = df[[xupper]] %>% unique() %>% sort()
  }
  if (!is.null(df[[yupper]] %>% levels())) {
    ygrps = df[[yupper]] %>% levels()
  } else {
    ygrps = df[[yupper]] %>% unique() %>% sort()
  }

  # If user specified to remove diagonal values, set all observations where ylower and xlower are identical to zero
  if (rm_diag) {
    df[df[[ylower]] == df[[xlower]],values] = NA
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

    # Create discrete scale based on these custom breaks
    df[[values]] = df[[values]] %>% findInterval(cbrks) %>% {. + 1} %>% addNA() %>%
      factor(levels = 1:length(brk_nms), labels = brk_nms)

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
  xglns = df %>% group_split(!!rlang::sym(xupper)) %>% purrr::map(~ .[[xlower]] %>% unique() %>% length()) %>% unlist()
  yglns = df %>% group_split(!!rlang::sym(yupper)) %>% purrr::map(~ .[[ylower]] %>% unique() %>% length()) %>% unlist()

  # Counter to keep track of interations of nested for loop
  i = 0

  # For each y-axis group
  for (ygrp in 1:length(ygrps)) {

    # For each x-axis group
    for (xgrp in 1:length(xgrps)) {

      # Increase interature counter by 1
      i = i + 1

      # Filter group-level migration data to only include origin and destination regions of interest
      sdf = df[df[[yupper]] == ygrps[ygrp] & df[[xupper]] == xgrps[xgrp],]

      # Order lower categories alphabetically
      sdf[[xlower]] = factor(sdf[[xlower]], levels = sdf[[xlower]] %>% unique() %>% sort()           )
      sdf[[ylower]] = factor(sdf[[ylower]], levels = sdf[[ylower]] %>% unique() %>% sort() %>% rev() )

      # Define main plot
      p = ggplot(sdf, aes(.data[[xlower]], .data[[ylower]])) +
        geom_tile(aes(fill = .data[[values]]), show.legend = TRUE) +
        theme(plot.margin = unit(rep(0,4), "cm"),
              axis.text.x  = element_text(angle = 90, hjust = 1.0, vjust = 0.3),
              axis.title.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
              axis.title.y = element_text(angle =  0, hjust = 0.5, vjust = 0.5),
              axis.ticks   = element_blank()) +
        labs(x = xgrps[xgrp], y = ygrps[ygrp])

      # Define colour scale
      if (!is.null(cbrks) && !is.null(cclrs)) {
        p = p + scale_fill_manual(name = lgttl, values = cclrs                                  , drop = F, na.value = "white")
      } else if (!is.null(cbrks) && is.null(cclrs)) {
        p = p + scale_fill_manual(name = lgttl, values = cg("white","#08306B",(length(cbrks)+1)), drop = F, na.value = "white")
      } else if (is.null(cbrks) && !is.null(cclrs) && norm_lgd) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cclrs                , limits = c(0,1) , na.value = "white")
      } else if (is.null(cbrks) && !is.null(cclrs) && norm_lgd == F) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cclrs                , limits = lg_lims, na.value = "white")
      } else if (is.null(cbrks) && is.null(cclrs) && norm_lgd) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cg("white","#08306B"), limits = c(0,1) , na.value = "white")
      } else if (is.null(cbrks) && is.null(cclrs) && norm_lgd == F) {
        p = p + scale_fill_gradientn(name = lgttl, colours = cg("white","#08306B"), limits = lg_lims, na.value = "white")
      }

      # To prevent legend showing NA values if rm_diag set to TRUE (in which case, diagonal set to NA), only show legend for plots that are not on the diagonal
      if (rm_diag && (sdf[[values]] %>% is.na() %>% sum() %>% {. > 0}) && (ygrp == xgrp) ) {
        p = p + theme(legend.position = "none")
      }

      # If bottom-left plot
      if (ygrp == length(ygrps) & xgrp == 1) {
        # Include provincia names on both axes
        p = p + theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank())
      } else if (ygrp < length(ygrps) & xgrp == 1) { # If left-hand plot
        # Include provincia names on y-axis
        p = p + theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.x  = element_blank())
      } else if (ygrp == length(ygrps) & xgrp > 1) { # If bottom plot
        # Include provincia names on x-axis
        p = p + theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y  = element_blank())
      } else { # If plot not on left of bottom edges of multiplot
        # Remove provincia names
        p = p + theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.x  = element_blank(),
                      axis.text.y  = element_blank())
      }

      # Add ggplot to plot list
      pl[[i]] = p

    }

  }

  # Define plot heights and widths (including group titles)
  wds = c((sum(xglns)*yttl_width),xglns)
  hts = c(yglns,(sum(yglns)*xttl_height))

  # Define plot spacer
  ps = plot_spacer()

  # Create empty lists to be populated with plot titles
  xttls = list()
  yttls = list()

  # Define plot titles
  for (xgrp in 1:length(xgrps)) {
    xttls[[xgrp]] = plt_ttl(xgrps[xgrp])
  }
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
    plts[(i+1):(i+length(xttls))] = pl[j:(j+length(xttls)-1)]

    # Adjust counters
    i = i + 1 + length(xttls)
    j = j +     length(xttls)
  }

  # Add x-axis plots
  plts[[length(plts)+1]] = ps
  plts[(length(plts)+1):(length(plts)+length(xttls))] = xttls

  # Define final plot
  plt = patchwork::wrap_plots(plts, widths = wds, heights = hts, guides = "collect")

  # Return final plot
  return(plt)
}
