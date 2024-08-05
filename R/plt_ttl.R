#' @title plt_ttl
#'
#' @description Creates plot containing the name of a given upper group. Used in combination with the
#' patchwork package to plot the names of the upper groups within the hhm function.
#'
#' @param ttl The name of the upper group.
#' @param axs The axis on which the name will appear (defaults to "x"). If `x`, the text will be
#' written at the top-centre of the plot. If `y`, the text will be written at the middle-right of the
#' plot.
#' @param rotate_title Whether the title should be rotate to be perpendicular to the axis (defaults
#' to TRUE). If TRUE, the title text on the x and y axes will be printed horizontally and vertically
#' respectively, with the reverse orientation if set to FALSE.
#'
#' @return A ggplot object containing the title of a given upper group, for use in the hhm function.
#'
#' @examples
#' plt_ttl("Group 1", axs = "y")
#' plt_ttl("Group 2")
#' plt_ttl("Group 1", axs = "y",rotate_title = FALSE)
#' plt_ttl("Group 2"           ,rotate_title = FALSE)
#'
#' @export
plt_ttl = function(ttl,axs="x",rotate_title=TRUE) {

  # If plotting on x-axis
  if (axs == "x") {

    # Place at top of plot
    p = ggplot(data.frame(x = 0:1, y = 0:1), aes(x = .data[["x"]], y = .data[["y"]])) +
      geom_point(col = "white") +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      theme_void() +
      theme(plot.margin = unit(rep(0,4), "cm"))

    # Speficy wehther title should be roated 90 degrees or not
    if (rotate_title) {
      p = p + geom_text(x = 0.495, y = 1.0, angle = 90, label = ttl, size = 4, hjust = 1)
    } else if (rotate_title == FALSE) {
      p = p + geom_text(x = 0.495, y = 0.9, angle =  0, label = ttl, size = 4, hjust = 1)
    } else {
      stop("rotate_title must be TRUE or FALSE.")
    }

  } else if (axs == "y") { # If plotting on y-axis

    # Place at left of plot
    p = ggplot(data.frame(x = 0:1, y = 0:1), aes(x = .data[["x"]], y = .data[["y"]])) +
      geom_point(col = "white") +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      theme_void() +
      theme(plot.margin = unit(rep(0,4), "cm"))

    # Speficy wehther title should be roated 90 degrees or not
    if (rotate_title) {
      p = p + geom_text(x = 1, y = 0.51, angle =  0, label = ttl, size = 4, hjust = 1)
    } else if (rotate_title == FALSE) {
      p = p + geom_text(x = 1, y = 0.51, angle = 90, label = ttl, size = 4, hjust = 1)
    } else {
      stop("rotate_title must be TRUE or FALSE.")
    }

  }

  # Return region name plot
  return(p)
}
