#' @title plt_ttl
#'
#' @description Creates plot containing the name of a given upper group. Used in combination with the
#' patchwork package to plot the names of the upper groups within the hhm function.
#'
#' @param ttl The name of the upper group.
#' @param axs The axis on which the name will appear. If `x`, the text will be written at the
#' top-centre of the plot. If `y`, the text will be written at the middle-right of the plot.
#'
#' @return A ggplot object containing the title of a given upper group, for use in the hhm function.
#'
#' @examples
#' plt_ttl("Group 1",axs="y")
#' plt_ttl("Group 2")
#'
#' @export
plt_ttl = function(ttl,axs="x") {

  # If plotting on x-axis
  if (axs == "x") {

    # Place at top of plot and rotate 90 degrees
    p = ggplot(data.frame(x = 0:1, y = 0:1), aes(x = .data[["x"]], y = .data[["y"]])) +
      geom_point(col = "white") +
      geom_text(x = 0.495, y = 1, angle = 90, label = ttl, size = 4, hjust = 1) +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      theme_void() +
      theme(plot.margin = unit(rep(0,4), "cm"))

  } else if (axs == "y") { # If plotting on y-axis

    # Place at left of plot
    p = ggplot(data.frame(x = 0:1, y = 0:1), aes(x = .data[["x"]], y = .data[["y"]])) +
      geom_point(col = "white") +
      geom_text(x = 1, y = 0.51, angle =  0, label = ttl, size = 4, hjust = 1) +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      theme_void() +
      theme(plot.margin = unit(rep(0,4), "cm"))

  }

  # Return region name plot
  return(p)
}
