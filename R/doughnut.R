#' @title Donut/Doughnut plot
#'
#' @description Generation of a donut/doughnut plot (equivalent of a pie chart)
#'
#' @references \url{https://magesblog.com/}
#'
#' @param x A vector containing the values to be plotted.
#' @param labels A string vector for the labels of the different sectors. By default as.character(x).
#' @param edges Number of edges of the shape. By default 200.
#' @param outer.radius Fraction of the area to dedicate to the outer circle. By default 0.8.
#' @param inner.radius Fraction of the area to dedicate to the inner circle. By default 0.4.
#' @param clockwise Logic value to define whether the values should be plotted in clockwise sense. By default \code{FALSE}.
#' @param init.angle Numeric value to define the starting angle for the data. By default if \code{clockwise = TRUE} 90, otherwise 0.
#' @param density A vector or single number to define de density of the lines in the filling color of each value plotted. By default \code{NULL}.
#' @param angle A vector or single number to define de angle of the lines in the filling color of each value plotted. By default 45.
#' @param col A vector of R standard colors for each value to be plotted. By default \code{NULL}.
#' @param border Logic value to define whether plot the border of the sectors. By default \code{FALSE}.
#' @param lty Numeric value to define the type of line for the borders. By default \code{NULL}.
#' @param main String to set the title of the plot. By default \code{NULL}.
#'
#' @return
#'
#' @examples
#' doughnut(x = c(3,5,9,12), inner.radius=0.5, col=c("red", "blue", "green", "yellow"))
#'
#' @export doughnut

doughnut =
  function (x,
            labels = as.character(x),
            edges = 200,
            outer.radius = 0.8,
            inner.radius = 0.4,
            clockwise = FALSE,
            init.angle = if (clockwise) 90 else 0,
            density = NULL,
            angle = 45,
            col = NULL,
            border = FALSE,
            lty = NULL,
            main = NULL,
            ...)
  {
    ### Open a null PDF to record the file
    grid::grid.newpage()
    pdf(NULL)
    dev.control(displaylist="enable")

    if (!is.numeric(x) || any(is.na(x) | x < 0))
      stop("'x' values must be positive.")
    if (is.null(labels))
      labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L])
      xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    plot.window(xlim, ylim, "", asp = 1)
    if (is.null(col))
      col <- if (is.null(density))
        palette()
    else par("fg")
    col <- rep(col, length.out = nx)
    border <- rep(border, length.out = nx)
    lty <- rep(lty, length.out = nx)
    angle <- rep(angle, length.out = nx)
    density <- rep(density, length.out = nx)
    twopi <- if (clockwise)
      -2 * pi
    else 2 * pi
    t2xy <- function(t, radius) {
      t2p <- twopi * t + init.angle * pi/180
      list(x = radius * cos(t2p),
           y = radius * sin(t2p))
    }
    for (i in 1L:nx) {
      n <- max(2, floor(edges * dx[i]))
      P <- t2xy(seq.int(x[i], x[i + 1], length.out = n),
                outer.radius)
      polygon(c(P$x, 0), c(P$y, 0), density = density[i],
              angle = angle[i], border = border[i],
              col = col[i], lty = lty[i])
      Pout <- t2xy(mean(x[i + 0:1]), outer.radius)
      lab <- as.character(labels[i])
      if (!is.na(lab) && nzchar(lab)) {
        lines(c(1, 1.05) * Pout$x, c(1, 1.05) * Pout$y)
        text(1.1 * Pout$x, 1.1 * Pout$y, labels[i],
             xpd = TRUE, adj = ifelse(Pout$x < 0, 1, 0),
             ...)
      }
      ## Add white disc
      Pin <- t2xy(seq.int(0, 1, length.out = n*nx),
                  inner.radius)
      polygon(Pin$x, Pin$y, density = density[i],
              angle = angle[i], border = border[i],
              col = "white", lty = lty[i])
    }

    title(main = main, ...)
    invisible(NULL)

    final_plot = recordPlot()
    invisible(dev.off())

    return(final_plot)
  }
