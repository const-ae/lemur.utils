
ScaleContinuous2 <- ggplot2::ggproto("ScaleContinuous2", ggplot2::ScaleContinuous,
  train = function(self, x){
    ggplot2::ggproto_parent(ggplot2::ScaleContinuous, self)$train(x)
    if(! is.null(self$qlimits)){
      qlims <- self$qlimits
      abs_limit <- sort(quantile(x, qlims))
      old_range <- self$range$range
      midpoint <- if(is.null(self$midpoint)){
        0
      }else{
        self$midpoint
      }
      self$range$range <- c(min(midpoint, max(abs_limit[1], old_range[1])),
                            max(midpoint, min(abs_limit[2], old_range[2])))
    }
  }
)



#' A divergent color scale designed to display differential expression values
#'
#' Creates a divergent colorscale around the `midpoint`. Unlike the
#' [ggplot2::scale_color_gradient2()] blue signals low values and red
#' are high values. In addition, `qlimits` filters out outliers that
#' distort the overall picture.
#'
#' @param qlimits the quantile limits. Alternative to limits, which allows
#'   specification of the quantiles to be displayed. By default, it is
#'   combined with `oob = scales::squish` to squish outliers into a
#'   the desired range. The argument has to be either `NULL`, in which
#'   case it is ignored, a numeric vector with two elements between 0 and 1,
#'   or a numeric vector with only a single element, in which case
#'   the limits are `c(qlimits, 1-qlimits)`.
#'
#' @seealso This function is based on [ggplot2::scale_color_gradient2()].
#'   See there, for details about the other parameter.
#' @examples
#'   library(ggplot2)
#'   library(tibble)
#'
#'   dat <- tibble(x = 1:10, y = seq(-3, 5, length.out = 10))
#'   ggplot(dat, aes(x = x, y = y)) +
#'     geom_col(aes(fill = y)) +
#'     scale_fill_de()
#'
#'   ggplot(dat, aes(x = x, y = y)) +
#'     geom_col(aes(fill = y)) +
#'     scale_fill_de(qlimits = c(0, 0.7))
#'
#'   dat <- tibble(x = 1:10, y = seq(-1, 3, length.out = 10))
#'   ggplot(dat, aes(x = x, y = exp(y))) +
#'     geom_col(aes(fill = y)) +
#'     scale_fill_de(midpoint = 1, qlimits = c(0.5, 1))
#'
#' @export
scale_color_de <- function(...){
  scale_de(..., aesthetics = "colour")
}

#' @rdname scale_color_de
#' @export
scale_colour_de <- function(...){
  scale_de(..., aesthetics = "colour")
}

#' @rdname scale_color_de
#' @export
scale_fill_de <- function(...){
  scale_de(..., aesthetics = "fill")
}

#' @rdname scale_color_de
#' @export
scale_de <- function(name = ggplot2::waiver(), ...,
                     qlimits = c(0, 1),
                     low = munsell::mnsl("10B 4/6"),
                     mid = munsell::mnsl("N 8/0"),
                     high = munsell::mnsl("10R 4/6"),
                     midpoint = 0,
                     space = "Lab",
                     na.value = "grey50",
                     oob = scales::squish,
                     transform = "identity",
                     guide = "colourbar",
                     aesthetics = "colour"){
  res <- ggplot2::continuous_scale(aesthetics = aesthetics, name = name,
                   palette = scales::pal_div_gradient(low, mid, high, space = space),
                   na.value = na.value, oob = oob, guide = guide , transform = transform, ...,
                   rescaler = mid_rescaler(mid = midpoint, transform = transform), super = ScaleContinuous2)
  res$qlimits <- if(is.null(qlimits)){
    qlimits
  }else if(length(qlimits) == 1){
    stopifnot(all(qlimits >= 0 & qlimits <= 1))
      c(qlimits, 1-qlimits)
  }else if(length(qlimits) == 2){
    stopifnot(all(qlimits >= 0 & qlimits <= 1))
    qlimits
  }else{
    stop("qlimits must be either NULL or a numeric vector with 1 or 2 elements")
  }
  stopifnot(length(midpoint) == 1 && is.numeric(midpoint))
  res$midpoint <- midpoint

  res
}


# See https://github.com/tidyverse/ggplot2/blob/2610840f7478af027294db0793df839199b8cb6b/R/scale-gradient.R#L146
mid_rescaler <- function(mid, transform = "identity"){
  transform <- scales::as.trans(transform)
  trans_mid <- transform$transform(mid)
  function(x, to = c(0, 1), from = range(x, na.rm=TRUE)){
    scales::rescale_mid(x, to, from, trans_mid)
  }
}



