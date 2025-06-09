#' Individual Tree Segmentation Algorithm
#'
#' This functions is made to be used in \link{segment_trees}. It implements an algorithm for tree
#' segmentation based on Li et al. (2012) (see reference). This method is a growing region
#' method working at the point cloud level. It is an implementation by lidR authors, from the original
#' paper, as close as possible from the original description. However we added a parameter \code{hmin}
#' to prevent over-segmentation for objects that are too low. This algorithm is known to be slow because
#' it has an algorithmic complexity worst that O(n^2).
#'
#' @param dt1 numeric. Threshold number 1. See reference page 79 in Li et al. (2012). Default is 1.5.
#' @param dt2 numeric. Threshold number 2. See reference page 79 in Li et al. (2012). Default is 2.
#' @param hmin numeric. Minimum height of a detected tree. Default is 2.
#' @param Zu numeric. If point elevation is greater than Zu, \code{dt2} is used, otherwise \code{dt1} is
#' used. See page 79 in Li et al. (2012). Default is 15.
#' @param R numeric, function, or \code{NULL}. Search radius. See page 79 in Li et al. (2012). Default is 2. If
#' \code{R = 0} all the points are automatically considered as local maxima and the search step
#' is skipped (much faster). If this value is \code{NULL}, search radius is calculated automatically
#' using the \link{lmfxauto} method (much slower). If it is a function, the function determines the size of
#' the window at any given location on the canopy. By default function takes the height of a given pixel or
#' point as its only argument and return the desired size of the search window when centered on that
#' pixel/point. This can be controlled with the 'R_args' parameter.
#' @param speed_up numeric. Maximum radius of a crown. Any value greater than a crown is
#' good because this parameter does not affect the result. However, it greatly affects the
#' computation speed by restricting the number of comparisons to perform.
#' The lower the value, the faster the method. Default is 10.
#' @param R_args list. Named list of argument for the function 'R' if 'R' is a function.
#'
#' @export
#'
#' @family individual tree segmentation algorithms
#' @family point-cloud based tree segmentation algorithms
#'
#' @references
#' Li, W., Guo, Q., Jakubowski, M. K., & Kelly, M. (2012). A new method for segmenting individual
#' trees from the lidar point cloud. Photogrammetric Engineering & Remote Sensing, 78(1), 75-84.
#'
#' @examples
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' poi <- "-drop_z_below 0 -inside 481280 3812940 481320 3812980"
#' las <- readLAS(LASfile, select = "xyz", filter = poi)
#' col <- pastel.colors(200)
#'
#' las <- segment_trees(las, li2012_auto(dt1 = 1.4))
#' #plot(las, color = "treeID", colorPalette = col)
#' @name its_li2012_auto
#' @md
li2012_auto = function(dt1 = 1.5, dt2 = 2, R = 2, Zu = 15, hmin = 2, speed_up = 10, ID = "treeID", R_args = "Z", plot = FALSE)
{
  lidR:::assert_is_a_number(dt1)
  lidR:::assert_is_a_number(dt2)
  lidR:::assert_is_a_number(Zu)
  lidR:::assert_is_a_number(hmin)
  lidR:::assert_is_a_number(speed_up)
  lidR:::assert_all_are_positive(dt1)
  lidR:::assert_all_are_positive(dt2)
  lidR:::assert_all_are_positive(Zu)
  lidR:::assert_all_are_positive(hmin)
  lidR:::assert_all_are_positive(speed_up)

  dt1      <- lazyeval::uq(dt1)
  dt2      <- lazyeval::uq(dt2)
  R        <- lazyeval::uq(R)
  Zu       <- lazyeval::uq(Zu)
  hmin     <- lazyeval::uq(hmin)
  speed_up <- lazyeval::uq(speed_up)
  R_args   <- lazyeval::uq(R_args)

  f = function(las)
  {
    lidR:::assert_is_valid_context(lidR:::LIDRCONTEXTITS, "li2012_auto")

    if (las[["Max Z"]] < hmin)
    {
      warning("'hmin' is higher than the highest point. No tree segmented.", call. = FALSE)
      return(rep(NA_integer_, nrow(las@data)))
    }
    else
    {
      if(is.null(R))
      {
        # Get R from the lmfauto_ws function
        ttop5 <- lidR::locate_trees(las, lmfx(5, hmin = hmin))
        d <- density(las)
        if(plot) 
        {
          A <- lidR::area(las)
          Aha   <- 10000/A
          ntop5 <- nrow(ttop5) * Aha
        }
        else
        {
          A     <- 400
          Aha   <- 10000/A
          x     <- sf::st_coordinates(ttop5)[, "X"]
          y     <- sf::st_coordinates(ttop5)[, "Y"]
          ntop5 <- C_count_in_disc(x, y, las@data$X, las@data$Y, sqrt(A/pi), lidR:::getThread())
          ntop5 <- ntop5 * Aha
        }
        R <- lmfauto_ws(las@data$Z, ntop5, d)
      }
      else if(is.function(R))
      {
        args <- lapply(R_args, function(x) if (x[1] %in% names(las)) las@data[[x]] else x)
        R <- do.call(R, args)
      }
      else
      {
        lidR:::assert_is_a_number(R)
      }
      return(C_li2012_auto(las, dt1, dt2, R, Zu, hmin, speed_up))
    }
  }

  f <- plugin_its(f, omp = FALSE, raster_based = FALSE)
  return(f)
}
