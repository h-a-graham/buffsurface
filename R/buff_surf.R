#' Buffer a point along a terrain surface.
#'
#' @param dtm A SpatRaster Digital Terrain Model
#' @param pnts A SpatVector of points
#' @param .dist Numeric - the buffer distance.
#' @param .res Numeric - approximate resampling resolution of the (internally
#' calculated) cost raster.
#' @param .smooth Logical default FALSE. Should the output buffer be smoothed using `smoothr::smooth`?
#' @param .method The smoothing method to use - see `smoothr::smooth` for details.
#' @param ... Passed to `smoothr::smooth` for eg. to control smoothness.
#'
#' @return
#' @export
#'
#' @examples
buff_surf <-
  function(dtm,
           pnts,
           .dist = 30,
           .res = 1,
           .smooth = FALSE,
           .method = c("ksmooth", "chaikin", "spline"),
           ...) {

    spv_list <- lapply(1:length(pnts), function(x) {
      dtm_crop <-
        terra::crop(dtm, terra::buffer(pnts[x], .dist + .dist * 0.1))

      .fac <- round(terra::res(dtm_crop)[1] / .res)

      dtm_crop_res <- terra::disagg(dtm_crop, .fac, method = "bilinear")

      slp <-
        terra::terrain(
          dtm_crop_res,
          v = "slope",
          unit = c("radians"),
          neighbors = 8,
          overwrite = TRUE
        )

      cost_rast <- (res(slp)[1] / cos(slp)) / res(slp)[1]

      extr_cells <- terra::extract(cost_rast, pnts[x], cells = TRUE)

      cost_rast[extr_cells$cell] <- 0

      cost_dist <- costDist(cost_rast)

      cost_dist[cost_dist >= (.dist + .res*0.5)] <- NA
      cost_dist[!is.na(cost_dist)] <- 1

      cost_ch  <- as.polygons(cost_dist)

      if (isTRUE(.smooth)) {
        cost_ch <- cost_ch |>
          sf::st_as_sf() |>
          smoothr::smooth(method = .method[1], ...) |>
          terra::vect()
      }

      cost_ch$PlotNo <- pnts[x]$Pl

      return((cost_ch))
    })
    do.call(rbind, spv_list)
  }


#' Buffer a point along a terrain surface.
#' Alias for `buff_surf`
#'
#' @param dtm A SpatRaster Digital Terrain Model
#' @param pnts A SpatVector of points
#' @param .dist Numeric - the buffer distance.
#' @param .res Numeric - approximate resampling resolution of the (internally
#' calculated) cost raster.
#' @param .smooth Logical default FALSE. Should the output buffer be smoothed using `smoothr::smooth`?
#' @param .method The smoothing method to use - see `smoothr::smooth` for details.
#' @param ... Passed to `smoothr::smooth` for eg. to control smoothness.
#'
#' @return
#' @export
#'
#' @examples
buff_surface <-
  function(dtm,
           pnts,
           .dist = 30,
           .res = 1,
           .smooth = FALSE,
           .method = c("ksmooth", "chaikin", "spline"),
           ...){
    buff_surf(dtm, pnts, .dist, .res, .smooth, .method, ...)
  }
