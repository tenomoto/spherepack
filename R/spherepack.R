.spherepack_env <- new.env(parent = emptyenv())
.spherepack_env$radius <- 6.37122e6

#' Get or Set the Sphere Radius Used by Spherical Harmonic Routines
#'
#' `get_sphere_radius()` returns the current sphere radius (in meters) used to
#' scale spatial differential and integral operators. `set_sphere_radius()` sets
#' the sphere radius (default is the Earth's mean radius, `6.37122e6` meters;
#' set to `1.0` for computations on the unit sphere).
#'
#' @param radius Numeric scalar giving the positive radius of the sphere in
#'   meters (or `1` for unit sphere).
#' @return `get_sphere_radius()` returns a single double-precision number.
#'   `set_sphere_radius()` returns `invisible(NULL)`.
#' @examples
#' orig_r <- get_sphere_radius()
#' orig_r
#' set_sphere_radius(1.0)
#' get_sphere_radius()
#' set_sphere_radius(orig_r)
#' @export
get_sphere_radius <- function() {
  .spherepack_env$radius
}

#' @rdname get_sphere_radius
#' @export
set_sphere_radius <- function(radius) {
  if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) || radius <= 0) {
    stop("'radius' must be a positive numeric scalar.")
  }
  .spherepack_env$radius <- as.double(radius)
  invisible(NULL)
}

#' Compute Gaussian Quadrature Latitudes, Colatitudes, and Weights
#'
#' Computes the `nlat` Gaussian colatitudes (in radians, North-to-South),
#' Gaussian latitudes (in degrees, South-to-North), and Gaussian quadrature
#' weights on the sphere using the Fourier-Newton method of Swarztrauber (2002).
#'
#' @param nlat Integer scalar (`>= 1`) specifying the number of Gaussian
#'   latitudes between the two poles.
#' @return A named list with components:
#' \describe{
#'   \item{`theta`}{Numeric vector of length `nlat` giving Gaussian colatitudes
#'     in radians, ordered North-to-South in `(0, pi)`.}
#'   \item{`lat`}{Numeric vector of length `nlat` giving Gaussian latitudes in
#'     degrees, ordered South-to-North in `(-90, 90)`.}
#'   \item{`wts`}{Numeric vector of length `nlat` giving Gaussian quadrature
#'     weights (summing to `2.0`).}
#' }
#' @examples
#' g <- gaqd(32)
#' length(g$lat)
#' sum(g$wts)
#' @export
gaqd <- function(nlat) {
  nlat <- as.integer(nlat)[1L]
  if (is.na(nlat) || nlat < 1L) {
    stop("'nlat' must be a positive integer.")
  }
  res <- .Call("c_gaqd", nlat, PACKAGE = "spherepack")
  lat_deg <- 90.0 - rev(res$theta) * (180.0 / pi)
  list(theta = res$theta, lat = lat_deg, wts = res$wts)
}

# Internal validation and reshaping helpers
.check_grid_array <- function(x, name = "x", min_nlon = 4L, min_nlat = 3L) {
  if (!is.numeric(x)) {
    stop(sprintf("'%s' must be a numeric array.", name))
  }
  d <- dim(x)
  if (is.null(d) || length(d) < 2L) {
    stop(sprintf("'%s' must be an array with at least 2 dimensions (nlon, nlat, ...).", name))
  }
  if (d[1L] < min_nlon) {
    stop(sprintf("First dimension (nlon) of '%s' must be >= %d.", name, min_nlon))
  }
  if (d[2L] < min_nlat) {
    stop(sprintf("Second dimension (nlat) of '%s' must be >= %d.", name, min_nlat))
  }
  if (anyNA(x)) {
    stop(sprintf("'%s' must not contain NA/NaN values.", name))
  }
  storage.mode(x) <- "double"
  nt <- if (length(d) == 2L) 1L else as.integer(prod(d[-c(1L, 2L)]))
  list(data = x, dim = d, nlon = d[1L], nlat = d[2L], nt = nt,
       extra_dim = if (length(d) > 2L) d[-c(1L, 2L)] else integer(0))
}

.check_grid_pair <- function(u, v, uname = "u", vname = "v", min_nlon = 4L) {
  iu <- .check_grid_array(u, uname, min_nlon = min_nlon)
  iv <- .check_grid_array(v, vname, min_nlon = min_nlon)
  if (!identical(iu$dim, iv$dim)) {
    stop(sprintf("'%s' and '%s' must have identical dimensions.", uname, vname))
  }
  list(u = iu$data, v = iv$data, dim = iu$dim, nlon = iu$nlon,
       nlat = iu$nlat, nt = iu$nt, extra_dim = iu$extra_dim)
}

.format_single <- function(vec, dim_vec) {
  dim(vec) <- dim_vec
  vec
}

.combine_Trailing <- function(list_vecs, base_dim) {
  k <- length(list_vecs)
  out <- unlist(list_vecs, use.names = FALSE)
  dim(out) <- c(base_dim, k)
  out
}

.copy_wrap_metadata <- function(res, template, long_name = NULL, units = NULL,
                                comp_names = NULL) {
  attrs <- attributes(template)
  if (!is.null(attrs)) {
    keep <- setdiff(names(attrs), c("dim", "dimnames", "names", "class"))
    for (nm in keep) {
      attr(res, nm) <- attrs[[nm]]
    }
  }
  tdn <- dimnames(template)
  if (!is.null(tdn)) {
    if (is.null(comp_names)) {
      dimnames(res) <- tdn
    } else if (length(dim(res)) == length(tdn) + 1L) {
      dimnames(res) <- c(tdn, list(component = comp_names))
    }
  } else if (!is.null(comp_names)) {
    dn <- vector("list", length(dim(res)))
    dn[[length(dn)]] <- comp_names
    dimnames(res) <- dn
  }
  if (!is.null(long_name)) {
    attr(res, "long_name") <- long_name
  }
  if (!is.null(units)) {
    attr(res, "units") <- units
  }
  res
}

.resolve_nlon <- function(nlon, mmax, nlat, is_vector = FALSE) {
  if (is.null(nlon)) {
    if (is_vector) {
      return(if (mmax == nlat) 2L * nlat else 2L * mmax)
    } else {
      return(if (mmax == nlat) 2L * nlat else 2L * (mmax - 1L))
    }
  }
  if (length(nlon) > 1L && !is.null(dim(nlon))) {
    return(as.integer(dim(nlon)[1L]))
  }
  nlon <- as.integer(nlon)[1L]
  if (is.na(nlon) || nlon < 4L) {
    stop("'nlon' must be an integer >= 4.")
  }
  nlon
}

#' Scalar Spherical Harmonic Analysis and Synthesis
#'
#' Functions for computing scalar spherical harmonic analysis (`shaec`, `shaeC`,
#' `shagc`, `shagC`) and synthesis (`shsec`, `shseC`, `shsgc`, `shsgC`,
#' `shsgc_R42`, `shsgc_R42_Wrap`) on global fixed (equally spaced) and Gaussian
#' grids.
#'
#' Spatial arrays `g` use column-major dimension ordering `(nlon, nlat, ...)`
#' with latitudes ordered South-to-North and no cyclic longitude point.
#' Spherical harmonic coefficient arrays `a` and `b` have dimension
#' `(N, nlat, ...)` where `N = min(nlat, (nlon + 2L) %/% 2L)`.
#' Uppercase variants (`shaeC`, `shagC`) return a combined array `ab` of
#' dimension `(N, nlat, ..., 2)`.
#'
#' @param g Numeric array of dimension `(nlon, nlat, ...)` containing the
#'   discrete scalar field to analyze.
#' @param a,b Numeric arrays of dimension `(N, nlat, ...)` containing real (`a`)
#'   and imaginary (`b`) spherical harmonic coefficients (or `(43, 43, ...)` for
#'   `shsgc_R42` and `shsgc_R42_Wrap`).
#' @param ab Numeric array of dimension `(N, nlat, ..., 2)` containing combined
#'   spherical harmonic coefficients as returned by `shaeC` or `shagC`.
#' @param nlon Optional integer specifying the longitude dimension of the
#'   synthesized spatial grid (or a template array whose first dimension is
#'   `nlon`). If `NULL`, `nlon` is inferred as `2 * (N - 1)` (or `2 * nlat` when
#'   `N == nlat`).
#' @return
#' \itemize{
#'   \item `shaec(g)`, `shagc(g)`: Named list `list(a = a, b = b)` where `a` and
#'     `b` have dimension `(N, nlat, ...)`.
#'   \item `shaeC(g)`, `shagC(g)`: Numeric array `ab` of dimension
#'     `(N, nlat, ..., 2)`.
#'   \item `shsec(a, b, nlon)`, `shsgc(a, b, nlon)`, `shseC(ab, nlon)`,
#'     `shsgC(ab, nlon)`: Synthesized spatial array of dimension
#'     `(nlon, nlat, ...)`.
#'   \item `shsgc_R42(a, b)`, `shsgc_R42_Wrap(a, b)`: Synthesized spatial array
#'     of dimension `(128, 108, ...)`.
#' }
#' @examples
#' nlon <- 36
#' nlat <- 19
#' lon <- seq(0, 360 - 360 / nlon, length.out = nlon) * pi / 180
#' lat <- seq(-90, 90, length.out = nlat) * pi / 180
#' g <- outer(cos(lon), cos(lat)^2)
#'
#' coefs <- shaec(g)
#' g_rec <- shsec(coefs$a, coefs$b, nlon)
#' max(abs(g - g_rec))
#'
#' ab <- shaeC(g)
#' g_rec2 <- shseC(ab, nlon)
#' max(abs(g - g_rec2))
#' @export
shaec <- function(g) {
  ig <- .check_grid_array(g, "g")
  mmax <- min(ig$nlat, (ig$nlon + 2L) %/% 2L)
  res <- .Call("c_sha", ig$data, ig$nlon, ig$nlat, ig$nt, 0L, PACKAGE = "spherepack")
  spec_dim <- c(mmax, ig$nlat, ig$extra_dim)
  list(a = .format_single(res$a, spec_dim),
       b = .format_single(res$b, spec_dim))
}

#' @rdname shaec
#' @export
shaeC <- function(g) {
  res <- shaec(g)
  .combine_Trailing(list(res$a, res$b), dim(res$a))
}

#' @rdname shaec
#' @export
shagc <- function(g) {
  ig <- .check_grid_array(g, "g")
  mmax <- min(ig$nlat, (ig$nlon + 2L) %/% 2L)
  res <- .Call("c_sha", ig$data, ig$nlon, ig$nlat, ig$nt, 1L, PACKAGE = "spherepack")
  spec_dim <- c(mmax, ig$nlat, ig$extra_dim)
  list(a = .format_single(res$a, spec_dim),
       b = .format_single(res$b, spec_dim))
}

#' @rdname shaec
#' @export
shagC <- function(g) {
  res <- shagc(g)
  .combine_Trailing(list(res$a, res$b), dim(res$a))
}

.shs_impl <- function(a, b, nlon, grid_type) {
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("'a' and 'b' must be numeric arrays.")
  }
  da <- dim(a)
  db <- dim(b)
  if (is.null(da) || length(da) < 2L || !identical(da, db)) {
    stop("'a' and 'b' must be arrays of identical dimensions (N, nlat, ...).")
  }
  mdab <- da[1L]
  nlat <- da[2L]
  nlon_val <- .resolve_nlon(nlon, mdab, nlat, is_vector = FALSE)
  mmax <- min(nlat, (nlon_val + 2L) %/% 2L)
  if (mdab < mmax) {
    stop(sprintf("First dimension of 'a' (%d) is smaller than required mmax (%d) for nlon = %d.",
                 mdab, mmax, nlon_val))
  }
  extra_dim <- if (length(da) > 2L) da[-c(1L, 2L)] else integer(0)
  nt <- if (length(da) == 2L) 1L else as.integer(prod(extra_dim))
  storage.mode(a) <- "double"
  storage.mode(b) <- "double"
  res <- .Call("c_shs", a, b, mdab, nlat, nlon_val, nt, as.integer(grid_type),
               PACKAGE = "spherepack")
  .format_single(res, c(nlon_val, nlat, extra_dim))
}

#' @rdname shaec
#' @export
shsec <- function(a, b, nlon = NULL) {
  .shs_impl(a, b, nlon, 0L)
}

#' @rdname shaec
#' @export
shsgc <- function(a, b, nlon = NULL) {
  .shs_impl(a, b, nlon, 1L)
}

.split_ab <- function(ab, name = "ab") {
  if (!is.numeric(ab)) {
    stop(sprintf("'%s' must be a numeric array.", name))
  }
  d <- dim(ab)
  if (is.null(d) || length(d) < 3L || d[length(d)] != 2L) {
    stop(sprintf("'%s' must be an array of dimension (N, nlat, ..., 2).", name))
  }
  base_dim <- d[-length(d)]
  half <- prod(base_dim)
  storage.mode(ab) <- "double"
  a <- ab[seq_len(half)]
  b <- ab[half + seq_len(half)]
  dim(a) <- base_dim
  dim(b) <- base_dim
  list(a = a, b = b)
}

#' @rdname shaec
#' @export
shseC <- function(ab, nlon = NULL) {
  sp <- .split_ab(ab, "ab")
  shsec(sp$a, sp$b, nlon = nlon)
}

#' @rdname shaec
#' @export
shsgC <- function(ab, nlon = NULL) {
  sp <- .split_ab(ab, "ab")
  shsgc(sp$a, sp$b, nlon = nlon)
}

#' @rdname shaec
#' @export
shsgc_R42 <- function(a, b) {
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("'a' and 'b' must be numeric arrays.")
  }
  da <- dim(a)
  db <- dim(b)
  if (is.null(da) || length(da) < 2L || !identical(da, db) ||
      da[1L] != 43L || da[2L] != 43L) {
    stop("'a' and 'b' must have first two dimensions (43, 43).")
  }
  extra_dim <- if (length(da) > 2L) da[-c(1L, 2L)] else integer(0)
  nt <- if (length(da) == 2L) 1L else as.integer(prod(extra_dim))
  storage.mode(a) <- "double"
  storage.mode(b) <- "double"
  dim(a) <- c(43L, 43L, nt)
  dim(b) <- c(43L, 43L, nt)

  nlat <- 108L
  nlon <- 128L
  mdab <- 65L
  a_full <- array(0.0, dim = c(mdab, nlat, nt))
  b_full <- array(0.0, dim = c(mdab, nlat, nt))

  # Detect whether (43, 43) is rhomboidal (m, n-m) or triangular (m, n)
  upper_idx <- outer(1:43, 1:43, `>`)
  lower_idx <- outer(2:43, 2:43, `<=`)
  has_upper <- any(a[, , 1L][upper_idx] != 0 | b[, , 1L][upper_idx] != 0)
  has_lower <- any(a[2:43, 2:43, 1L][lower_idx] != 0 | b[2:43, 2:43, 1L][lower_idx] != 0)
  is_triangular <- !has_upper && has_lower

  if (is_triangular) {
    a_full[1:43, 1:43, ] <- a
    b_full[1:43, 1:43, ] <- b
  } else {
    for (m_idx in 1:43) {
      for (j_idx in 1:43) {
        n_idx <- m_idx + j_idx - 1L
        a_full[m_idx, n_idx, ] <- a[m_idx, j_idx, ]
        b_full[m_idx, n_idx, ] <- b[m_idx, j_idx, ]
      }
    }
  }

  res <- .Call("c_shs", a_full, b_full, mdab, nlat, nlon, nt, 1L,
               PACKAGE = "spherepack")
  .format_single(res, c(nlon, nlat, extra_dim))
}

#' @rdname shaec
#' @export
shsgc_R42_Wrap <- function(a, b) {
  res <- shsgc_R42(a, b)
  attrs <- attributes(a)
  if (!is.null(attrs)) {
    keep <- setdiff(names(attrs), c("dim", "dimnames", "names", "class"))
    for (nm in keep) {
      attr(res, nm) <- attrs[[nm]]
    }
  }
  g <- gaqd(108L)
  lon <- seq(0.0, 360.0 - 360.0 / 128.0, length.out = 128L)
  attr(res, "lon") <- lon
  attr(res, "lat") <- g$lat
  if (is.null(attr(res, "long_name"))) {
    attr(res, "long_name") <- "R42 spherical harmonic synthesis"
  }
  res
}

#' Vector Spherical Harmonic Analysis and Synthesis
#'
#' Functions for computing vector spherical harmonic analysis (`vhaec`, `vhaeC`,
#' `vhagc`, `vhagC`) and synthesis (`vhsec`, `vhseC`, `vhsgc`, `vhsgC`) on global
#' fixed (equally spaced) and Gaussian grids.
#'
#' Spatial vector components `u` (zonal, eastward) and `v` (meridional,
#' northward) use column-major dimension ordering `(nlon, nlat, ...)` with
#' latitudes ordered South-to-North. Vector spherical harmonic coefficient
#' arrays `br`, `bi`, `cr`, `ci` have dimension `(N, nlat, ...)` where
#' `N = min(nlat, (nlon + 1L) %/% 2L)`. Uppercase analysis functions (`vhaeC`,
#' `vhagC`) return a combined array `bc` of dimension `(N, nlat, ..., 4)`, and
#' uppercase synthesis functions (`vhseC`, `vhsgC`) return a combined array of
#' dimension `(nlon, nlat, ..., 2)` containing `(u, v)`.
#'
#' @param u,v Numeric arrays of dimension `(nlon, nlat, ...)` giving zonal (`u`)
#'   and meridional (`v`) vector components.
#' @param br,bi,cr,ci Numeric arrays of dimension `(N, nlat, ...)` giving vector
#'   spherical harmonic coefficients in mathematical coordinates.
#' @param bc Numeric array of dimension `(N, nlat, ..., 4)` containing combined
#'   vector spherical harmonic coefficients `(br, bi, cr, ci)`.
#' @param nlon Optional integer specifying the longitude dimension of the
#'   synthesized spatial grid (or a template array).
#' @return
#' \itemize{
#'   \item `vhaec(u, v)`, `vhagc(u, v)`: Named list
#'     `list(br = br, bi = bi, cr = cr, ci = ci)`.
#'   \item `vhaeC(u, v)`, `vhagC(u, v)`: Combined array `bc` of dimension
#'     `(N, nlat, ..., 4)`.
#'   \item `vhsec(br, bi, cr, ci, nlon)`, `vhsgc(br, bi, cr, ci, nlon)`: Named
#'     list `list(u = u, v = v)`.
#'   \item `vhseC(bc, nlon)`, `vhsgC(bc, nlon)`: Combined array of dimension
#'     `(nlon, nlat, ..., 2)`.
#' }
#' @examples
#' nlon <- 36
#' nlat <- 19
#' lon <- seq(0, 360 - 360 / nlon, length.out = nlon) * pi / 180
#' lat <- seq(-90, 90, length.out = nlat) * pi / 180
#' u <- outer(cos(lon), cos(lat))
#' v <- outer(-sin(lon), sin(lat) * cos(lat))
#'
#' vc <- vhaec(u, v)
#' uv <- vhsec(vc$br, vc$bi, vc$cr, vc$ci, nlon)
#' max(abs(u - uv$u))
#' max(abs(v - uv$v))
#' @export
vhaec <- function(u, v) {
  iv <- .check_grid_pair(u, v, "u", "v", min_nlon = 1L)
  mmax <- min(iv$nlat, (iv$nlon + 1L) %/% 2L)
  res <- .Call("c_vha", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt, 0L, PACKAGE = "spherepack")
  spec_dim <- c(mmax, iv$nlat, iv$extra_dim)
  list(br = .format_single(res$br, spec_dim),
       bi = .format_single(res$bi, spec_dim),
       cr = .format_single(res$cr, spec_dim),
       ci = .format_single(res$ci, spec_dim))
}

#' @rdname vhaec
#' @export
vhaeC <- function(u, v) {
  res <- vhaec(u, v)
  .combine_Trailing(list(res$br, res$bi, res$cr, res$ci), dim(res$br))
}

#' @rdname vhaec
#' @export
vhagc <- function(u, v) {
  iv <- .check_grid_pair(u, v, "u", "v", min_nlon = 1L)
  mmax <- min(iv$nlat, (iv$nlon + 1L) %/% 2L)
  res <- .Call("c_vha", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt, 1L, PACKAGE = "spherepack")
  spec_dim <- c(mmax, iv$nlat, iv$extra_dim)
  list(br = .format_single(res$br, spec_dim),
       bi = .format_single(res$bi, spec_dim),
       cr = .format_single(res$cr, spec_dim),
       ci = .format_single(res$ci, spec_dim))
}

#' @rdname vhaec
#' @export
vhagC <- function(u, v) {
  res <- vhagc(u, v)
  .combine_Trailing(list(res$br, res$bi, res$cr, res$ci), dim(res$br))
}

.vhs_impl <- function(br, bi, cr, ci, nlon, grid_type) {
  if (!is.numeric(br) || !is.numeric(bi) || !is.numeric(cr) || !is.numeric(ci)) {
    stop("'br', 'bi', 'cr', and 'ci' must be numeric arrays.")
  }
  d <- dim(br)
  if (is.null(d) || length(d) < 2L ||
      !identical(d, dim(bi)) || !identical(d, dim(cr)) || !identical(d, dim(ci))) {
    stop("'br', 'bi', 'cr', and 'ci' must have identical dimensions (N, nlat, ...).")
  }
  mdab <- d[1L]
  nlat <- d[2L]
  nlon_val <- .resolve_nlon(nlon, mdab, nlat, is_vector = TRUE)
  mmax <- min(nlat, (nlon_val + 1L) %/% 2L)
  if (mdab < mmax) {
    stop(sprintf("First dimension of 'br' (%d) is smaller than required mmax (%d) for nlon = %d.",
                 mdab, mmax, nlon_val))
  }
  extra_dim <- if (length(d) > 2L) d[-c(1L, 2L)] else integer(0)
  nt <- if (length(d) == 2L) 1L else as.integer(prod(extra_dim))
  storage.mode(br) <- "double"
  storage.mode(bi) <- "double"
  storage.mode(cr) <- "double"
  storage.mode(ci) <- "double"
  res <- .Call("c_vhs", br, bi, cr, ci, mdab, nlat, nlon_val, nt, as.integer(grid_type),
               PACKAGE = "spherepack")
  out_dim <- c(nlon_val, nlat, extra_dim)
  list(u = .format_single(res$u, out_dim),
       v = .format_single(res$v, out_dim))
}

#' @rdname vhaec
#' @export
vhsec <- function(br, bi, cr, ci, nlon = NULL) {
  .vhs_impl(br, bi, cr, ci, nlon, 0L)
}

#' @rdname vhaec
#' @export
vhsgc <- function(br, bi, cr, ci, nlon = NULL) {
  .vhs_impl(br, bi, cr, ci, nlon, 1L)
}

.split_bc <- function(bc, name = "bc") {
  if (!is.numeric(bc)) {
    stop(sprintf("'%s' must be a numeric array.", name))
  }
  d <- dim(bc)
  if (is.null(d) || length(d) < 3L || d[length(d)] != 4L) {
    stop(sprintf("'%s' must be an array of dimension (N, nlat, ..., 4).", name))
  }
  base_dim <- d[-length(d)]
  q <- prod(base_dim)
  storage.mode(bc) <- "double"
  br <- bc[seq_len(q)]
  bi <- bc[q + seq_len(q)]
  cr <- bc[2L * q + seq_len(q)]
  ci <- bc[3L * q + seq_len(q)]
  dim(br) <- base_dim
  dim(bi) <- base_dim
  dim(cr) <- base_dim
  dim(ci) <- base_dim
  list(br = br, bi = bi, cr = cr, ci = ci)
}

#' @rdname vhaec
#' @export
vhseC <- function(bc, nlon = NULL) {
  sp <- .split_bc(bc, "bc")
  res <- vhsec(sp$br, sp$bi, sp$cr, sp$ci, nlon = nlon)
  .combine_Trailing(list(res$u, res$v), dim(res$u))
}

#' @rdname vhaec
#' @export
vhsgC <- function(bc, nlon = NULL) {
  sp <- .split_bc(bc, "bc")
  res <- vhsgc(sp$br, sp$bi, sp$cr, sp$ci, nlon = nlon)
  .combine_Trailing(list(res$u, res$v), dim(res$u))
}

#' Divergence and Vorticity from Wind Components (and Inverse Operators)
#'
#' Computes divergence (`uv2dv*`), relative vorticity (`uv2vr*`), or both
#' simultaneously (`uv2vrdv*`) from wind components `(u, v)`, as well as the
#' inverse operators computing divergent wind from divergence (`dv2uv*`),
#' rotational (non-divergent) wind from vorticity (`vr2uv*`), and total wind
#' from vorticity and divergence (`vrdv2uv*`) on fixed (`*f`, `*F`) and
#' Gaussian (`*g`, `*G`) global grids.
#'
#' @param u,v Numeric arrays of dimension `(nlon, nlat, ...)` containing zonal
#'   (`u`) and meridional (`v`) wind components in ascending latitude order
#'   (South-to-North).
#' @param dv Numeric array of dimension `(nlon, nlat, ...)` containing
#'   divergence values.
#' @param vr Numeric array of dimension `(nlon, nlat, ...)` containing relative
#'   vorticity values.
#' @return
#' \itemize{
#'   \item `uv2dvf`, `uv2dvg`, `uv2dvF`, `uv2dvG`, `uv2dvF_Wrap`, `uv2dvG_Wrap`:
#'     Divergence array of the same dimension as `u`.
#'   \item `uv2vrf`, `uv2vrg`, `uv2vrF`, `uv2vrG`, `uv2vrF_Wrap`, `uv2vrG_Wrap`:
#'     Vorticity array of the same dimension as `u`.
#'   \item `uv2vrdvf`, `uv2vrdvg`: Named list `list(vr = vr, dv = dv)`.
#'   \item `uv2vrdvF`, `uv2vrdvG`: Combined array of dimension
#'     `(nlon, nlat, ..., 2)` containing `(vr, dv)`.
#'   \item `dv2uvf`, `dv2uvg`, `vr2uvf`, `vr2uvg`, `vrdv2uvf`, `vrdv2uvg`:
#'     Named list `list(u = u, v = v)`.
#'   \item `dv2uvF`, `dv2uvG`, `dv2uvF_Wrap`, `dv2uvG_Wrap`, `vr2uvF`, `vr2uvG`,
#'     `vr2uvF_Wrap`, `vr2uvG_Wrap`, `vrdv2uvF`, `vrdv2uvG`: Combined array of
#'     dimension `(nlon, nlat, ..., 2)` containing `(u, v)`.
#' }
#' @examples
#' nlon <- 36
#' nlat <- 19
#' lon <- seq(0, 360 - 360 / nlon, length.out = nlon) * pi / 180
#' lat <- seq(-90, 90, length.out = nlat) * pi / 180
#' u <- outer(cos(lon), cos(lat))
#' v <- outer(-sin(lon), sin(lat) * cos(lat))
#'
#' vrdv <- uv2vrdvf(u, v)
#' uv_rec <- vrdv2uvf(vrdv$vr, vrdv$dv)
#' max(abs(u - uv_rec$u))
#' max(abs(v - uv_rec$v))
#' @export
uv2dvf <- function(u, v) {
  iv <- .check_grid_pair(u, v)
  res <- .Call("c_uv2vrdv_op", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               0L, 1L, get_sphere_radius(), PACKAGE = "spherepack")
  .format_single(res, iv$dim)
}

#' @rdname uv2dvf
#' @export
uv2dvF <- function(u, v) {
  uv2dvf(u, v)
}

#' @rdname uv2dvf
#' @export
uv2dvF_Wrap <- function(u, v) {
  res <- uv2dvf(u, v)
  .copy_wrap_metadata(res, u, long_name = "divergence", units = "1/s")
}

#' @rdname uv2dvf
#' @export
uv2dvg <- function(u, v) {
  iv <- .check_grid_pair(u, v)
  res <- .Call("c_uv2vrdv_op", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               1L, 1L, get_sphere_radius(), PACKAGE = "spherepack")
  .format_single(res, iv$dim)
}

#' @rdname uv2dvf
#' @export
uv2dvG <- function(u, v) {
  uv2dvg(u, v)
}

#' @rdname uv2dvf
#' @export
uv2dvG_Wrap <- function(u, v) {
  res <- uv2dvg(u, v)
  .copy_wrap_metadata(res, u, long_name = "divergence", units = "1/s")
}

#' @rdname uv2dvf
#' @export
uv2vrf <- function(u, v) {
  iv <- .check_grid_pair(u, v)
  res <- .Call("c_uv2vrdv_op", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               0L, 2L, get_sphere_radius(), PACKAGE = "spherepack")
  .format_single(res, iv$dim)
}

#' @rdname uv2dvf
#' @export
uv2vrF <- function(u, v) {
  uv2vrf(u, v)
}

#' @rdname uv2dvf
#' @export
uv2vrF_Wrap <- function(u, v) {
  res <- uv2vrf(u, v)
  .copy_wrap_metadata(res, u, long_name = "vorticity", units = "1/s")
}

#' @rdname uv2dvf
#' @export
uv2vrg <- function(u, v) {
  iv <- .check_grid_pair(u, v)
  res <- .Call("c_uv2vrdv_op", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               1L, 2L, get_sphere_radius(), PACKAGE = "spherepack")
  .format_single(res, iv$dim)
}

#' @rdname uv2dvf
#' @export
uv2vrG <- function(u, v) {
  uv2vrg(u, v)
}

#' @rdname uv2dvf
#' @export
uv2vrG_Wrap <- function(u, v) {
  res <- uv2vrg(u, v)
  .copy_wrap_metadata(res, u, long_name = "vorticity", units = "1/s")
}

#' @rdname uv2dvf
#' @export
uv2vrdvf <- function(u, v) {
  iv <- .check_grid_pair(u, v)
  res <- .Call("c_uv2vrdv_op", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               0L, 3L, get_sphere_radius(), PACKAGE = "spherepack")
  list(vr = .format_single(res$vr, iv$dim),
       dv = .format_single(res$dv, iv$dim))
}

#' @rdname uv2dvf
#' @export
uv2vrdvF <- function(u, v) {
  res <- uv2vrdvf(u, v)
  .combine_Trailing(list(res$vr, res$dv), dim(res$vr))
}

#' @rdname uv2dvf
#' @export
uv2vrdvg <- function(u, v) {
  iv <- .check_grid_pair(u, v)
  res <- .Call("c_uv2vrdv_op", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               1L, 3L, get_sphere_radius(), PACKAGE = "spherepack")
  list(vr = .format_single(res$vr, iv$dim),
       dv = .format_single(res$dv, iv$dim))
}

#' @rdname uv2dvf
#' @export
uv2vrdvG <- function(u, v) {
  res <- uv2vrdvg(u, v)
  .combine_Trailing(list(res$vr, res$dv), dim(res$vr))
}

#' @rdname uv2dvf
#' @export
dv2uvf <- function(dv) {
  idv <- .check_grid_array(dv, "dv")
  res <- .Call("c_vrdv2uv_op", NULL, idv$data, idv$nlon, idv$nlat, idv$nt,
               0L, 1L, get_sphere_radius(), PACKAGE = "spherepack")
  list(u = .format_single(res$u, idv$dim),
       v = .format_single(res$v, idv$dim))
}

#' @rdname uv2dvf
#' @export
dv2uvF <- function(dv) {
  res <- dv2uvf(dv)
  .combine_Trailing(list(res$u, res$v), dim(res$u))
}

#' @rdname uv2dvf
#' @export
dv2uvF_Wrap <- function(dv) {
  res <- dv2uvF(dv)
  .copy_wrap_metadata(res, dv, long_name = "divergent wind components",
                      units = "m/s", comp_names = c("u", "v"))
}

#' @rdname uv2dvf
#' @export
dv2uvg <- function(dv) {
  idv <- .check_grid_array(dv, "dv")
  res <- .Call("c_vrdv2uv_op", NULL, idv$data, idv$nlon, idv$nlat, idv$nt,
               1L, 1L, get_sphere_radius(), PACKAGE = "spherepack")
  list(u = .format_single(res$u, idv$dim),
       v = .format_single(res$v, idv$dim))
}

#' @rdname uv2dvf
#' @export
dv2uvG <- function(dv) {
  res <- dv2uvg(dv)
  .combine_Trailing(list(res$u, res$v), dim(res$u))
}

#' @rdname uv2dvf
#' @export
dv2uvG_Wrap <- function(dv) {
  res <- dv2uvG(dv)
  .copy_wrap_metadata(res, dv, long_name = "divergent wind components",
                      units = "m/s", comp_names = c("u", "v"))
}

#' @rdname uv2dvf
#' @export
vr2uvf <- function(vr) {
  ivr <- .check_grid_array(vr, "vr")
  res <- .Call("c_vrdv2uv_op", ivr$data, NULL, ivr$nlon, ivr$nlat, ivr$nt,
               0L, 2L, get_sphere_radius(), PACKAGE = "spherepack")
  list(u = .format_single(res$u, ivr$dim),
       v = .format_single(res$v, ivr$dim))
}

#' @rdname uv2dvf
#' @export
vr2uvF <- function(vr) {
  res <- vr2uvf(vr)
  .combine_Trailing(list(res$u, res$v), dim(res$u))
}

#' @rdname uv2dvf
#' @export
vr2uvF_Wrap <- function(vr) {
  res <- vr2uvF(vr)
  .copy_wrap_metadata(res, vr, long_name = "rotational wind components",
                      units = "m/s", comp_names = c("u", "v"))
}

#' @rdname uv2dvf
#' @export
vr2uvg <- function(vr) {
  ivr <- .check_grid_array(vr, "vr")
  res <- .Call("c_vrdv2uv_op", ivr$data, NULL, ivr$nlon, ivr$nlat, ivr$nt,
               1L, 2L, get_sphere_radius(), PACKAGE = "spherepack")
  list(u = .format_single(res$u, ivr$dim),
       v = .format_single(res$v, ivr$dim))
}

#' @rdname uv2dvf
#' @export
vr2uvG <- function(vr) {
  res <- vr2uvg(vr)
  .combine_Trailing(list(res$u, res$v), dim(res$u))
}

#' @rdname uv2dvf
#' @export
vr2uvG_Wrap <- function(vr) {
  res <- vr2uvG(vr)
  .copy_wrap_metadata(res, vr, long_name = "rotational wind components",
                      units = "m/s", comp_names = c("u", "v"))
}

#' @rdname uv2dvf
#' @export
vrdv2uvf <- function(vr, dv) {
  ivd <- .check_grid_pair(vr, dv, "vr", "dv")
  res <- .Call("c_vrdv2uv_op", ivd$u, ivd$v, ivd$nlon, ivd$nlat, ivd$nt,
               0L, 3L, get_sphere_radius(), PACKAGE = "spherepack")
  list(u = .format_single(res$u, ivd$dim),
       v = .format_single(res$v, ivd$dim))
}

#' @rdname uv2dvf
#' @export
vrdv2uvF <- function(vr, dv) {
  res <- vrdv2uvf(vr, dv)
  .combine_Trailing(list(res$u, res$v), dim(res$u))
}

#' @rdname uv2dvf
#' @export
vrdv2uvg <- function(vr, dv) {
  ivd <- .check_grid_pair(vr, dv, "vr", "dv")
  res <- .Call("c_vrdv2uv_op", ivd$u, ivd$v, ivd$nlon, ivd$nlat, ivd$nt,
               1L, 3L, get_sphere_radius(), PACKAGE = "spherepack")
  list(u = .format_single(res$u, ivd$dim),
       v = .format_single(res$v, ivd$dim))
}

#' @rdname uv2dvf
#' @export
vrdv2uvG <- function(vr, dv) {
  res <- vrdv2uvg(vr, dv)
  .combine_Trailing(list(res$u, res$v), dim(res$u))
}

#' Streamfunction and Velocity Potential from Wind Components (and Inverse)
#'
#' `uv2sfvpf`, `uv2sfvpF`, `uv2sfvpg`, and `uv2sfvpG` compute the streamfunction
#' (`sf`) and velocity potential (`vp`) from wind components `(u, v)` on fixed
#' or Gaussian global grids. `sfvp2uvf` and `sfvp2uvg` perform the inverse
#' operation, reconstructing `(u, v)` from `(sf, vp)`.
#'
#' @param u,v Numeric arrays of dimension `(nlon, nlat, ...)` containing zonal
#'   (`u`) and meridional (`v`) wind components.
#' @param sf,vp Numeric arrays of dimension `(nlon, nlat, ...)` containing
#'   streamfunction (`sf`) and velocity potential (`vp`).
#' @return
#' \itemize{
#'   \item `uv2sfvpf`, `uv2sfvpg`: Named list `list(sf = sf, vp = vp)`.
#'   \item `uv2sfvpF`, `uv2sfvpG`: Combined array of dimension
#'     `(nlon, nlat, ..., 2)` containing `(sf, vp)`.
#'   \item `sfvp2uvf`, `sfvp2uvg`: Named list `list(u = u, v = v)`.
#' }
#' @examples
#' nlon <- 36
#' nlat <- 19
#' lon <- seq(0, 360 - 360 / nlon, length.out = nlon) * pi / 180
#' lat <- seq(-90, 90, length.out = nlat) * pi / 180
#' u <- outer(cos(lon), cos(lat))
#' v <- outer(-sin(lon), sin(lat) * cos(lat))
#'
#' sv <- uv2sfvpf(u, v)
#' uv_rec <- sfvp2uvf(sv$sf, sv$vp)
#' max(abs(u - uv_rec$u))
#' max(abs(v - uv_rec$v))
#' @export
uv2sfvpf <- function(u, v) {
  iv <- .check_grid_pair(u, v)
  res <- .Call("c_uv2sfvp", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               0L, get_sphere_radius(), PACKAGE = "spherepack")
  list(sf = .format_single(res$sf, iv$dim),
       vp = .format_single(res$vp, iv$dim))
}

#' @rdname uv2sfvpf
#' @export
uv2sfvpF <- function(u, v) {
  res <- uv2sfvpf(u, v)
  .combine_Trailing(list(res$sf, res$vp), dim(res$sf))
}

#' @rdname uv2sfvpf
#' @export
uv2sfvpg <- function(u, v) {
  iv <- .check_grid_pair(u, v)
  res <- .Call("c_uv2sfvp", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               1L, get_sphere_radius(), PACKAGE = "spherepack")
  list(sf = .format_single(res$sf, iv$dim),
       vp = .format_single(res$vp, iv$dim))
}

#' @rdname uv2sfvpf
#' @export
uv2sfvpG <- function(u, v) {
  res <- uv2sfvpg(u, v)
  .combine_Trailing(list(res$sf, res$vp), dim(res$sf))
}

#' @rdname uv2sfvpf
#' @export
sfvp2uvf <- function(sf, vp) {
  isv <- .check_grid_pair(sf, vp, "sf", "vp")
  res <- .Call("c_sfvp2uv", isv$u, isv$v, isv$nlon, isv$nlat, isv$nt,
               0L, get_sphere_radius(), PACKAGE = "spherepack")
  list(u = .format_single(res$u, isv$dim),
       v = .format_single(res$v, isv$dim))
}

#' @rdname uv2sfvpf
#' @export
sfvp2uvg <- function(sf, vp) {
  isv <- .check_grid_pair(sf, vp, "sf", "vp")
  res <- .Call("c_sfvp2uv", isv$u, isv$v, isv$nlon, isv$nlat, isv$nt,
               1L, get_sphere_radius(), PACKAGE = "spherepack")
  list(u = .format_single(res$u, isv$dim),
       v = .format_single(res$v, isv$dim))
}

#' Spherical Gradient, Inverse Gradient, and Latitudinal Vector Derivative
#'
#' `gradsf` and `gradsg` compute the zonal (`gzx`) and meridional (`gzy`)
#' components of the horizontal gradient of a scalar field `z` on fixed and
#' Gaussian grids. `igradsf`, `igradsg`, `igradsF`, and `igradsG` invert the
#' gradient to recover the zero-mean scalar field `z`. `lderuvf` and `lderuvg`
#' compute the northward latitudinal derivatives `(uy, vy)` of a vector field
#' `(u, v)`.
#'
#' @param z Numeric array of dimension `(nlon, nlat, ...)` containing a scalar
#'   field on a global grid.
#' @param gzx,gzy Numeric arrays of dimension `(nlon, nlat, ...)` containing the
#'   zonal (eastward) and meridional (northward) gradient components.
#' @param u,v Numeric arrays of dimension `(nlon, nlat, ...)` containing zonal
#'   and meridional vector components.
#' @return
#' \itemize{
#'   \item `gradsf`, `gradsg`: Named list `list(gzx = gzx, gzy = gzy)`.
#'   \item `igradsf`, `igradsg`, `igradsF`, `igradsG`: Scalar array `z` of
#'     dimension `(nlon, nlat, ...)`.
#'   \item `lderuvf`, `lderuvg`: Named list `list(uy = uy, vy = vy)`.
#' }
#' @examples
#' nlon <- 36
#' nlat <- 19
#' lon <- seq(0, 360 - 360 / nlon, length.out = nlon) * pi / 180
#' lat <- seq(-90, 90, length.out = nlat) * pi / 180
#' z <- outer(cos(lon), cos(lat))
#'
#' gr <- gradsf(z)
#' z_rec <- igradsF(gr$gzx, gr$gzy)
#' max(abs(z - z_rec))
#' @export
gradsf <- function(z) {
  iz <- .check_grid_array(z, "z")
  res <- .Call("c_grads", iz$data, iz$nlon, iz$nlat, iz$nt,
               0L, get_sphere_radius(), PACKAGE = "spherepack")
  list(gzx = .format_single(res$gzx, iz$dim),
       gzy = .format_single(res$gzy, iz$dim))
}

#' @rdname gradsf
#' @export
gradsg <- function(z) {
  iz <- .check_grid_array(z, "z")
  res <- .Call("c_grads", iz$data, iz$nlon, iz$nlat, iz$nt,
               1L, get_sphere_radius(), PACKAGE = "spherepack")
  list(gzx = .format_single(res$gzx, iz$dim),
       gzy = .format_single(res$gzy, iz$dim))
}

#' @rdname gradsf
#' @export
igradsf <- function(gzx, gzy) {
  ig <- .check_grid_pair(gzx, gzy, "gzx", "gzy")
  res <- .Call("c_igrads", ig$u, ig$v, ig$nlon, ig$nlat, ig$nt,
               0L, get_sphere_radius(), PACKAGE = "spherepack")
  .format_single(res, ig$dim)
}

#' @rdname gradsf
#' @export
igradsF <- function(gzx, gzy) {
  igradsf(gzx, gzy)
}

#' @rdname gradsf
#' @export
igradsg <- function(gzx, gzy) {
  ig <- .check_grid_pair(gzx, gzy, "gzx", "gzy")
  res <- .Call("c_igrads", ig$u, ig$v, ig$nlon, ig$nlat, ig$nt,
               1L, get_sphere_radius(), PACKAGE = "spherepack")
  .format_single(res, ig$dim)
}

#' @rdname gradsf
#' @export
igradsG <- function(gzx, gzy) {
  igradsg(gzx, gzy)
}

#' @rdname gradsf
#' @export
lderuvf <- function(u, v) {
  iv <- .check_grid_pair(u, v, min_nlon = 1L)
  res <- .Call("c_lderuv", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               0L, get_sphere_radius(), PACKAGE = "spherepack")
  list(uy = .format_single(res$uy, iv$dim),
       vy = .format_single(res$vy, iv$dim))
}

#' @rdname gradsf
#' @export
lderuvg <- function(u, v) {
  iv <- .check_grid_pair(u, v, min_nlon = 1L)
  res <- .Call("c_lderuv", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               1L, get_sphere_radius(), PACKAGE = "spherepack")
  list(uy = .format_single(res$uy, iv$dim),
       vy = .format_single(res$vy, iv$dim))
}

#' Scalar and Vector Laplacian and Inverse Laplacian Operators
#'
#' Computes the scalar Laplacian (`lapsf`, `lapsF`, `lapsg`, `lapsG`), inverse
#' scalar Laplacian / Helmholtz solver (`ilapsf`, `ilapsF`, `ilapsF_Wrap`,
#' `ilapsg`, `ilapsG`, `ilapsG_Wrap`), vector Laplacian (`lapvf`, `lapvg`), and
#' inverse vector Laplacian (`ilapvf`, `ilapvg`) on fixed and Gaussian global
#' grids.
#'
#' @param z Numeric array of dimension `(nlon, nlat, ...)` containing a scalar
#'   field.
#' @param zlap Numeric array of dimension `(nlon, nlat, ...)` containing the
#'   scalar Laplacian field to invert.
#' @param zlmbda Numeric scalar or array of Helmholtz constants (default `0.0`
#'   solves the Poisson equation). If `zlap` has extra dimensions beyond
#'   `(nlon, nlat)`, `zlmbda` may be a scalar or an array matching those extra
#'   dimensions.
#' @param u,v Numeric arrays of dimension `(nlon, nlat, ...)` containing vector
#'   components.
#' @param ulap,vlap Numeric arrays of dimension `(nlon, nlat, ...)` containing
#'   vector Laplacian components to invert.
#' @return
#' \itemize{
#'   \item `lapsf`, `lapsF`, `lapsg`, `lapsG`: Scalar Laplacian array of the
#'     same dimension as `z`.
#'   \item `ilapsf`, `ilapsF`, `ilapsF_Wrap`, `ilapsg`, `ilapsG`, `ilapsG_Wrap`:
#'     Inverted scalar array of the same dimension as `zlap`.
#'   \item `lapvf`, `lapvg`: Named list `list(ulap = ulap, vlap = vlap)`.
#'   \item `ilapvf`, `ilapvg`: Named list `list(u = u, v = v)`.
#' }
#' @examples
#' nlon <- 36
#' nlat <- 19
#' lon <- seq(0, 360 - 360 / nlon, length.out = nlon) * pi / 180
#' lat <- seq(-90, 90, length.out = nlat) * pi / 180
#' z <- outer(cos(lon), cos(lat))
#'
#' zlap <- lapsF(z)
#' z_rec <- ilapsF(zlap, 0)
#' max(abs(z - z_rec))
#' @export
lapsf <- function(z) {
  iz <- .check_grid_array(z, "z")
  res <- .Call("c_laps", iz$data, iz$nlon, iz$nlat, iz$nt,
               0L, get_sphere_radius(), PACKAGE = "spherepack")
  .format_single(res, iz$dim)
}

#' @rdname lapsf
#' @export
lapsF <- function(z) {
  lapsf(z)
}

#' @rdname lapsf
#' @export
lapsg <- function(z) {
  iz <- .check_grid_array(z, "z")
  res <- .Call("c_laps", iz$data, iz$nlon, iz$nlat, iz$nt,
               1L, get_sphere_radius(), PACKAGE = "spherepack")
  .format_single(res, iz$dim)
}

#' @rdname lapsf
#' @export
lapsG <- function(z) {
  lapsg(z)
}

.prepare_zlmbda <- function(zlmbda, nt) {
  if (!is.numeric(zlmbda) || anyNA(zlmbda)) {
    stop("'zlmbda' must be numeric without NA values.")
  }
  z <- as.double(zlmbda)
  if (length(z) == 1L) {
    rep(z, nt)
  } else if (length(z) == nt) {
    z
  } else {
    stop(sprintf("'zlmbda' must have length 1 or %d (matching extra dimensions of 'zlap').", nt))
  }
}

#' @rdname lapsf
#' @export
ilapsf <- function(zlap, zlmbda = 0.0) {
  iz <- .check_grid_array(zlap, "zlap")
  zlm <- .prepare_zlmbda(zlmbda, iz$nt)
  res <- .Call("c_ilaps", iz$data, zlm, iz$nlon, iz$nlat, iz$nt,
               0L, get_sphere_radius(), PACKAGE = "spherepack")
  .format_single(res, iz$dim)
}

#' @rdname lapsf
#' @export
ilapsF <- function(zlap, zlmbda = 0.0) {
  ilapsf(zlap, zlmbda)
}

#' @rdname lapsf
#' @export
ilapsF_Wrap <- function(zlap, zlmbda = 0.0) {
  res <- ilapsf(zlap, zlmbda)
  .copy_wrap_metadata(res, zlap, long_name = "inverse Laplacian")
}

#' @rdname lapsf
#' @export
ilapsg <- function(zlap, zlmbda = 0.0) {
  iz <- .check_grid_array(zlap, "zlap")
  zlm <- .prepare_zlmbda(zlmbda, iz$nt)
  res <- .Call("c_ilaps", iz$data, zlm, iz$nlon, iz$nlat, iz$nt,
               1L, get_sphere_radius(), PACKAGE = "spherepack")
  .format_single(res, iz$dim)
}

#' @rdname lapsf
#' @export
ilapsG <- function(zlap, zlmbda = 0.0) {
  ilapsg(zlap, zlmbda)
}

#' @rdname lapsf
#' @export
ilapsG_Wrap <- function(zlap, zlmbda = 0.0) {
  res <- ilapsg(zlap, zlmbda)
  .copy_wrap_metadata(res, zlap, long_name = "inverse Laplacian")
}

#' @rdname lapsf
#' @export
lapvf <- function(u, v) {
  iv <- .check_grid_pair(u, v, min_nlon = 1L)
  res <- .Call("c_vlap_op", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               0L, 0L, get_sphere_radius(), PACKAGE = "spherepack")
  list(ulap = .format_single(res$ulap, iv$dim),
       vlap = .format_single(res$vlap, iv$dim))
}

#' @rdname lapsf
#' @export
lapvg <- function(u, v) {
  iv <- .check_grid_pair(u, v, min_nlon = 1L)
  res <- .Call("c_vlap_op", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               1L, 0L, get_sphere_radius(), PACKAGE = "spherepack")
  list(ulap = .format_single(res$ulap, iv$dim),
       vlap = .format_single(res$vlap, iv$dim))
}

#' @rdname lapsf
#' @export
ilapvf <- function(ulap, vlap) {
  iv <- .check_grid_pair(ulap, vlap, "ulap", "vlap", min_nlon = 1L)
  res <- .Call("c_vlap_op", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               0L, 1L, get_sphere_radius(), PACKAGE = "spherepack")
  list(u = .format_single(res$u, iv$dim),
       v = .format_single(res$v, iv$dim))
}

#' @rdname lapsf
#' @export
ilapvg <- function(ulap, vlap) {
  iv <- .check_grid_pair(ulap, vlap, "ulap", "vlap", min_nlon = 1L)
  res <- .Call("c_vlap_op", iv$u, iv$v, iv$nlon, iv$nlat, iv$nt,
               1L, 1L, get_sphere_radius(), PACKAGE = "spherepack")
  list(u = .format_single(res$u, iv$dim),
       v = .format_single(res$v, iv$dim))
}

#' Spectral Truncation and Exponential Tapering of Spherical Harmonics
#'
#' Performs triangular truncation (`tri_trunc`, `tri_trunC`), rhomboidal
#' truncation (`rhomb_trunc`, `rhomb_trunC`), or isotropic exponential tapering
#' (`exp_tapersh_wgts`, `exp_tapersh`, `exp_tapershC`) of spherical harmonic
#' coefficients following Sardeshmukh and Hoskins (1984).
#'
#' @param a,b Numeric arrays of dimension `(N_m, nlat, ...)` containing real and
#'   imaginary spherical harmonic coefficients.
#' @param ab Numeric array of dimension `(N_m, nlat, ..., 2)` containing
#'   combined spherical harmonic coefficients.
#' @param T Non-negative integer truncation wavenumber.
#' @param M Integer specifying the number of total-wavenumber weights to
#'   generate (typically `nlat`).
#' @param N Numeric mode at which the exponential taper weight equals
#'   `exp(-1) = 0.3678794`.
#' @param r Positive integer exponent controlling the steepness of the
#'   exponential filter.
#' @return
#' \itemize{
#'   \item `tri_trunc`, `rhomb_trunc`, `exp_tapersh`: Named list
#'     `list(a = a, b = b)` with truncated or tapered coefficients.
#'   \item `tri_trunC`, `rhomb_trunC`, `exp_tapershC`: Truncated or tapered
#'     array of the same dimension as `ab`.
#'   \item `exp_tapersh_wgts`: Numeric vector of length `M` containing the
#'     isotropic spectral weights for degrees `n = 0, 1, ..., M - 1`.
#' }
#' @examples
#' w <- exp_tapersh_wgts(45, 29, 4)
#' round(w[1:10], 5)
#'
#' nlon <- 36
#' nlat <- 19
#' g <- matrix(rnorm(nlon * nlat), nlon, nlat)
#' ab <- shaeC(g)
#' ab_t10 <- tri_trunC(ab, 10)
#' ab_r10 <- rhomb_trunC(ab, 10)
#' ab_tap <- exp_tapershC(ab, 10, 4)
#' @export
tri_trunc <- function(a, b, T) {
  if (!is.numeric(a) || !is.numeric(b) || !identical(dim(a), dim(b)) || length(dim(a)) < 2L) {
    stop("'a' and 'b' must be numeric arrays of identical dimensions (N, nlat, ...).")
  }
  T <- as.integer(T)[1L]
  if (is.na(T) || T < 0L) {
    stop("'T' must be a non-negative integer.")
  }
  d <- dim(a)
  nm <- d[1L]
  nlat <- d[2L]
  nt <- if (length(d) == 2L) 1L else as.integer(prod(d[-c(1L, 2L)]))
  storage.mode(a) <- "double"
  storage.mode(b) <- "double"
  dim(a) <- c(nm, nlat, nt)
  dim(b) <- c(nm, nlat, nt)

  mask <- outer(seq_len(nm), seq_len(nlat), function(m_idx, n_idx) {
    (n_idx < m_idx) | ((n_idx - 1L) > T)
  })
  for (k in seq_len(nt)) {
    ak <- a[, , k]
    bk <- b[, , k]
    ak[mask] <- 0.0
    bk[mask] <- 0.0
    a[, , k] <- ak
    b[, , k] <- bk
  }
  dim(a) <- d
  dim(b) <- d
  list(a = a, b = b)
}

#' @rdname tri_trunc
#' @export
tri_trunC <- function(ab, T) {
  sp <- .split_ab(ab, "ab")
  res <- tri_trunc(sp$a, sp$b, T)
  .combine_Trailing(list(res$a, res$b), dim(res$a))
}

#' @rdname tri_trunc
#' @export
rhomb_trunc <- function(a, b, T) {
  if (!is.numeric(a) || !is.numeric(b) || !identical(dim(a), dim(b)) || length(dim(a)) < 2L) {
    stop("'a' and 'b' must be numeric arrays of identical dimensions (N, nlat, ...).")
  }
  T <- as.integer(T)[1L]
  if (is.na(T) || T < 0L) {
    stop("'T' must be a non-negative integer.")
  }
  d <- dim(a)
  nm <- d[1L]
  nlat <- d[2L]
  nt <- if (length(d) == 2L) 1L else as.integer(prod(d[-c(1L, 2L)]))
  storage.mode(a) <- "double"
  storage.mode(b) <- "double"
  dim(a) <- c(nm, nlat, nt)
  dim(b) <- c(nm, nlat, nt)

  mask <- outer(seq_len(nm), seq_len(nlat), function(m_idx, n_idx) {
    m <- m_idx - 1L
    n <- n_idx - 1L
    (n < m) | (m > T) | ((n - m) > T)
  })
  for (k in seq_len(nt)) {
    ak <- a[, , k]
    bk <- b[, , k]
    ak[mask] <- 0.0
    bk[mask] <- 0.0
    a[, , k] <- ak
    b[, , k] <- bk
  }
  dim(a) <- d
  dim(b) <- d
  list(a = a, b = b)
}

#' @rdname tri_trunc
#' @export
rhomb_trunC <- function(ab, T) {
  sp <- .split_ab(ab, "ab")
  res <- rhomb_trunc(sp$a, sp$b, T)
  .combine_Trailing(list(res$a, res$b), dim(res$a))
}

#' @rdname tri_trunc
#' @export
exp_tapersh_wgts <- function(M, N, r) {
  M <- as.integer(M)[1L]
  N <- as.double(N)[1L]
  r <- as.integer(r)[1L]
  if (is.na(M) || M < 1L) stop("'M' must be a positive integer.")
  if (is.na(N) || N <= 0) stop("'N' must be a positive number.")
  if (is.na(r) || r < 1L) stop("'r' must be a positive integer.")
  n <- as.double(0:(M - 1L))
  denom <- N * (N + 1.0)
  exp(-((n * (n + 1.0)) / denom)^r)
}

#' @rdname tri_trunc
#' @export
exp_tapersh <- function(a, b, N, r) {
  if (!is.numeric(a) || !is.numeric(b) || !identical(dim(a), dim(b)) || length(dim(a)) < 2L) {
    stop("'a' and 'b' must be numeric arrays of identical dimensions (N_m, nlat, ...).")
  }
  d <- dim(a)
  nm <- d[1L]
  nlat <- d[2L]
  nt <- if (length(d) == 2L) 1L else as.integer(prod(d[-c(1L, 2L)]))
  w <- exp_tapersh_wgts(nlat, N, r)
  w_mat <- matrix(rep(w, each = nm), nrow = nm, ncol = nlat)
  storage.mode(a) <- "double"
  storage.mode(b) <- "double"
  dim(a) <- c(nm, nlat, nt)
  dim(b) <- c(nm, nlat, nt)
  for (k in seq_len(nt)) {
    a[, , k] <- a[, , k] * w_mat
    b[, , k] <- b[, , k] * w_mat
  }
  dim(a) <- d
  dim(b) <- d
  list(a = a, b = b)
}

#' @rdname tri_trunc
#' @export
exp_tapershC <- function(ab, N, r) {
  sp <- .split_ab(ab, "ab")
  res <- exp_tapersh(sp$a, sp$b, N, r)
  .combine_Trailing(list(res$a, res$b), dim(res$a))
}
