# Comprehensive Test Suite for spherepack (all 88 exported functions)
# Following AGENTS.md reproducibility guidelines
seed <- 514
set.seed(seed)

library(spherepack)

assert_near <- function(actual, expected, tol = 1e-10, label = "") {
  if (!is.numeric(actual) || !is.numeric(expected) || length(actual) != length(expected)) {
    stop(sprintf("FAIL [%s]: invalid type or length mismatch", label))
  }
  err <- max(abs(actual - expected))
  if (is.na(err) || err > tol) {
    stop(sprintf("FAIL [%s]: max abs error = %.3e > tol = %.3e", label, err, tol))
  }
}

# 1. Test sphere radius getter/setter (get_sphere_radius, set_sphere_radius)
r0 <- get_sphere_radius()
assert_near(r0, 6.37122e6, tol = 1e-2, label = "default sphere radius")
res_set <- set_sphere_radius(1.0)
stopifnot(is.null(res_set))
assert_near(get_sphere_radius(), 1.0, tol = 1e-15, label = "sphere radius set to 1.0")
set_sphere_radius(r0)
assert_near(get_sphere_radius(), r0, tol = 1e-2, label = "sphere radius restored")

# 2. Test Gaussian quadrature (gaqd)
nlat_g <- 18L
nlon_g <- 36L
g <- gaqd(nlat_g)
stopifnot(length(g$theta) == nlat_g, length(g$lat) == nlat_g, length(g$wts) == nlat_g)
assert_near(sum(g$wts), 2.0, tol = 1e-14, label = "sum of Gaussian weights == 2")
stopifnot(all(diff(g$lat) > 0)) # South-to-North

# Set up test grids (fixed and Gaussian) with unit radius for exact analytic checks
set_sphere_radius(1.0)

nlat_f <- 19L
nlon_f <- 36L
lat_f <- seq(-90, 90, length.out = nlat_f) * pi / 180
lon_f <- (0:(nlon_f - 1L)) * (2 * pi / nlon_f)

lat_g <- g$lat * pi / 180
lon_g <- (0:(nlon_g - 1L)) * (2 * pi / nlon_g)

# Analytic scalar field: Psi(lon, lat) = sin(lat) + cos^2(lat) * cos(2 * lon)
# Notice:
# - P_1^0(sin lat) = sin(lat), eigenvalue of Laplacian = -1*(1+1) = -2
# - P_2^2(sin lat) cos(2 lon) = 3 cos^2(lat) cos(2 lon), eigenvalue of Laplacian = -2*(2+1) = -6
# Therefore:
#   Laplacian(Psi) = -2 * sin(lat) - 6 * cos^2(lat) * cos(2 * lon)
#   dPsi/dlon / cos(lat) = -2 * cos(lat) * sin(2 * lon)
#   dPsi/dlat            = cos(lat) - 2 * sin(lat) * cos(lat) * cos(2 * lon)
make_scalar_grid <- function(lon, lat) {
  nlon <- length(lon)
  nlat <- length(lat)
  psi <- matrix(0.0, nrow = nlon, ncol = nlat)
  lap_psi <- matrix(0.0, nrow = nlon, ncol = nlat)
  grad_lon <- matrix(0.0, nrow = nlon, ncol = nlat)
  grad_lat <- matrix(0.0, nrow = nlon, ncol = nlat)
  for (j in seq_len(nlat)) {
    slat <- sin(lat[j])
    clat <- cos(lat[j])
    for (i in seq_len(nlon)) {
      c2lon <- cos(2 * lon[i])
      s2lon <- sin(2 * lon[i])
      psi[i, j] <- slat + (clat^2) * c2lon
      lap_psi[i, j] <- -2.0 * slat - 6.0 * (clat^2) * c2lon
      grad_lon[i, j] <- -2.0 * clat * s2lon
      grad_lat[i, j] <- clat - 2.0 * slat * clat * c2lon
    }
  }
  list(psi = psi, lap = lap_psi, gx = grad_lon, gy = grad_lat)
}

sf_f <- make_scalar_grid(lon_f, lat_f)
sf_g <- make_scalar_grid(lon_g, lat_g)

# 3. Scalar Harmonic Analysis & Synthesis
# (shaec, shsec, shagc, shsgc, shaeC, shseC, shagC, shsgC, shsgc_R42, shsgc_R42_Wrap)
ab_f <- shaec(sf_f$psi)
psi_rec_f <- shsec(ab_f$a, ab_f$b, nlon = nlon_f)
assert_near(psi_rec_f, sf_f$psi, tol = 1e-13, label = "shaec/shsec round-trip")

ab_fC <- shaeC(sf_f$psi)
stopifnot(identical(dim(ab_fC), c(min(nlat_f, (nlon_f + 2L) %/% 2L), nlat_f, 2L)))
psi_rec_fC <- shseC(ab_fC, nlon = nlon_f)
assert_near(psi_rec_fC, sf_f$psi, tol = 1e-13, label = "shaeC/shseC round-trip")

ab_g <- shagc(sf_g$psi)
psi_rec_g <- shsgc(ab_g$a, ab_g$b, nlon = nlon_g)
assert_near(psi_rec_g, sf_g$psi, tol = 1e-13, label = "shagc/shsgc round-trip")

ab_gC <- shagC(sf_g$psi)
psi_rec_gC <- shsgC(ab_gC, nlon = nlon_g)
assert_near(psi_rec_gC, sf_g$psi, tol = 1e-13, label = "shagC/shsgC round-trip")

# Test R42 synthesis (128 x 108 Gaussian grid from 43 x 43 spectral coefficients)
a_r42 <- matrix(0.0, nrow = 43L, ncol = 43L)
b_r42 <- matrix(0.0, nrow = 43L, ncol = 43L)
a_r42[1L, 2L] <- 1.0
attr(a_r42, "units") <- "K"
g_r42 <- shsgc_R42(a_r42, b_r42)
g_r42_w <- shsgc_R42_Wrap(a_r42, b_r42)
stopifnot(identical(dim(g_r42), c(128L, 108L)))
assert_near(as.numeric(g_r42_w), as.numeric(g_r42), tol = 1e-15, label = "shsgc_R42_Wrap values")
stopifnot(identical(attr(g_r42_w, "units"), "K"))

# 4. Gradients and Inverse Gradients (gradsf, gradsg, igradsf, igradsF, igradsg, igradsG)
gr_f <- gradsf(sf_f$psi)
assert_near(gr_f$gzx, sf_f$gx, tol = 1e-12, label = "gradsf gzx")
assert_near(gr_f$gzy, sf_f$gy, tol = 1e-12, label = "gradsf gzy")

gr_g <- gradsg(sf_g$psi)
assert_near(gr_g$gzx, sf_g$gx, tol = 1e-12, label = "gradsg gzx")
assert_near(gr_g$gzy, sf_g$gy, tol = 1e-12, label = "gradsg gzy")

assert_near(igradsf(gr_f$gzx, gr_f$gzy), sf_f$psi, tol = 1e-12, label = "igradsf")
assert_near(igradsF(gr_f$gzx, gr_f$gzy), sf_f$psi, tol = 1e-12, label = "igradsF")
assert_near(igradsg(gr_g$gzx, gr_g$gzy), sf_g$psi, tol = 1e-12, label = "igradsg")
assert_near(igradsG(gr_g$gzx, gr_g$gzy), sf_g$psi, tol = 1e-12, label = "igradsG")

# 5. Scalar Laplacian & Inverse Laplacian
# (lapsf, lapsg, lapsF, lapsG, ilapsf, ilapsg, ilapsF, ilapsG, ilapsF_Wrap, ilapsG_Wrap)
lap_f <- lapsf(sf_f$psi)
assert_near(lap_f, sf_f$lap, tol = 1e-12, label = "lapsf analytic")
assert_near(lapsF(sf_f$psi), sf_f$lap, tol = 1e-12, label = "lapsF analytic")

lap_g <- lapsg(sf_g$psi)
assert_near(lap_g, sf_g$lap, tol = 1e-12, label = "lapsg analytic")
assert_near(lapsG(sf_g$psi), sf_g$lap, tol = 1e-12, label = "lapsG analytic")

assert_near(ilapsf(lap_f), sf_f$psi, tol = 1e-12, label = "ilapsf")
assert_near(ilapsF(lap_f), sf_f$psi, tol = 1e-12, label = "ilapsF")
assert_near(ilapsF_Wrap(lap_f), sf_f$psi, tol = 1e-12, label = "ilapsF_Wrap")
assert_near(ilapsg(lap_g), sf_g$psi, tol = 1e-12, label = "ilapsg")
assert_near(ilapsG(lap_g), sf_g$psi, tol = 1e-12, label = "ilapsG")
assert_near(ilapsG_Wrap(lap_g), sf_g$psi, tol = 1e-12, label = "ilapsG_Wrap")

# 6. Vector Harmonics, Divergence, Vorticity, Streamfunction & Velocity Potential
# Construct a wind field from known streamfunction (sf) and velocity potential (vp):
# Let sf = sf_f$psi, vp = 0.5 * sf_f$psi
# Then vorticity = Laplacian(sf), divergence = Laplacian(vp) = 0.5 * Laplacian(sf)
sf_true_f <- sf_f$psi
vp_true_f <- 0.5 * sf_f$psi
vrt_true_f <- sf_f$lap
div_true_f <- 0.5 * sf_f$lap

uv_f <- sfvp2uvf(sf_true_f, vp_true_f)

sf_true_g <- sf_g$psi
vp_true_g <- 0.5 * sf_g$psi
vrt_true_g <- sf_g$lap
div_true_g <- 0.5 * sf_g$lap

uv_g <- sfvp2uvg(sf_true_g, vp_true_g)

# Vector harmonic analysis & synthesis (vhaec, vhsec, vhagc, vhsgc, vhaeC, vhseC, vhagC, vhsgC)
vh_f <- vhaec(uv_f$u, uv_f$v)
uv_rec_f <- vhsec(vh_f$br, vh_f$bi, vh_f$cr, vh_f$ci, nlon = nlon_f)
assert_near(uv_rec_f$u, uv_f$u, tol = 1e-12, label = "vhaec/vhsec u")
assert_near(uv_rec_f$v, uv_f$v, tol = 1e-12, label = "vhaec/vhsec v")

vh_fC <- vhaeC(uv_f$u, uv_f$v)
stopifnot(identical(dim(vh_fC), c(min(nlat_f, (nlon_f + 1L) %/% 2L), nlat_f, 4L)))
uv_rec_fC <- vhseC(vh_fC, nlon = nlon_f)
assert_near(uv_rec_fC[, , 1L], uv_f$u, tol = 1e-12, label = "vhaeC/vhseC u")
assert_near(uv_rec_fC[, , 2L], uv_f$v, tol = 1e-12, label = "vhaeC/vhseC v")

vh_g <- vhagc(uv_g$u, uv_g$v)
uv_rec_g <- vhsgc(vh_g$br, vh_g$bi, vh_g$cr, vh_g$ci, nlon = nlon_g)
assert_near(uv_rec_g$u, uv_g$u, tol = 1e-12, label = "vhagc/vhsgc u")
assert_near(uv_rec_g$v, uv_g$v, tol = 1e-12, label = "vhagc/vhsgc v")

vh_gC <- vhagC(uv_g$u, uv_g$v)
uv_rec_gC <- vhsgC(vh_gC, nlon = nlon_g)
assert_near(uv_rec_gC[, , 1L], uv_g$u, tol = 1e-12, label = "vhagC/vhsgC u")
assert_near(uv_rec_gC[, , 2L], uv_g$v, tol = 1e-12, label = "vhagC/vhsgC v")

# Divergence & Vorticity (uv2dvf/g/F/G/_Wrap, uv2vrf/g/F/G/_Wrap, uv2vrdvf/g/F/G, uv2sfvpf/g/F/G)
assert_near(uv2dvf(uv_f$u, uv_f$v), div_true_f, tol = 1e-12, label = "uv2dvf")
assert_near(uv2dvF(uv_f$u, uv_f$v), div_true_f, tol = 1e-12, label = "uv2dvF")
assert_near(uv2dvF_Wrap(uv_f$u, uv_f$v), div_true_f, tol = 1e-12, label = "uv2dvF_Wrap")
assert_near(uv2dvg(uv_g$u, uv_g$v), div_true_g, tol = 1e-12, label = "uv2dvg")
assert_near(uv2dvG(uv_g$u, uv_g$v), div_true_g, tol = 1e-12, label = "uv2dvG")
assert_near(uv2dvG_Wrap(uv_g$u, uv_g$v), div_true_g, tol = 1e-12, label = "uv2dvG_Wrap")

assert_near(uv2vrf(uv_f$u, uv_f$v), vrt_true_f, tol = 1e-12, label = "uv2vrf")
assert_near(uv2vrF(uv_f$u, uv_f$v), vrt_true_f, tol = 1e-12, label = "uv2vrF")
assert_near(uv2vrF_Wrap(uv_f$u, uv_f$v), vrt_true_f, tol = 1e-12, label = "uv2vrF_Wrap")
assert_near(uv2vrg(uv_g$u, uv_g$v), vrt_true_g, tol = 1e-12, label = "uv2vrg")
assert_near(uv2vrG(uv_g$u, uv_g$v), vrt_true_g, tol = 1e-12, label = "uv2vrG")
assert_near(uv2vrG_Wrap(uv_g$u, uv_g$v), vrt_true_g, tol = 1e-12, label = "uv2vrG_Wrap")

vrdv_f <- uv2vrdvf(uv_f$u, uv_f$v)
assert_near(vrdv_f$vr, vrt_true_f, tol = 1e-12, label = "uv2vrdvf vr")
assert_near(vrdv_f$dv, div_true_f, tol = 1e-12, label = "uv2vrdvf dv")
vrdv_fC <- uv2vrdvF(uv_f$u, uv_f$v)
assert_near(vrdv_fC[, , 1L], vrt_true_f, tol = 1e-12, label = "uv2vrdvF vr")
assert_near(vrdv_fC[, , 2L], div_true_f, tol = 1e-12, label = "uv2vrdvF dv")

vrdv_g <- uv2vrdvg(uv_g$u, uv_g$v)
assert_near(vrdv_g$vr, vrt_true_g, tol = 1e-12, label = "uv2vrdvg vr")
assert_near(vrdv_g$dv, div_true_g, tol = 1e-12, label = "uv2vrdvg dv")
vrdv_gC <- uv2vrdvG(uv_g$u, uv_g$v)
assert_near(vrdv_gC[, , 1L], vrt_true_g, tol = 1e-12, label = "uv2vrdvG vr")
assert_near(vrdv_gC[, , 2L], div_true_g, tol = 1e-12, label = "uv2vrdvG dv")

sfvp_f <- uv2sfvpf(uv_f$u, uv_f$v)
assert_near(sfvp_f$sf, sf_true_f, tol = 1e-12, label = "uv2sfvpf sf")
assert_near(sfvp_f$vp, vp_true_f, tol = 1e-12, label = "uv2sfvpf vp")
sfvp_fC <- uv2sfvpF(uv_f$u, uv_f$v)
assert_near(sfvp_fC[, , 1L], sf_true_f, tol = 1e-12, label = "uv2sfvpF sf")
assert_near(sfvp_fC[, , 2L], vp_true_f, tol = 1e-12, label = "uv2sfvpF vp")

sfvp_g <- uv2sfvpg(uv_g$u, uv_g$v)
assert_near(sfvp_g$sf, sf_true_g, tol = 1e-12, label = "uv2sfvpg sf")
assert_near(sfvp_g$vp, vp_true_g, tol = 1e-12, label = "uv2sfvpg vp")
sfvp_gC <- uv2sfvpG(uv_g$u, uv_g$v)
assert_near(sfvp_gC[, , 1L], sf_true_g, tol = 1e-12, label = "uv2sfvpG sf")
assert_near(sfvp_gC[, , 2L], vp_true_g, tol = 1e-12, label = "uv2sfvpG vp")

# Inverse Vorticity/Divergence to Winds (dv2uvf/g/F/G/_Wrap, vr2uvf/g/F/G/_Wrap, vrdv2uvf/g/F/G)
uv_from_vrdv_f <- vrdv2uvf(vrt_true_f, div_true_f)
assert_near(uv_from_vrdv_f$u, uv_f$u, tol = 1e-12, label = "vrdv2uvf u")
assert_near(uv_from_vrdv_f$v, uv_f$v, tol = 1e-12, label = "vrdv2uvf v")
uv_from_vrdv_fC <- vrdv2uvF(vrt_true_f, div_true_f)
assert_near(uv_from_vrdv_fC[, , 1L], uv_f$u, tol = 1e-12, label = "vrdv2uvF u")
assert_near(uv_from_vrdv_fC[, , 2L], uv_f$v, tol = 1e-12, label = "vrdv2uvF v")

uv_from_vrdv_g <- vrdv2uvg(vrt_true_g, div_true_g)
assert_near(uv_from_vrdv_g$u, uv_g$u, tol = 1e-12, label = "vrdv2uvg u")
assert_near(uv_from_vrdv_g$v, uv_g$v, tol = 1e-12, label = "vrdv2uvg v")
uv_from_vrdv_gC <- vrdv2uvG(vrt_true_g, div_true_g)
assert_near(uv_from_vrdv_gC[, , 1L], uv_g$u, tol = 1e-12, label = "vrdv2uvG u")
assert_near(uv_from_vrdv_gC[, , 2L], uv_g$v, tol = 1e-12, label = "vrdv2uvG v")

# Irrotational + non-divergent decomposition sum equals total wind
uv_div_f <- dv2uvf(div_true_f)
uv_rot_f <- vr2uvf(vrt_true_f)
assert_near(uv_div_f$u + uv_rot_f$u, uv_f$u, tol = 1e-12, label = "dv2uvf + vr2uvf == u")
assert_near(uv_div_f$v + uv_rot_f$v, uv_f$v, tol = 1e-12, label = "dv2uvf + vr2uvf == v")
assert_near(dv2uvF(div_true_f)[, , 1L] + vr2uvF(vrt_true_f)[, , 1L], uv_f$u, tol = 1e-12, label = "dv2uvF + vr2uvF == u")
assert_near(dv2uvF_Wrap(div_true_f)[, , 1L] + vr2uvF_Wrap(vrt_true_f)[, , 1L], uv_f$u, tol = 1e-12, label = "dv2uvF_Wrap + vr2uvF_Wrap == u")

uv_div_g <- dv2uvg(div_true_g)
uv_rot_g <- vr2uvg(vrt_true_g)
assert_near(uv_div_g$u + uv_rot_g$u, uv_g$u, tol = 1e-12, label = "dv2uvg + vr2uvg == u")
assert_near(uv_div_g$v + uv_rot_g$v, uv_g$v, tol = 1e-12, label = "dv2uvg + vr2uvg == v")
assert_near(dv2uvG(div_true_g)[, , 1L] + vr2uvG(vrt_true_g)[, , 1L], uv_g$u, tol = 1e-12, label = "dv2uvG + vr2uvG == u")
assert_near(dv2uvG_Wrap(div_true_g)[, , 1L] + vr2uvG_Wrap(vrt_true_g)[, , 1L], uv_g$u, tol = 1e-12, label = "dv2uvG_Wrap + vr2uvG_Wrap == v")

# 7. Vector Laplacian and Inverse Vector Laplacian (lapvf, lapvg, ilapvf, ilapvg)
lv_f <- lapvf(uv_f$u, uv_f$v)
uv_ilv_f <- ilapvf(lv_f$ulap, lv_f$vlap)
assert_near(uv_ilv_f$u, uv_f$u, tol = 1e-12, label = "lapvf/ilapvf u")
assert_near(uv_ilv_f$v, uv_f$v, tol = 1e-12, label = "lapvf/ilapvf v")

lv_g <- lapvg(uv_g$u, uv_g$v)
uv_ilv_g <- ilapvg(lv_g$ulap, lv_g$vlap)
assert_near(uv_ilv_g$u, uv_g$u, tol = 1e-12, label = "lapvg/ilapvg u")
assert_near(uv_ilv_g$v, uv_g$v, tol = 1e-12, label = "lapvg/ilapvg v")

# 8. Latitudinal Derivatives of Winds (lderuvf, lderuvg)
# Test with pure zonal wind u(lat) = sin(lat)*cos(lat), v(lat) = 0 -> du/dlat = cos^2(lat) - sin^2(lat) = cos(2*lat)
u_zonal_f <- matrix(rep(sin(lat_f) * cos(lat_f), each = nlon_f), nrow = nlon_f, ncol = nlat_f)
v_zonal_f <- matrix(0.0, nrow = nlon_f, ncol = nlat_f)
du_dlat_f <- matrix(rep(cos(2 * lat_f), each = nlon_f), nrow = nlon_f, ncol = nlat_f)

ld_f <- lderuvf(u_zonal_f, v_zonal_f)
assert_near(ld_f$uy, du_dlat_f, tol = 1e-12, label = "lderuvf uy")
assert_near(ld_f$vy, v_zonal_f, tol = 1e-12, label = "lderuvf vy")

u_zonal_g <- matrix(rep(sin(lat_g) * cos(lat_g), each = nlon_g), nrow = nlon_g, ncol = nlat_g)
v_zonal_g <- matrix(0.0, nrow = nlon_g, ncol = nlat_g)
du_dlat_g <- matrix(rep(cos(2 * lat_g), each = nlon_g), nrow = nlon_g, ncol = nlat_g)

ld_g <- lderuvg(u_zonal_g, v_zonal_g)
assert_near(ld_g$uy, du_dlat_g, tol = 1e-12, label = "lderuvg uy")
assert_near(ld_g$vy, v_zonal_g, tol = 1e-12, label = "lderuvg vy")

# 9. Spectral Truncation and Tapering (tri_trunc, rhomb_trunc, tri_trunC, rhomb_trunC, exp_tapersh_wgts, exp_tapersh, exp_tapershC)
# Since sf_f$psi only has total wavenumbers n=1 and n=2:
# - T2 truncation keeps n <= 2 intact
# - T1 truncation zeroes out n = 2, leaving only sin(lat)
tr_t2 <- tri_trunc(ab_f$a, ab_f$b, T = 2L)
assert_near(shsec(tr_t2$a, tr_t2$b, nlon = nlon_f), sf_f$psi, tol = 1e-13, label = "tri_trunc T=2 preserves n<=2")

tr_t1 <- tri_trunc(ab_f$a, ab_f$b, T = 1L)
sin_lat_f <- matrix(rep(sin(lat_f), each = nlon_f), nrow = nlon_f, ncol = nlat_f)
assert_near(shsec(tr_t1$a, tr_t1$b, nlon = nlon_f), sin_lat_f, tol = 1e-13, label = "tri_trunc T=1 leaves n=1")

tr_t1C <- tri_trunC(ab_fC, T = 1L)
assert_near(shseC(tr_t1C, nlon = nlon_f), sin_lat_f, tol = 1e-13, label = "tri_trunC T=1 leaves n=1")

rh_t2 <- rhomb_trunc(ab_f$a, ab_f$b, T = 2L)
assert_near(shsec(rh_t2$a, rh_t2$b, nlon = nlon_f), sf_f$psi, tol = 1e-13, label = "rhomb_trunc T=2 preserves n<=2")
rh_t2C <- rhomb_trunC(ab_fC, T = 2L)
assert_near(shseC(rh_t2C, nlon = nlon_f), sf_f$psi, tol = 1e-13, label = "rhomb_trunC T=2 preserves n<=2")

wgts <- exp_tapersh_wgts(M = nlat_f, N = 5, r = 2)
stopifnot(length(wgts) == nlat_f)
assert_near(wgts[1L], 1.0, tol = 1e-15, label = "exp_tapersh_wgts n=0 is 1")
assert_near(wgts[6L], exp(-1.0), tol = 1e-15, label = "exp_tapersh_wgts n=N is exp(-1)")

tap <- exp_tapersh(ab_f$a, ab_f$b, N = 5, r = 2)
tapC <- exp_tapershC(ab_fC, N = 5, r = 2)
assert_near(tapC[, , 1L], tap$a, tol = 1e-15, label = "exp_tapershC a")
assert_near(tapC[, , 2L], tap$b, tol = 1e-15, label = "exp_tapershC b")

# 10. Multi-dimensional array support (nlon, nlat, ntime)
psi_3d <- array(c(sf_f$psi, 2.0 * sf_f$psi), dim = c(nlon_f, nlat_f, 2L))
lap_3d <- lapsf(psi_3d)
stopifnot(identical(dim(lap_3d), c(nlon_f, nlat_f, 2L)))
assert_near(lap_3d[, , 1L], sf_f$lap, tol = 1e-12, label = "3D lapsf slice 1")
assert_near(lap_3d[, , 2L], 2.0 * sf_f$lap, tol = 2e-12, label = "3D lapsf slice 2")

set_sphere_radius(r0)
cat("ALL 88 SPHEREPACK FUNCTION TESTS PASSED SUCCESSFULLY!\n")
