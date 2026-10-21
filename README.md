# spherepack

**Spherical Harmonic Analysis, Synthesis, and Differential Operators for R**

`spherepack` provides a 64-bit double-precision R interface to the **SPHEREPACK 3.2** Fortran 77 library for spherical harmonic analysis and synthesis of scalar and vector fields, as well as divergence, vorticity, streamfunction, velocity potential, gradient, scalar and vector Laplacians, inverse operators, latitudinal derivatives, spectral truncation, and exponential tapering on global equally spaced (fixed) and Gaussian grids.

All 88 exported functions follow the naming and calling conventions of the [NCAR Command Language (NCL) Spherical Harmonic Routines](https://www.ncl.ucar.edu/Document/Functions/spherical.shtml) while leveraging the column-major array ordering `(nlon, nlat, ...)` shared between Fortran and R for zero-transposition performance.

## Features

- **64-bit Double Precision**: Core SPHEREPACK 3.2 Fortran 77 routines converted to 64-bit `DOUBLE PRECISION`.
- **Two Global Grid Types**:
  - **Equally spaced (fixed) grids** (`f` / `F` / `ec` / `eC`): Uniform latitude spacing from $-90^\circ$ (South Pole) to $+90^\circ$ (North Pole).
  - **Gaussian grids** (`g` / `G` / `gc` / `gC`): Gaussian latitudes ordered South-to-North (excluding poles), computed via `gaqd()`.
- **NCL-Compatible Interfaces**:
  - Lowercase procedural variants (e.g., `shaec()`, `uv2vrdvg()`, `uv2sfvpg()`) returning named lists for multi-field outputs.
  - Uppercase functional variants (e.g., `shaeC()`, `uv2vrdvG()`, `uv2dvG()`) returning stacked numeric arrays.
  - `_Wrap` variants (e.g., `uv2dvG_Wrap()`, `ilapsG_Wrap()`) preserving array metadata (`dimnames`, `units`, `long_name`).
- **Multi-Dimensional Broadcasting**: Automatically iterates over arbitrary trailing dimensions such as vertical levels, time steps, or ensemble members `(nlon, nlat, ...)`.

## Installation

Install from a local source checkout or tarball:

```sh
R CMD INSTALL .
```

Or from within R:

```r
install.packages(".", repos = NULL, type = "source")
```

Alternatively using `pak`:

```r
install.packages("pak") # skip if already installed
pak::pak("tenomoto/spherepack")
```

## Quickstart

### 1. Gaussian Grid Generation and Scalar Harmonic Transform

```r
library(spherepack)

# Generate a 64 x 32 Gaussian grid (South-to-North latitudes)
nlon <- 64L
nlat <- 32L
g <- gaqd(nlat)

lon <- (0:(nlon - 1L)) * (2 * pi / nlon)
lat <- g$lat * pi / 180

# Analytic scalar field: sin(phi) + cos^2(phi) * cos(2 * lambda)
psi <- outer(lon, lat, function(lam, phi) {
  sin(phi) + (cos(phi)^2) * cos(2 * lam)
})

# Spherical harmonic analysis and synthesis on the Gaussian grid
ab <- shagc(psi)
psi_rec <- shsgc(ab$a, ab$b, nlon = nlon)
max(abs(psi_rec - psi))
#> [1] 4.440892e-15

# Triangular truncation at total wavenumber T = 1 (retains sin(phi) mode)
ab_t1 <- tri_trunc(ab$a, ab$b, T = 1L)
psi_t1 <- shsgc(ab_t1$a, ab_t1$b, nlon = nlon)
```

### 2. Scalar Gradient, Laplacian, and Inverse Operators

```r
# Use unit sphere (R = 1) for comparison with analytic eigenvalues
set_sphere_radius(1.0)

# Gradient and inverse gradient
grad <- gradsg(psi)
psi_from_grad <- igradsg(grad$gzx, grad$gzy)
max(abs(psi_from_grad - psi))
#> [1] 4.218847e-15

# Scalar Laplacian and inverse Laplacian
lap <- lapsg(psi)
psi_from_lap <- ilapsg(lap)
max(abs(psi_from_lap - psi))
#> [1] 4.662937e-15
```

### 3. Helmholtz Wind Decomposition (Vorticity, Divergence, Streamfunction, Velocity Potential)

```r
# Restore mean Earth radius (6.37122e6 m)
set_sphere_radius(6.37122e6)

# Define streamfunction (sf) and velocity potential (vp) in m^2/s
sf <- 1.0e7 * psi
vp <- 0.4e7 * outer(lon, lat, function(lam, phi) sin(2 * phi) * sin(lam))

# Synthesize horizontal wind (u, v) from (sf, vp)
uv <- sfvp2uvg(sf, vp)

# Compute relative vorticity (vr) and divergence (dv) from (u, v)
vrdv <- uv2vrdvg(uv$u, uv$v)

# Recover streamfunction and velocity potential from (u, v)
sfvp <- uv2sfvpg(uv$u, uv$v)

# Decompose total wind into rotational and divergent components
uv_rot <- vr2uvg(vrdv$vr)
uv_div <- dv2uvg(vrdv$dv)
max(abs((uv_rot$u + uv_div$u) - uv$u))
#> [1] 1.554312e-14
```

## Summary of Exported Functions

| Category | Fixed Grid (`f` / `F` / `ec`) | Gaussian Grid (`g` / `G` / `gc`) |
| :--- | :--- | :--- |
| **Sphere Radius & Quadrature** | `get_sphere_radius`, `set_sphere_radius` | `gaqd` |
| **Scalar Analysis** | `shaec`, `shaeC` | `shagc`, `shagC` |
| **Scalar Synthesis** | `shsec`, `shseC` | `shsgc`, `shsgC`, `shsgc_R42`, `shsgc_R42_Wrap` |
| **Vector Analysis** | `vhaec`, `vhaeC` | `vhagc`, `vhagC` |
| **Vector Synthesis** | `vhsec`, `vhseC` | `vhsgc`, `vhsgC` |
| **Divergence from Wind** | `uv2dvf`, `uv2dvF`, `uv2dvF_Wrap` | `uv2dvg`, `uv2dvG`, `uv2dvG_Wrap` |
| **Vorticity from Wind** | `uv2vrf`, `uv2vrF`, `uv2vrF_Wrap` | `uv2vrg`, `uv2vrG`, `uv2vrG_Wrap` |
| **Vorticity & Divergence** | `uv2vrdvf`, `uv2vrdvF` | `uv2vrdvg`, `uv2vrdvG` |
| **Divergent Wind from Divergence** | `dv2uvf`, `dv2uvF`, `dv2uvF_Wrap` | `dv2uvg`, `dv2uvG`, `dv2uvG_Wrap` |
| **Rotational Wind from Vorticity** | `vr2uvf`, `vr2uvF`, `vr2uvF_Wrap` | `vr2uvg`, `vr2uvG`, `vr2uvG_Wrap` |
| **Wind from Vorticity & Divergence** | `vrdv2uvf`, `vrdv2uvF` | `vrdv2uvg`, `vrdv2uvG` |
| **Streamfunction & Velocity Potential** | `uv2sfvpf`, `uv2sfvpF`, `sfvp2uvf` | `uv2sfvpg`, `uv2sfvpG`, `sfvp2uvg` |
| **Gradient & Inverse Gradient** | `gradsf`, `igradsf`, `igradsF` | `gradsg`, `igradsg`, `igradsG` |
| **Latitudinal Derivative of Wind** | `lderuvf` | `lderuvg` |
| **Scalar Laplacian & Inverse** | `lapsf`, `lapsF`, `ilapsf`, `ilapsF`, `ilapsF_Wrap` | `lapsg`, `lapsG`, `ilapsg`, `ilapsG`, `ilapsG_Wrap` |
| **Vector Laplacian & Inverse** | `lapvf`, `ilapvf` | `lapvg`, `ilapvg` |
| **Spectral Truncation & Tapering** | `tri_trunc`, `tri_trunC`, `rhomb_trunc`, `rhomb_trunC`, `exp_tapersh_wgts`, `exp_tapersh`, `exp_tapershC` | *(applies to both grid types)* |

## References and Links

- **Original SPHEREPACK Repository (NCAR)**:  
  <https://github.com/NCAR/NCAR-Classic-Libraries-for-Geophysics/tree/main/SpherePack>
- **NCAR Command Language (NCL) Spherical Harmonic Functions**:  
  <https://www.ncl.ucar.edu/Document/Functions/spherical.shtml>
- **Reference Paper**:  
  Adams, J. C., and P. N. Swarztrauber (1999): SPHEREPACK 3.0: A model development facility. *Monthly Weather Review*, **127**(8), 1872–1878, [https://doi.org/10.1175/1520-0493(1999)127<1872:SAMDF>2.0.CO;2](https://doi.org/10.1175/1520-0493%281999%29127%3C1872:SAMDF%3E2.0.CO;2)

## License

2-Clause BSD License. See [LICENSE.md](LICENSE.md) and [inst/COPYRIGHTS](inst/COPYRIGHTS) for full copyright and license details.
