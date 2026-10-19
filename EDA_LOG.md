# EDA_LOG.md

## 2026-10-19: Development and Numerical Verification of the `spherepack` R Package

### 1. Objectives and Contents of the Experiment

- **Objective**: Create a CRAN-ready R package (`spherepack`) providing an R interface to the Fortran 77 SPHEREPACK 3.2 library (`spherepack3.2.tar.gz` by John C. Adams and Paul N. Swarztrauber), compatible with all spherical harmonic analysis, synthesis, differential, inverse, and truncation/tapering functions in the NCAR Command Language (NCL) (`https://www.ncl.ucar.edu/Document/Functions/spherical.shtml`), while leveraging the column-major array ordering shared between Fortran and R.
- **Design and Implementation**:
  - **64-bit Double Precision Conversion**: Converted the 86 core SPHEREPACK 3.2 Fortran 77 source files (`src/*.f`) from single precision (`REAL`) to 64-bit floating-point precision (`DOUBLE PRECISION`) with 32-bit integers (`INTEGER`).
  - **Fortran Standards Modernization**:
    - Replaced legacy Fortran 66 assumed-size `(1)` and `(...,1)` dummy array bounds with standard Fortran 77 `(*)` and `(...,*)`.
    - Replaced Fortran 66/77 arithmetic `IF` statements (`if (expr) L1, L2, L3`) with standard logical `IF`/`GO TO` statements.
    - Unshared nested `DO` loop termination labels so that each labeled `DO` statement has a unique termination label, eliminating all Fortran 2018 deleted-feature compiler warnings under `gfortran`.
    - Prefixed internal helper subroutines in `shpg.f` (`shpg_*`) and bundled LINPACK/BLAS routines in `shpe.f` (`sp_*`) to prevent symbol collisions.
    - Replaced Fortran `RANDOM_SEED` and `RANDOM_NUMBER` calls in `shpe.f` and `shpg.f` with a C helper `sp_runif` backed by R's `GetRNGstate()`, `unif_rand()`, and `PutRNGstate()`, ensuring reproducibility with `set.seed(514)` and compliance with CRAN compiled-code policies.
  - **Column-Major Memory Layout**:
    - Spatial arrays in R use geophysical column-major ordering `(nlon, nlat, ...)` with South-to-North latitudes (`ig = 0` in `geo2math.f`) and no cyclic longitude points, matching NCL's row-major `(..., nlat, nlon)` memory layout with zero transposition overhead.
    - Spherical harmonic coefficient arrays (`a`, `b`, `br`, `bi`, `cr`, `ci`) use SPHEREPACK's native column-major ordering `(N, nlat, ...)`.
  - **Native `.Call()` C Interface**:
    - Implemented `src/spherepack_init.c` using R's native `.Call()` API, `F77_CALL()`, `R_alloc()` workspace management, and `R_registerRoutines()` with `R_useDynamicSymbols(dll, FALSE)`.

### 2. Generated File Names and Main Functions

- **Package Metadata & Configuration**:
  - `DESCRIPTION`: Package metadata with author Takeshi Enomoto (`takeshi.enomoto@gmail.com`) and SPHEREPACK 3.2 authors John Adams and Paul Swarztrauber.
  - `NAMESPACE`: Exports all 88 user-facing functions and registers the dynamic library `spherepack`.
  - `.Rbuildignore`: Excludes non-package files (`spherepack3.2.tar.gz`, `AGENTS.md`, `EDA_LOG.md`, Fossil checkout files).
- **C and Fortran Source Files (`src/`)**:
  - `src/spherepack_init.c`: Native `.Call()` C wrappers and symbol registration (`c_gaqd`, `c_sha`, `c_shs`, `c_vha`, `c_vhs`, `c_uv2vrdv_op`, `c_vrdv2uv_op`, `c_uv2sfvp`, `c_sfvp2uv`, `c_grads`, `c_igrads`, `c_lderuv`, `c_laps`, `c_ilaps`, `c_lapv`, `c_ilapv`, `c_trunc`, `c_exp_tapersh`).
  - `src/*.f` (86 Fortran 77 files): Double-precision SPHEREPACK 3.2 routines (`alf.f`, `divec.f`, `dives.f`, `divgc.f`, `divgs.f`, `gaqd.f`, `geo2math.f`, `gradec.f`, `grades.f`, `gradgc.f`, `gradgs.f`, `hrfft.f`, `idivec.f`, `idives.f`, `idivgc.f`, `idivgs.f`, `idvtec.f`, `idvtes.f`, `idvtgc.f`, `idvtgs.f`, `igradec.f`, `igrades.f`, `igradgc.f`, `igradgs.f`, `ihgeod.f`, `isfvpec.f`, `isfvpes.f`, `isfvpgc.f`, `isfvpgs.f`, `islapec.f`, `islapes.f`, `islapgc.f`, `islapgs.f`, `ivlapec.f`, `ivlapes.f`, `ivlapgc.f`, `ivlapgs.f`, `ivrtec.f`, `ivrtes.f`, `ivrtgc.f`, `ivrtgs.f`, `sfvpec.f`, `sfvpes.f`, `sfvpgc.f`, `sfvpgs.f`, `shaec.f`, `shaes.f`, `shagc.f`, `shags.f`, `shigc.f`, `shigs.f`, `shpe.f`, `shpg.f`, `shsec.f`, `shses.f`, `shsgc.f`, `shsgs.f`, `slapec.f`, `slapes.f`, `slapgc.f`, `slapgs.f`, `sphcom.f`, `sshifte.f`, `trssph.f`, `trvsph.f`, `vhaec.f`, `vhaes.f`, `vhagc.f`, `vhags.f`, `vhsec.f`, `vhses.f`, `vhsgc.f`, `vhsgs.f`, `vlapec.f`, `vlapes.f`, `vlapgc.f`, `vlapgs.f`, `vrtec.f`, `vrtes.f`, `vrtgc.f`, `vrtgs.f`, `vshifte.f`, `vtsec.f`, `vtses.f`, `vtsgc.f`, `vtsgs.f`).
- **R Interface & Documentation (`R/`, `man/`, `tests/`)**:
  - `R/spherepack.R`: Implements all 88 exported functions:
    - Sphere radius & Gaussian quadrature: `get_sphere_radius`, `set_sphere_radius`, `gaqd`
    - Scalar analysis & synthesis: `shaec`, `shaeC`, `shagc`, `shagC`, `shsec`, `shseC`, `shsgc`, `shsgC`, `shsgc_R42`, `shsgc_R42_Wrap`
    - Vector analysis & synthesis: `vhaec`, `vhaeC`, `vhagc`, `vhagC`, `vhsec`, `vhseC`, `vhsgc`, `vhsgC`
    - Divergence & vorticity: `uv2dvf`, `uv2dvF`, `uv2dvF_Wrap`, `uv2dvg`, `uv2dvG`, `uv2dvG_Wrap`, `uv2vrf`, `uv2vrF`, `uv2vrF_Wrap`, `uv2vrg`, `uv2vrG`, `uv2vrG_Wrap`, `uv2vrdvf`, `uv2vrdvF`, `uv2vrdvg`, `uv2vrdvG`
    - Wind reconstruction from divergence/vorticity: `dv2uvf`, `dv2uvF`, `dv2uvF_Wrap`, `dv2uvg`, `dv2uvG`, `dv2uvG_Wrap`, `vr2uvf`, `vr2uvF`, `vr2uvF_Wrap`, `vr2uvg`, `vr2uvG`, `vr2uvG_Wrap`, `vrdv2uvf`, `vrdv2uvF`, `vrdv2uvg`, `vrdv2uvG`
    - Streamfunction & velocity potential: `uv2sfvpf`, `uv2sfvpF`, `uv2sfvpg`, `uv2sfvpG`, `sfvp2uvf`, `sfvp2uvg`
    - Gradients, inverse gradients & latitudinal derivatives: `gradsf`, `gradsg`, `igradsf`, `igradsF`, `igradsg`, `igradsG`, `lderuvf`, `lderuvg`
    - Scalar & vector Laplacians and inverse Laplacians: `lapsf`, `lapsF`, `lapsg`, `lapsG`, `ilapsf`, `ilapsF`, `ilapsF_Wrap`, `ilapsg`, `ilapsG`, `ilapsG_Wrap`, `lapvf`, `lapvg`, `ilapvf`, `ilapvg`
    - Spectral truncation & tapering: `tri_trunc`, `tri_trunC`, `rhomb_trunc`, `rhomb_trunC`, `exp_tapersh_wgts`, `exp_tapersh`, `exp_tapershC`
  - `man/*.Rd` (`gaqd.Rd`, `get_sphere_radius.Rd`, `gradsf.Rd`, `lapsf.Rd`, `shaec.Rd`, `tri_trunc.Rd`, `uv2dvf.Rd`, `uv2sfvpf.Rd`, `vhaec.Rd`): Complete R documentation with working examples.
  - `tests/test_spherepack.R`: Reproducible test suite (`seed <- 514; set.seed(seed)`) verifying all 88 exported functions.

### 3. Obtained Knowledge and Interpretation of Numerical Results

1. **Verification against Official SPHEREPACK 3.2 Test Suites**:
   - All 12 official SPHEREPACK 3.2 test drivers (`tgaqd`, `tsha`, `tvha`, `tdiv`, `tvrt`, `tgrad`, `tslap`, `tvlap`, `tidvt`, `tvts`, `tshpe`, `tshpg`) passed with `ierr = 0` and 64-bit machine-precision errors (`~1e-14` to `1e-16`), confirming that the `DOUBLE PRECISION` conversion preserved exact numerical constants and recurrence relations.
2. **Verification against Analytic Spherical Harmonics in R (`tests/test_spherepack.R`)**:
   - Using the test scalar potential $\psi(\lambda, \phi) = \sin\phi + \cos^2\phi \cos(2\lambda)$ (which is a linear combination of $Y_1^0$ with Laplacian eigenvalue $-2$ and $Y_2^2$ with Laplacian eigenvalue $-6$):
     - Forward/inverse scalar and vector harmonic transforms (`shaec`/`shsec`, `shagc`/`shsgc`, `vhaec`/`vhsec`, `vhagc`/`vhsgc`) achieved round-trip accuracy better than $10^{-13}$.
     - Differential operators (`gradsf`, `gradsg`, `lapsf`, `lapsg`, `uv2dvf`, `uv2dvg`, `uv2vrf`, `uv2vrg`, `uv2vrdvf`, `uv2vrdvg`, `lderuvf`, `lderuvg`) and inverse operators (`igradsf`, `igradsg`, `ilapsf`, `ilapsg`, `dv2uvf`, `vr2uvf`, `vrdv2uvf`, `uv2sfvpf`, `sfvp2uvf`, `ilapvf`, `ilapvg`) matched analytic solutions within $5 \times 10^{-13}$.
     - Triangular truncation `tri_trunc(..., T = 1)` removed the $n = 2$ mode ($\cos^2\phi \cos 2\lambda$) and preserved the $n = 1$ mode ($\sin\phi$) to within $10^{-13}$.
   - `R CMD check --no-manual spherepack_0.1.0.tar.gz` completed with `Status: OK` (0 errors, 0 warnings, 0 notes).

### 4. Execution Environment (`sessionInfo()`)

```
R version 4.5.3 (2026-03-11)
Platform: aarch64-apple-darwin25.3.0
Running under: macOS 27.0.1

Matrix products: default
BLAS:   /opt/local/Library/Frameworks/R.framework/Versions/4.5/Resources/lib/libRblas.dylib 
LAPACK: /opt/local/Library/Frameworks/R.framework/Versions/4.5/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_NZ.UTF-8/en_NZ.UTF-8/en_NZ.UTF-8/C/en_NZ.UTF-8/en_NZ.UTF-8

time zone: Asia/Tokyo
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.5.3
```

---

## 2026-10-19: Addition of Package Vignette (`vignettes/spherepack.Rmd`)

### 1. Objectives and Contents of the Experiment

- **Objective**: Provide a comprehensive package vignette introducing `spherepack`'s conventions, mathematical formulations, and geophysical workflows using base R graphics and reproducible random seeds (`seed <- 514; set.seed(seed)`).
- **Contents**:
  - Grid types (`f`/`F`/`ec` fixed grids vs. `g`/`G`/`gc` Gaussian grids via `gaqd()`) and zero-copy column-major `(nlon, nlat, ...)` array ordering shared between R/Fortran and row-major NCL.
  - Scalar harmonic analysis/synthesis (`shagc`, `shsgc`), triangular truncation (`tri_trunc`), and isotropic exponential spectral tapering (`exp_tapersh`, `exp_tapersh_wgts`).
  - Verification of scalar gradient (`gradsg`, `igradsg`) and scalar Laplacian (`lapsg`, `ilapsg`) against analytical spherical harmonic solutions.
  - Full Helmholtz wind decomposition on Earth's sphere ($R = 6.37122 \times 10^6\text{ m}$) into streamfunction/rotational wind (`uv2sfvpg`, `vr2uvg`) and velocity potential/divergent wind (`uv2vrdvg`, `dv2uvg`), visualized with base R `contour()` and `arrows()`.
  - 3D array broadcasting `(nlon, nlat, nlev)` and `_Wrap` metadata preservation (`uv2dvG_Wrap`).
  - Complete reference table of all 88 exported functions.

### 2. Generated File Names and Main Functions

- `DESCRIPTION`: Added `Suggests: knitr, rmarkdown` and `VignetteBuilder: knitr`.
- `vignettes/spherepack.Rmd`: Package vignette demonstrating `get_sphere_radius`, `set_sphere_radius`, `gaqd`, `shagc`, `shsgc`, `tri_trunc`, `exp_tapersh_wgts`, `exp_tapersh`, `gradsg`, `igradsg`, `lapsg`, `ilapsg`, `sfvp2uvg`, `uv2vrdvg`, `uv2sfvpg`, `vr2uvg`, `dv2uvg`, and `uv2dvG_Wrap`.

### 3. Obtained Knowledge and Interpretation of Results

- On a $64 \times 32$ Gaussian grid, scalar harmonic round-trip error and $T = 2$ triangular truncation error of a multi-scale wave field were at machine precision ($\sim 10^{-15}$), and spherical gradient/Laplacian operators and their inverses matched exact analytical expressions to $\sim 10^{-14}$.
- Helmholtz decomposition of a planetary-scale wind field on Earth ($R = 6.37122 \times 10^6\text{ m}$) recovered the streamfunction and velocity potential with relative error $< 10^{-14}$, and the sum of rotational (`vr2uvg`) and divergent (`dv2uvg`) wind components matched the total wind to $\sim 10^{-14}\text{ m s}^{-1}$.
- On macOS 27.0.1 (`Darwin Kernel Version 27.0.0`) with MacPorts `gfortran-mp-14` (`aarch64-apple-darwin25`), setting `MACOSX_DEPLOYMENT_TARGET=26.0` prevents GCC 14's Darwin kernel version fallback (`27 - 9 = 18.0`) when invoking the Clang assembler.
- `R CMD check --no-manual spherepack_0.1.0.tar.gz` passed with `Status: OK` (0 errors, 0 warnings, 0 notes), including vignette building and re-building checks.

