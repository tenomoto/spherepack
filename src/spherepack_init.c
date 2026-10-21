#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Random.h>
#include <math.h>
#include <string.h>

/* Provide R's uniform RNG to Fortran 77 shpe.f and shpg.f */
void F77_SUB(sp_runif)(int *n, double *x) {
  int i, len = *n;
  GetRNGstate();
  for (i = 0; i < len; i++) {
    x[i] = unif_rand();
  }
  PutRNGstate();
}

/* Fortran 77 declarations from SPHEREPACK 3.2 (converted to DOUBLE PRECISION) */
extern void F77_NAME(gaqd)(int *nlat, double *theta, double *wts,
                           double *w, int *lwork, int *ierror);

extern void F77_NAME(geo2maths)(int *ig, int *nlon, int *nlat,
                                double *sg, double *sm, double *work);
extern void F77_NAME(math2geos)(int *ig, int *nlat, int *nlon,
                                double *sm, double *sg, double *work);
extern void F77_NAME(geo2mathv)(int *ig, int *nlon, int *nlat,
                                double *ug, double *vg,
                                double *vm, double *wm, double *work);
extern void F77_NAME(math2geov)(int *ig, int *nlat, int *nlon,
                                double *vm, double *wm,
                                double *ug, double *vg, double *work);

extern void F77_NAME(shaeci)(int *nlat, int *nlon, double *wshaec, int *lshaec,
                             double *dwork, int *ldwork, int *ierror);
extern void F77_NAME(shaec)(int *nlat, int *nlon, int *isym, int *nt,
                            double *g, int *idg, int *jdg,
                            double *a, double *b, int *mdab, int *ndab,
                            double *wshaec, int *lshaec,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(shagci)(int *nlat, int *nlon, double *wshagc, int *lshagc,
                             double *dwork, int *ldwork, int *ierror);
extern void F77_NAME(shagc)(int *nlat, int *nlon, int *isym, int *nt,
                            double *g, int *idg, int *jdg,
                            double *a, double *b, int *mdab, int *ndab,
                            double *wshagc, int *lshagc,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(shseci)(int *nlat, int *nlon, double *wshsec, int *lshsec,
                             double *dwork, int *ldwork, int *ierror);
extern void F77_NAME(shsec)(int *nlat, int *nlon, int *isym, int *nt,
                            double *g, int *idg, int *jdg,
                            double *a, double *b, int *mdab, int *ndab,
                            double *wshsec, int *lshsec,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(shsgci)(int *nlat, int *nlon, double *wshsgc, int *lshsgc,
                             double *dwork, int *ldwork, int *ierror);
extern void F77_NAME(shsgc)(int *nlat, int *nlon, int *isym, int *nt,
                            double *g, int *idg, int *jdg,
                            double *a, double *b, int *mdab, int *ndab,
                            double *wshsgc, int *lshsgc,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(vhaeci)(int *nlat, int *nlon, double *wvhaec, int *lvhaec,
                             double *dwork, int *ldwork, int *ierror);
extern void F77_NAME(vhaec)(int *nlat, int *nlon, int *ityp, int *nt,
                            double *v, double *w, int *idvw, int *jdvw,
                            double *br, double *bi, double *cr, double *ci,
                            int *mdab, int *ndab,
                            double *wvhaec, int *lvhaec,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(vhagci)(int *nlat, int *nlon, double *wvhagc, int *lvhagc,
                             double *dwork, int *ldwork, int *ierror);
extern void F77_NAME(vhagc)(int *nlat, int *nlon, int *ityp, int *nt,
                            double *v, double *w, int *idvw, int *jdvw,
                            double *br, double *bi, double *cr, double *ci,
                            int *mdab, int *ndab,
                            double *wvhagc, int *lvhagc,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(vhseci)(int *nlat, int *nlon, double *wvhsec, int *lvhsec,
                             double *dwork, int *ldwork, int *ierror);
extern void F77_NAME(vhsec)(int *nlat, int *nlon, int *ityp, int *nt,
                            double *v, double *w, int *idvw, int *jdvw,
                            double *br, double *bi, double *cr, double *ci,
                            int *mdab, int *ndab,
                            double *wvhsec, int *lvhsec,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(vhsgci)(int *nlat, int *nlon, double *wvhsgc, int *lvhsgc,
                             double *dwork, int *ldwork, int *ierror);
extern void F77_NAME(vhsgc)(int *nlat, int *nlon, int *ityp, int *nt,
                            double *v, double *w, int *idvw, int *jdvw,
                            double *br, double *bi, double *cr, double *ci,
                            int *mdab, int *ndab,
                            double *wvhsgc, int *lvhsgc,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(divec)(int *nlat, int *nlon, int *isym, int *nt,
                            double *dv, int *idv, int *jdv,
                            double *br, double *bi, int *mdb, int *ndb,
                            double *wshsec, int *lshsec,
                            double *work, int *lwork, int *ierror);
extern void F77_NAME(divgc)(int *nlat, int *nlon, int *isym, int *nt,
                            double *dv, int *idv, int *jdv,
                            double *br, double *bi, int *mdb, int *ndb,
                            double *wshsgc, int *lshsgc,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(vrtec)(int *nlat, int *nlon, int *isym, int *nt,
                            double *vort, int *ivrt, int *jvrt,
                            double *cr, double *ci, int *mdc, int *ndc,
                            double *wshsec, int *lshsec,
                            double *work, int *lwork, int *ierror);
extern void F77_NAME(vrtgc)(int *nlat, int *nlon, int *isym, int *nt,
                            double *vort, int *ivrt, int *jvrt,
                            double *cr, double *ci, int *mdc, int *ndc,
                            double *wshsgc, int *lshsgc,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(idivec)(int *nlat, int *nlon, int *isym, int *nt,
                             double *v, double *w, int *idvw, int *jdvw,
                             double *a, double *b, int *mdab, int *ndab,
                             double *wvhsec, int *lvhsec,
                             double *work, int *lwork,
                             double *pertrb, int *ierror);
extern void F77_NAME(idivgc)(int *nlat, int *nlon, int *isym, int *nt,
                             double *v, double *w, int *idvw, int *jdvw,
                             double *a, double *b, int *mdab, int *ndab,
                             double *wvhsgc, int *lvhsgc,
                             double *work, int *lwork,
                             double *pertrb, int *ierror);

extern void F77_NAME(ivrtec)(int *nlat, int *nlon, int *isym, int *nt,
                             double *v, double *w, int *idvw, int *jdvw,
                             double *a, double *b, int *mdab, int *ndab,
                             double *wvhsec, int *lvhsec,
                             double *work, int *lwork,
                             double *pertrb, int *ierror);
extern void F77_NAME(ivrtgc)(int *nlat, int *nlon, int *isym, int *nt,
                             double *v, double *w, int *idvw, int *jdvw,
                             double *a, double *b, int *mdab, int *ndab,
                             double *wvhsgc, int *lvhsgc,
                             double *work, int *lwork,
                             double *pertrb, int *ierror);

extern void F77_NAME(idvtec)(int *nlat, int *nlon, int *isym, int *nt,
                             double *v, double *w, int *idvw, int *jdvw,
                             double *ad, double *bd, double *av, double *bv,
                             int *mdab, int *ndab,
                             double *wvhsec, int *lvhsec,
                             double *work, int *lwork,
                             double *pertbd, double *pertbv, int *ierror);
extern void F77_NAME(idvtgc)(int *nlat, int *nlon, int *isym, int *nt,
                             double *v, double *w, int *idvw, int *jdvw,
                             double *ad, double *bd, double *av, double *bv,
                             int *mdab, int *ndab,
                             double *wvhsgc, int *lvhsgc,
                             double *work, int *lwork,
                             double *pertbd, double *pertbv, int *ierror);

extern void F77_NAME(sfvpec)(int *nlat, int *nlon, int *isym, int *nt,
                             double *sf, double *vp, int *idv, int *jdv,
                             double *br, double *bi, double *cr, double *ci,
                             int *mdb, int *ndb,
                             double *wshsec, int *lshsec,
                             double *work, int *lwork, int *ierror);
extern void F77_NAME(sfvpgc)(int *nlat, int *nlon, int *isym, int *nt,
                             double *sf, double *vp, int *idv, int *jdv,
                             double *br, double *bi, double *cr, double *ci,
                             int *mdb, int *ndb,
                             double *wshsgc, int *lshsgc,
                             double *work, int *lwork, int *ierror);

extern void F77_NAME(isfvpec)(int *nlat, int *nlon, int *isym, int *nt,
                              double *v, double *w, int *idv, int *jdv,
                              double *as, double *bs, double *av, double *bv,
                              int *mdb, int *ndb,
                              double *wvhsec, int *lvhsec,
                              double *work, int *lwork, int *ierror);
extern void F77_NAME(isfvpgc)(int *nlat, int *nlon, int *isym, int *nt,
                              double *v, double *w, int *idv, int *jdv,
                              double *as, double *bs, double *av, double *bv,
                              int *mdb, int *ndb,
                              double *wvhsgc, int *lvhsgc,
                              double *work, int *lwork, int *ierror);

extern void F77_NAME(gradec)(int *nlat, int *nlon, int *isym, int *nt,
                             double *v, double *w, int *idvw, int *jdvw,
                             double *a, double *b, int *mdab, int *ndab,
                             double *wvhsec, int *lvhsec,
                             double *work, int *lwork, int *ierror);
extern void F77_NAME(gradgc)(int *nlat, int *nlon, int *isym, int *nt,
                             double *v, double *w, int *idvw, int *jdvw,
                             double *a, double *b, int *mdab, int *ndab,
                             double *wvhsgc, int *lvhsgc,
                             double *work, int *lwork, int *ierror);

extern void F77_NAME(igradec)(int *nlat, int *nlon, int *isym, int *nt,
                              double *sf, int *isf, int *jsf,
                              double *br, double *bi, int *mdb, int *ndb,
                              double *wshsec, int *lshsec,
                              double *work, int *lwork, int *ierror);
extern void F77_NAME(igradgc)(int *nlat, int *nlon, int *isym, int *nt,
                              double *sf, int *isf, int *jsf,
                              double *br, double *bi, int *mdb, int *ndb,
                              double *wshsgc, int *lshsgc,
                              double *work, int *lwork, int *ierror);

extern void F77_NAME(slapec)(int *nlat, int *nlon, int *isym, int *nt,
                             double *slap, int *ids, int *jds,
                             double *a, double *b, int *mdab, int *ndab,
                             double *wshsec, int *lshsec,
                             double *work, int *lwork, int *ierror);
extern void F77_NAME(slapgc)(int *nlat, int *nlon, int *isym, int *nt,
                             double *slap, int *ids, int *jds,
                             double *a, double *b, int *mdab, int *ndab,
                             double *wshsgc, int *lshsgc,
                             double *work, int *lwork, int *ierror);

extern void F77_NAME(islapec)(int *nlat, int *nlon, int *isym, int *nt,
                              double *xlmbda, double *sf, int *ids, int *jds,
                              double *a, double *b, int *mdab, int *ndab,
                              double *wshsec, int *lshsec,
                              double *work, int *lwork,
                              double *pertrb, int *ierror);
extern void F77_NAME(islapgc)(int *nlat, int *nlon, int *isym, int *nt,
                              double *xlmbda, double *sf, int *ids, int *jds,
                              double *a, double *b, int *mdab, int *ndab,
                              double *wshsgc, int *lshsgc,
                              double *work, int *lwork,
                              double *pertrb, int *ierror);

extern void F77_NAME(vlapec)(int *nlat, int *nlon, int *ityp, int *nt,
                             double *vlap, double *wlap, int *idvw, int *jdvw,
                             double *br, double *bi, double *cr, double *ci,
                             int *mdbc, int *ndbc,
                             double *wvhsec, int *lvhsec,
                             double *work, int *lwork, int *ierror);
extern void F77_NAME(vlapgc)(int *nlat, int *nlon, int *ityp, int *nt,
                             double *vlap, double *wlap, int *idvw, int *jdvw,
                             double *br, double *bi, double *cr, double *ci,
                             int *mdbc, int *ndbc,
                             double *wvhsgc, int *lvhsgc,
                             double *work, int *lwork, int *ierror);

extern void F77_NAME(ivlapec)(int *nlat, int *nlon, int *ityp, int *nt,
                              double *v, double *w, int *idvw, int *jdvw,
                              double *br, double *bi, double *cr, double *ci,
                              int *mdbc, int *ndbc,
                              double *wvhsec, int *lvhsec,
                              double *work, int *lwork, int *ierror);
extern void F77_NAME(ivlapgc)(int *nlat, int *nlon, int *ityp, int *nt,
                              double *v, double *w, int *idvw, int *jdvw,
                              double *br, double *bi, double *cr, double *ci,
                              int *mdbc, int *ndbc,
                              double *wvhsgc, int *lvhsgc,
                              double *work, int *lwork, int *ierror);

extern void F77_NAME(vtseci)(int *nlat, int *nlon, double *wvts, int *lwvts,
                             double *dwork, int *ldwork, int *ierror);
extern void F77_NAME(vtsec)(int *nlat, int *nlon, int *ityp, int *nt,
                            double *vt, double *wt, int *idvw, int *jdvw,
                            double *br, double *bi, double *cr, double *ci,
                            int *mdab, int *ndab,
                            double *wvts, int *lwvts,
                            double *work, int *lwork, int *ierror);

extern void F77_NAME(vtsgci)(int *nlat, int *nlon, double *wvts, int *lwvts,
                             double *dwork, int *ldwork, int *ierror);
extern void F77_NAME(vtsgc)(int *nlat, int *nlon, int *ityp, int *nt,
                            double *vt, double *wt, int *idvw, int *jdvw,
                            double *br, double *bi, double *cr, double *ci,
                            int *mdab, int *ndab,
                            double *wvts, int *lwvts,
                            double *work, int *lwork, int *ierror);

/* Workspace size helpers */
static int imin(int a, int b) { return a < b ? a : b; }
static int imax(int a, int b) { return a > b ? a : b; }

static int calc_lshaec(int nlat, int nlon) {
  int l1 = imin(nlat, (nlon + 2) / 2);
  int l2 = (nlat + 1) / 2;
  return 2 * nlat * l2 + 3 * ((l1 * (2 * nlat - l1 + 1)) / 2) + nlon + 64;
}

static int calc_lshagc(int nlat, int nlon) {
  int l1 = imin(nlat, (nlon + 2) / 2);
  int l2 = (nlat + 1) / 2;
  return nlat * (2 * l2 + 3 * l1) + 3 * ((l1 * (2 * nlat - l1 + 1)) / 2) + nlon + 64;
}

static int calc_lvhaec(int nlat, int nlon) {
  int l1 = imin(nlat, (nlon + 2) / 2);
  int l2 = (nlat + 1) / 2;
  int idz = (l1 * (2 * nlat - l1 + 1)) / 2;
  int lzimn = idz * l2;
  int lw1 = 4 * nlat * l2 + 3 * l1 * (2 * nlat - l1 + 2) + nlon + 64;
  int lw2 = 2 * lzimn + nlon + 64;
  return imax(lw1, lw2);
}

static int calc_ldwork(int nlat) {
  return 4 * nlat * (nlat + 4) + 64;
}

static int calc_lwork(int nlat, int nlon, int nt) {
  int l1 = imin(nlat, (nlon + 2) / 2);
  int m6 = imax(6 * nlat, nlon);
  return nlat * (4 * nt * nlon + m6 + 6 * l1 * nt + 12) + 64;
}

static void check_ier(int ierror, const char *routine) {
  if (ierror != 0) {
    Rf_error("SPHEREPACK error %d in subroutine %s", ierror, routine);
  }
}

static void geo2math_scalar_batch(int nlon, int nlat, int nt,
                                  double *sg, double *sm, double *work) {
  int ig = 0;
  size_t slice = (size_t)nlon * (size_t)nlat;
  for (int k = 0; k < nt; k++) {
    F77_CALL(geo2maths)(&ig, &nlon, &nlat, sg + k * slice, sm + k * slice, work);
  }
}

static void math2geo_scalar_batch(int nlat, int nlon, int nt,
                                  double *sm, double *sg, double *work,
                                  double scale) {
  int ig = 0;
  size_t slice = (size_t)nlon * (size_t)nlat;
  for (int k = 0; k < nt; k++) {
    F77_CALL(math2geos)(&ig, &nlat, &nlon, sm + k * slice, sg + k * slice, work);
  }
  if (scale != 1.0) {
    size_t total = slice * (size_t)nt;
    for (size_t i = 0; i < total; i++) {
      sg[i] *= scale;
    }
  }
}

static void geo2math_vector_batch(int nlon, int nlat, int nt,
                                  double *ug, double *vg,
                                  double *vm, double *wm, double *work) {
  int ig = 0;
  size_t slice = (size_t)nlon * (size_t)nlat;
  for (int k = 0; k < nt; k++) {
    F77_CALL(geo2mathv)(&ig, &nlon, &nlat,
                        ug + k * slice, vg + k * slice,
                        vm + k * slice, wm + k * slice, work);
  }
}

static void math2geo_vector_batch(int nlat, int nlon, int nt,
                                  double *vm, double *wm,
                                  double *ug, double *vg, double *work,
                                  double scale) {
  int ig = 0;
  size_t slice = (size_t)nlon * (size_t)nlat;
  for (int k = 0; k < nt; k++) {
    F77_CALL(math2geov)(&ig, &nlat, &nlon,
                        vm + k * slice, wm + k * slice,
                        ug + k * slice, vg + k * slice, work);
  }
  if (scale != 1.0) {
    size_t total = slice * (size_t)nt;
    for (size_t i = 0; i < total; i++) {
      ug[i] *= scale;
      vg[i] *= scale;
    }
  }
}

/* 1. Gaussian quadrature points and weights */
SEXP c_gaqd(SEXP s_nlat) {
  int nlat = Rf_asInteger(s_nlat);
  if (nlat < 1) {
    Rf_error("nlat must be >= 1");
  }
  SEXP s_theta = PROTECT(Rf_allocVector(REALSXP, nlat));
  SEXP s_wts = PROTECT(Rf_allocVector(REALSXP, nlat));
  double dummy_w = 0.0;
  int lwork = 1;
  int ierror = 0;
  F77_CALL(gaqd)(&nlat, REAL(s_theta), REAL(s_wts), &dummy_w, &lwork, &ierror);
  check_ier(ierror, "gaqd");

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("theta"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("wts"));
  SET_VECTOR_ELT(res, 0, s_theta);
  SET_VECTOR_ELT(res, 1, s_wts);
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(4);
  return res;
}

/* 2. Scalar spherical harmonic analysis (grid_type: 0 = fixed, 1 = gaussian) */
SEXP c_sha(SEXP s_g, SEXP s_nlon, SEXP s_nlat, SEXP s_nt, SEXP s_grid_type) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  int isym = 0;
  int mmax = imin(nlat, (nlon + 2) / 2);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  size_t n_spec = (size_t)mmax * (size_t)nlat * (size_t)nt;

  double *gm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));
  geo2math_scalar_batch(nlon, nlat, nt, REAL(s_g), gm, gwork);

  SEXP s_a = PROTECT(Rf_allocVector(REALSXP, n_spec));
  SEXP s_b = PROTECT(Rf_allocVector(REALSXP, n_spec));
  memset(REAL(s_a), 0, n_spec * sizeof(double));
  memset(REAL(s_b), 0, n_spec * sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    int lshaec = calc_lshaec(nlat, nlon);
    double *wshaec = (double *)R_alloc(lshaec, sizeof(double));
    F77_CALL(shaeci)(&nlat, &nlon, wshaec, &lshaec, dwork, &ldwork, &ierror);
    check_ier(ierror, "shaeci");
    F77_CALL(shaec)(&nlat, &nlon, &isym, &nt, gm, &nlat, &nlon,
                    REAL(s_a), REAL(s_b), &mmax, &nlat,
                    wshaec, &lshaec, work, &lwork, &ierror);
    check_ier(ierror, "shaec");
  } else {
    int lshagc = calc_lshagc(nlat, nlon);
    double *wshagc = (double *)R_alloc(lshagc, sizeof(double));
    F77_CALL(shagci)(&nlat, &nlon, wshagc, &lshagc, dwork, &ldwork, &ierror);
    check_ier(ierror, "shagci");
    F77_CALL(shagc)(&nlat, &nlon, &isym, &nt, gm, &nlat, &nlon,
                    REAL(s_a), REAL(s_b), &mmax, &nlat,
                    wshagc, &lshagc, work, &lwork, &ierror);
    check_ier(ierror, "shagc");
  }

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("a"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("b"));
  SET_VECTOR_ELT(res, 0, s_a);
  SET_VECTOR_ELT(res, 1, s_b);
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(4);
  return res;
}

/* 3. Scalar spherical harmonic synthesis (grid_type: 0 = fixed, 1 = gaussian) */
SEXP c_shs(SEXP s_a, SEXP s_b, SEXP s_mdab, SEXP s_nlat, SEXP s_nlon,
           SEXP s_nt, SEXP s_grid_type) {
  int mdab = Rf_asInteger(s_mdab);
  int nlat = Rf_asInteger(s_nlat);
  int nlon = Rf_asInteger(s_nlon);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  int isym = 0;

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *gm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    int lshsec = calc_lshaec(nlat, nlon);
    double *wshsec = (double *)R_alloc(lshsec, sizeof(double));
    F77_CALL(shseci)(&nlat, &nlon, wshsec, &lshsec, dwork, &ldwork, &ierror);
    check_ier(ierror, "shseci");
    F77_CALL(shsec)(&nlat, &nlon, &isym, &nt, gm, &nlat, &nlon,
                    REAL(s_a), REAL(s_b), &mdab, &nlat,
                    wshsec, &lshsec, work, &lwork, &ierror);
    check_ier(ierror, "shsec");
  } else {
    int lshsgc = calc_lshagc(nlat, nlon);
    double *wshsgc = (double *)R_alloc(lshsgc, sizeof(double));
    F77_CALL(shsgci)(&nlat, &nlon, wshsgc, &lshsgc, dwork, &ldwork, &ierror);
    check_ier(ierror, "shsgci");
    F77_CALL(shsgc)(&nlat, &nlon, &isym, &nt, gm, &nlat, &nlon,
                    REAL(s_a), REAL(s_b), &mdab, &nlat,
                    wshsgc, &lshsgc, work, &lwork, &ierror);
    check_ier(ierror, "shsgc");
  }

  SEXP s_g = PROTECT(Rf_allocVector(REALSXP, n_grid));
  math2geo_scalar_batch(nlat, nlon, nt, gm, REAL(s_g), gwork, 1.0);
  UNPROTECT(1);
  return s_g;
}

/* 4. Vector spherical harmonic analysis (grid_type: 0 = fixed, 1 = gaussian) */
SEXP c_vha(SEXP s_u, SEXP s_v, SEXP s_nlon, SEXP s_nlat, SEXP s_nt, SEXP s_grid_type) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  int ityp = 0;
  int mmax = imin(nlat, (nlon + 1) / 2);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  size_t n_spec = (size_t)mmax * (size_t)nlat * (size_t)nt;

  double *vm = (double *)R_alloc(n_grid, sizeof(double));
  double *wm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));
  geo2math_vector_batch(nlon, nlat, nt, REAL(s_u), REAL(s_v), vm, wm, gwork);

  SEXP s_br = PROTECT(Rf_allocVector(REALSXP, n_spec));
  SEXP s_bi = PROTECT(Rf_allocVector(REALSXP, n_spec));
  SEXP s_cr = PROTECT(Rf_allocVector(REALSXP, n_spec));
  SEXP s_ci = PROTECT(Rf_allocVector(REALSXP, n_spec));
  memset(REAL(s_br), 0, n_spec * sizeof(double));
  memset(REAL(s_bi), 0, n_spec * sizeof(double));
  memset(REAL(s_cr), 0, n_spec * sizeof(double));
  memset(REAL(s_ci), 0, n_spec * sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int lvha = calc_lvhaec(nlat, nlon);
  double *wvha = (double *)R_alloc(lvha, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    F77_CALL(vhaeci)(&nlat, &nlon, wvha, &lvha, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhaeci");
    F77_CALL(vhaec)(&nlat, &nlon, &ityp, &nt, vm, wm, &nlat, &nlon,
                    REAL(s_br), REAL(s_bi), REAL(s_cr), REAL(s_ci),
                    &mmax, &nlat, wvha, &lvha, work, &lwork, &ierror);
    check_ier(ierror, "vhaec");
  } else {
    F77_CALL(vhagci)(&nlat, &nlon, wvha, &lvha, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhagci");
    F77_CALL(vhagc)(&nlat, &nlon, &ityp, &nt, vm, wm, &nlat, &nlon,
                    REAL(s_br), REAL(s_bi), REAL(s_cr), REAL(s_ci),
                    &mmax, &nlat, wvha, &lvha, work, &lwork, &ierror);
    check_ier(ierror, "vhagc");
  }

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(nms, 0, Rf_mkChar("br"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("bi"));
  SET_STRING_ELT(nms, 2, Rf_mkChar("cr"));
  SET_STRING_ELT(nms, 3, Rf_mkChar("ci"));
  SET_VECTOR_ELT(res, 0, s_br);
  SET_VECTOR_ELT(res, 1, s_bi);
  SET_VECTOR_ELT(res, 2, s_cr);
  SET_VECTOR_ELT(res, 3, s_ci);
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(6);
  return res;
}

/* 5. Vector spherical harmonic synthesis (grid_type: 0 = fixed, 1 = gaussian) */
SEXP c_vhs(SEXP s_br, SEXP s_bi, SEXP s_cr, SEXP s_ci,
           SEXP s_mdab, SEXP s_nlat, SEXP s_nlon, SEXP s_nt, SEXP s_grid_type) {
  int mdab = Rf_asInteger(s_mdab);
  int nlat = Rf_asInteger(s_nlat);
  int nlon = Rf_asInteger(s_nlon);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  int ityp = 0;

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *vm = (double *)R_alloc(n_grid, sizeof(double));
  double *wm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int lvhs = calc_lvhaec(nlat, nlon);
  double *wvhs = (double *)R_alloc(lvhs, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    F77_CALL(vhseci)(&nlat, &nlon, wvhs, &lvhs, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhseci");
    F77_CALL(vhsec)(&nlat, &nlon, &ityp, &nt, vm, wm, &nlat, &nlon,
                    REAL(s_br), REAL(s_bi), REAL(s_cr), REAL(s_ci),
                    &mdab, &nlat, wvhs, &lvhs, work, &lwork, &ierror);
    check_ier(ierror, "vhsec");
  } else {
    F77_CALL(vhsgci)(&nlat, &nlon, wvhs, &lvhs, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhsgci");
    F77_CALL(vhsgc)(&nlat, &nlon, &ityp, &nt, vm, wm, &nlat, &nlon,
                    REAL(s_br), REAL(s_bi), REAL(s_cr), REAL(s_ci),
                    &mdab, &nlat, wvhs, &lvhs, work, &lwork, &ierror);
    check_ier(ierror, "vhsgc");
  }

  SEXP s_u = PROTECT(Rf_allocVector(REALSXP, n_grid));
  SEXP s_v = PROTECT(Rf_allocVector(REALSXP, n_grid));
  math2geo_vector_batch(nlat, nlon, nt, vm, wm, REAL(s_u), REAL(s_v), gwork, 1.0);

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("u"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("v"));
  SET_VECTOR_ELT(res, 0, s_u);
  SET_VECTOR_ELT(res, 1, s_v);
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(4);
  return res;
}

/* Internal helper: analyze vector (u, v) into (br, bi, cr, ci) */
static void analyze_vector(int nlon, int nlat, int nt, int grid_type,
                           double *ug, double *vg,
                           double **p_br, double **p_bi,
                           double **p_cr, double **p_ci,
                           int *p_mmax) {
  int ityp = 0;
  int mmax = imin(nlat, (nlon + 1) / 2);
  *p_mmax = mmax;
  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  size_t n_spec = (size_t)mmax * (size_t)nlat * (size_t)nt;

  double *vm = (double *)R_alloc(n_grid, sizeof(double));
  double *wm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));
  geo2math_vector_batch(nlon, nlat, nt, ug, vg, vm, wm, gwork);

  double *br = (double *)R_alloc(n_spec, sizeof(double));
  double *bi = (double *)R_alloc(n_spec, sizeof(double));
  double *cr = (double *)R_alloc(n_spec, sizeof(double));
  double *ci = (double *)R_alloc(n_spec, sizeof(double));
  memset(br, 0, n_spec * sizeof(double));
  memset(bi, 0, n_spec * sizeof(double));
  memset(cr, 0, n_spec * sizeof(double));
  memset(ci, 0, n_spec * sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int lvha = calc_lvhaec(nlat, nlon);
  double *wvha = (double *)R_alloc(lvha, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    F77_CALL(vhaeci)(&nlat, &nlon, wvha, &lvha, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhaeci");
    F77_CALL(vhaec)(&nlat, &nlon, &ityp, &nt, vm, wm, &nlat, &nlon,
                    br, bi, cr, ci, &mmax, &nlat, wvha, &lvha, work, &lwork, &ierror);
    check_ier(ierror, "vhaec");
  } else {
    F77_CALL(vhagci)(&nlat, &nlon, wvha, &lvha, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhagci");
    F77_CALL(vhagc)(&nlat, &nlon, &ityp, &nt, vm, wm, &nlat, &nlon,
                    br, bi, cr, ci, &mmax, &nlat, wvha, &lvha, work, &lwork, &ierror);
    check_ier(ierror, "vhagc");
  }

  *p_br = br;
  *p_bi = bi;
  *p_cr = cr;
  *p_ci = ci;
}

/* Internal helper: analyze scalar g into (a, b) */
static void analyze_scalar(int nlon, int nlat, int nt, int grid_type,
                           double *sg,
                           double **p_a, double **p_b, int *p_mmax) {
  int isym = 0;
  int mmax = imin(nlat, (nlon + 2) / 2);
  *p_mmax = mmax;
  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  size_t n_spec = (size_t)mmax * (size_t)nlat * (size_t)nt;

  double *sm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));
  geo2math_scalar_batch(nlon, nlat, nt, sg, sm, gwork);

  double *a = (double *)R_alloc(n_spec, sizeof(double));
  double *b = (double *)R_alloc(n_spec, sizeof(double));
  memset(a, 0, n_spec * sizeof(double));
  memset(b, 0, n_spec * sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    int lshaec = calc_lshaec(nlat, nlon);
    double *wshaec = (double *)R_alloc(lshaec, sizeof(double));
    F77_CALL(shaeci)(&nlat, &nlon, wshaec, &lshaec, dwork, &ldwork, &ierror);
    check_ier(ierror, "shaeci");
    F77_CALL(shaec)(&nlat, &nlon, &isym, &nt, sm, &nlat, &nlon,
                    a, b, &mmax, &nlat, wshaec, &lshaec, work, &lwork, &ierror);
    check_ier(ierror, "shaec");
  } else {
    int lshagc = calc_lshagc(nlat, nlon);
    double *wshagc = (double *)R_alloc(lshagc, sizeof(double));
    F77_CALL(shagci)(&nlat, &nlon, wshagc, &lshagc, dwork, &ldwork, &ierror);
    check_ier(ierror, "shagci");
    F77_CALL(shagc)(&nlat, &nlon, &isym, &nt, sm, &nlat, &nlon,
                    a, b, &mmax, &nlat, wshagc, &lshagc, work, &lwork, &ierror);
    check_ier(ierror, "shagc");
  }

  *p_a = a;
  *p_b = b;
}

/* 6. uv2dv (mode = 1), uv2vr (mode = 2), uv2vrdv (mode = 3) */
SEXP c_uv2vrdv_op(SEXP s_u, SEXP s_v, SEXP s_nlon, SEXP s_nlat, SEXP s_nt,
                  SEXP s_grid_type, SEXP s_mode, SEXP s_radius) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  int mode = Rf_asInteger(s_mode); /* 1 = dv, 2 = vr, 3 = both */
  double radius = Rf_asReal(s_radius);
  int isym = 0;

  double *br, *bi, *cr, *ci;
  int mmax;
  analyze_vector(nlon, nlat, nt, grid_type, REAL(s_u), REAL(s_v),
                 &br, &bi, &cr, &ci, &mmax);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));
  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int ierror = 0;

  int lshs = (grid_type == 0) ? calc_lshaec(nlat, nlon) : calc_lshagc(nlat, nlon);
  double *wshs = (double *)R_alloc(lshs, sizeof(double));
  if (grid_type == 0) {
    F77_CALL(shseci)(&nlat, &nlon, wshs, &lshs, dwork, &ldwork, &ierror);
    check_ier(ierror, "shseci");
  } else {
    F77_CALL(shsgci)(&nlat, &nlon, wshs, &lshs, dwork, &ldwork, &ierror);
    check_ier(ierror, "shsgci");
  }

  SEXP s_dv = R_NilValue, s_vr = R_NilValue;
  int nprot = 0;

  if (mode == 1 || mode == 3) {
    double *dvm = (double *)R_alloc(n_grid, sizeof(double));
    if (grid_type == 0) {
      F77_CALL(divec)(&nlat, &nlon, &isym, &nt, dvm, &nlat, &nlon,
                      br, bi, &mmax, &nlat, wshs, &lshs, work, &lwork, &ierror);
      check_ier(ierror, "divec");
    } else {
      F77_CALL(divgc)(&nlat, &nlon, &isym, &nt, dvm, &nlat, &nlon,
                      br, bi, &mmax, &nlat, wshs, &lshs, work, &lwork, &ierror);
      check_ier(ierror, "divgc");
    }
    s_dv = PROTECT(Rf_allocVector(REALSXP, n_grid));
    nprot++;
    math2geo_scalar_batch(nlat, nlon, nt, dvm, REAL(s_dv), gwork, 1.0 / radius);
  }

  if (mode == 2 || mode == 3) {
    double *vrm = (double *)R_alloc(n_grid, sizeof(double));
    if (grid_type == 0) {
      F77_CALL(vrtec)(&nlat, &nlon, &isym, &nt, vrm, &nlat, &nlon,
                      cr, ci, &mmax, &nlat, wshs, &lshs, work, &lwork, &ierror);
      check_ier(ierror, "vrtec");
    } else {
      F77_CALL(vrtgc)(&nlat, &nlon, &isym, &nt, vrm, &nlat, &nlon,
                      cr, ci, &mmax, &nlat, wshs, &lshs, work, &lwork, &ierror);
      check_ier(ierror, "vrtgc");
    }
    s_vr = PROTECT(Rf_allocVector(REALSXP, n_grid));
    nprot++;
    math2geo_scalar_batch(nlat, nlon, nt, vrm, REAL(s_vr), gwork, 1.0 / radius);
  }

  if (mode == 1) {
    UNPROTECT(nprot);
    return s_dv;
  } else if (mode == 2) {
    UNPROTECT(nprot);
    return s_vr;
  } else {
    SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(nms, 0, Rf_mkChar("vr"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("dv"));
    SET_VECTOR_ELT(res, 0, s_vr);
    SET_VECTOR_ELT(res, 1, s_dv);
    Rf_setAttrib(res, R_NamesSymbol, nms);
    UNPROTECT(nprot + 2);
    return res;
  }
}

/* 7. dv2uv (mode = 1), vr2uv (mode = 2), vrdv2uv (mode = 3) */
SEXP c_vrdv2uv_op(SEXP s_vr, SEXP s_dv, SEXP s_nlon, SEXP s_nlat, SEXP s_nt,
                  SEXP s_grid_type, SEXP s_mode, SEXP s_radius) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  int mode = Rf_asInteger(s_mode); /* 1 = dv2uv, 2 = vr2uv, 3 = vrdv2uv */
  double radius = Rf_asReal(s_radius);
  int isym = 0;

  double *ad = NULL, *bd = NULL, *av = NULL, *bv = NULL;
  int mmax = 0;
  if (mode == 1 || mode == 3) {
    analyze_scalar(nlon, nlat, nt, grid_type, REAL(s_dv), &ad, &bd, &mmax);
  }
  if (mode == 2 || mode == 3) {
    analyze_scalar(nlon, nlat, nt, grid_type, REAL(s_vr), &av, &bv, &mmax);
  }

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *vm = (double *)R_alloc(n_grid, sizeof(double));
  double *wm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));
  double *pertbd = (double *)R_alloc(nt, sizeof(double));
  double *pertbv = (double *)R_alloc(nt, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int lvhs = calc_lvhaec(nlat, nlon);
  double *wvhs = (double *)R_alloc(lvhs, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    F77_CALL(vhseci)(&nlat, &nlon, wvhs, &lvhs, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhseci");
    if (mode == 1) {
      F77_CALL(idivec)(&nlat, &nlon, &isym, &nt, vm, wm, &nlat, &nlon,
                       ad, bd, &mmax, &nlat, wvhs, &lvhs, work, &lwork, pertbd, &ierror);
      check_ier(ierror, "idivec");
    } else if (mode == 2) {
      F77_CALL(ivrtec)(&nlat, &nlon, &isym, &nt, vm, wm, &nlat, &nlon,
                       av, bv, &mmax, &nlat, wvhs, &lvhs, work, &lwork, pertbv, &ierror);
      check_ier(ierror, "ivrtec");
    } else {
      F77_CALL(idvtec)(&nlat, &nlon, &isym, &nt, vm, wm, &nlat, &nlon,
                       ad, bd, av, bv, &mmax, &nlat, wvhs, &lvhs, work, &lwork,
                       pertbd, pertbv, &ierror);
      check_ier(ierror, "idvtec");
    }
  } else {
    F77_CALL(vhsgci)(&nlat, &nlon, wvhs, &lvhs, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhsgci");
    if (mode == 1) {
      F77_CALL(idivgc)(&nlat, &nlon, &isym, &nt, vm, wm, &nlat, &nlon,
                       ad, bd, &mmax, &nlat, wvhs, &lvhs, work, &lwork, pertbd, &ierror);
      check_ier(ierror, "idivgc");
    } else if (mode == 2) {
      F77_CALL(ivrtgc)(&nlat, &nlon, &isym, &nt, vm, wm, &nlat, &nlon,
                       av, bv, &mmax, &nlat, wvhs, &lvhs, work, &lwork, pertbv, &ierror);
      check_ier(ierror, "ivrtgc");
    } else {
      F77_CALL(idvtgc)(&nlat, &nlon, &isym, &nt, vm, wm, &nlat, &nlon,
                       ad, bd, av, bv, &mmax, &nlat, wvhs, &lvhs, work, &lwork,
                       pertbd, pertbv, &ierror);
      check_ier(ierror, "idvtgc");
    }
  }

  SEXP s_u = PROTECT(Rf_allocVector(REALSXP, n_grid));
  SEXP s_v = PROTECT(Rf_allocVector(REALSXP, n_grid));
  math2geo_vector_batch(nlat, nlon, nt, vm, wm, REAL(s_u), REAL(s_v), gwork, radius);

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("u"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("v"));
  SET_VECTOR_ELT(res, 0, s_u);
  SET_VECTOR_ELT(res, 1, s_v);
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(4);
  return res;
}

/* 8. uv2sfvp */
SEXP c_uv2sfvp(SEXP s_u, SEXP s_v, SEXP s_nlon, SEXP s_nlat, SEXP s_nt,
               SEXP s_grid_type, SEXP s_radius) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  double radius = Rf_asReal(s_radius);
  int isym = 0;

  double *br, *bi, *cr, *ci;
  int mmax;
  analyze_vector(nlon, nlat, nt, grid_type, REAL(s_u), REAL(s_v),
                 &br, &bi, &cr, &ci, &mmax);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *sfm = (double *)R_alloc(n_grid, sizeof(double));
  double *vpm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    int lshsec = calc_lshaec(nlat, nlon);
    double *wshsec = (double *)R_alloc(lshsec, sizeof(double));
    F77_CALL(shseci)(&nlat, &nlon, wshsec, &lshsec, dwork, &ldwork, &ierror);
    check_ier(ierror, "shseci");
    F77_CALL(sfvpec)(&nlat, &nlon, &isym, &nt, sfm, vpm, &nlat, &nlon,
                     br, bi, cr, ci, &mmax, &nlat, wshsec, &lshsec, work, &lwork, &ierror);
    check_ier(ierror, "sfvpec");
  } else {
    int lshsgc = calc_lshagc(nlat, nlon);
    double *wshsgc = (double *)R_alloc(lshsgc, sizeof(double));
    F77_CALL(shsgci)(&nlat, &nlon, wshsgc, &lshsgc, dwork, &ldwork, &ierror);
    check_ier(ierror, "shsgci");
    F77_CALL(sfvpgc)(&nlat, &nlon, &isym, &nt, sfm, vpm, &nlat, &nlon,
                     br, bi, cr, ci, &mmax, &nlat, wshsgc, &lshsgc, work, &lwork, &ierror);
    check_ier(ierror, "sfvpgc");
  }

  SEXP s_sf = PROTECT(Rf_allocVector(REALSXP, n_grid));
  SEXP s_vp = PROTECT(Rf_allocVector(REALSXP, n_grid));
  math2geo_scalar_batch(nlat, nlon, nt, sfm, REAL(s_sf), gwork, radius);
  math2geo_scalar_batch(nlat, nlon, nt, vpm, REAL(s_vp), gwork, radius);

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("sf"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("vp"));
  SET_VECTOR_ELT(res, 0, s_sf);
  SET_VECTOR_ELT(res, 1, s_vp);
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(4);
  return res;
}

/* 9. sfvp2uv */
SEXP c_sfvp2uv(SEXP s_sf, SEXP s_vp, SEXP s_nlon, SEXP s_nlat, SEXP s_nt,
               SEXP s_grid_type, SEXP s_radius) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  double radius = Rf_asReal(s_radius);
  int isym = 0;

  double *as, *bs, *av, *bv;
  int mmax;
  analyze_scalar(nlon, nlat, nt, grid_type, REAL(s_sf), &as, &bs, &mmax);
  analyze_scalar(nlon, nlat, nt, grid_type, REAL(s_vp), &av, &bv, &mmax);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *vm = (double *)R_alloc(n_grid, sizeof(double));
  double *wm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int lvhs = calc_lvhaec(nlat, nlon);
  double *wvhs = (double *)R_alloc(lvhs, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    F77_CALL(vhseci)(&nlat, &nlon, wvhs, &lvhs, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhseci");
    F77_CALL(isfvpec)(&nlat, &nlon, &isym, &nt, vm, wm, &nlat, &nlon,
                      as, bs, av, bv, &mmax, &nlat, wvhs, &lvhs, work, &lwork, &ierror);
    check_ier(ierror, "isfvpec");
  } else {
    F77_CALL(vhsgci)(&nlat, &nlon, wvhs, &lvhs, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhsgci");
    F77_CALL(isfvpgc)(&nlat, &nlon, &isym, &nt, vm, wm, &nlat, &nlon,
                      as, bs, av, bv, &mmax, &nlat, wvhs, &lvhs, work, &lwork, &ierror);
    check_ier(ierror, "isfvpgc");
  }

  SEXP s_u = PROTECT(Rf_allocVector(REALSXP, n_grid));
  SEXP s_v = PROTECT(Rf_allocVector(REALSXP, n_grid));
  math2geo_vector_batch(nlat, nlon, nt, vm, wm, REAL(s_u), REAL(s_v), gwork, 1.0 / radius);

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("u"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("v"));
  SET_VECTOR_ELT(res, 0, s_u);
  SET_VECTOR_ELT(res, 1, s_v);
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(4);
  return res;
}

/* 10. gradsf / gradsg */
SEXP c_grads(SEXP s_z, SEXP s_nlon, SEXP s_nlat, SEXP s_nt,
             SEXP s_grid_type, SEXP s_radius) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  double radius = Rf_asReal(s_radius);
  int isym = 0;

  double *a, *b;
  int mmax;
  analyze_scalar(nlon, nlat, nt, grid_type, REAL(s_z), &a, &b, &mmax);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *vm = (double *)R_alloc(n_grid, sizeof(double));
  double *wm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int lvhs = calc_lvhaec(nlat, nlon);
  double *wvhs = (double *)R_alloc(lvhs, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    F77_CALL(vhseci)(&nlat, &nlon, wvhs, &lvhs, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhseci");
    F77_CALL(gradec)(&nlat, &nlon, &isym, &nt, vm, wm, &nlat, &nlon,
                     a, b, &mmax, &nlat, wvhs, &lvhs, work, &lwork, &ierror);
    check_ier(ierror, "gradec");
  } else {
    F77_CALL(vhsgci)(&nlat, &nlon, wvhs, &lvhs, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhsgci");
    F77_CALL(gradgc)(&nlat, &nlon, &isym, &nt, vm, wm, &nlat, &nlon,
                     a, b, &mmax, &nlat, wvhs, &lvhs, work, &lwork, &ierror);
    check_ier(ierror, "gradgc");
  }

  SEXP s_gzx = PROTECT(Rf_allocVector(REALSXP, n_grid));
  SEXP s_gzy = PROTECT(Rf_allocVector(REALSXP, n_grid));
  math2geo_vector_batch(nlat, nlon, nt, vm, wm, REAL(s_gzx), REAL(s_gzy), gwork, 1.0 / radius);

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("gzx"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("gzy"));
  SET_VECTOR_ELT(res, 0, s_gzx);
  SET_VECTOR_ELT(res, 1, s_gzy);
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(4);
  return res;
}

/* 11. igradsf / igradsg */
SEXP c_igrads(SEXP s_gzx, SEXP s_gzy, SEXP s_nlon, SEXP s_nlat, SEXP s_nt,
              SEXP s_grid_type, SEXP s_radius) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  double radius = Rf_asReal(s_radius);
  int isym = 0;

  double *br, *bi, *cr, *ci;
  int mmax;
  analyze_vector(nlon, nlat, nt, grid_type, REAL(s_gzx), REAL(s_gzy),
                 &br, &bi, &cr, &ci, &mmax);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *zm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    int lshsec = calc_lshaec(nlat, nlon);
    double *wshsec = (double *)R_alloc(lshsec, sizeof(double));
    F77_CALL(shseci)(&nlat, &nlon, wshsec, &lshsec, dwork, &ldwork, &ierror);
    check_ier(ierror, "shseci");
    F77_CALL(igradec)(&nlat, &nlon, &isym, &nt, zm, &nlat, &nlon,
                      br, bi, &mmax, &nlat, wshsec, &lshsec, work, &lwork, &ierror);
    check_ier(ierror, "igradec");
  } else {
    int lshsgc = calc_lshagc(nlat, nlon);
    double *wshsgc = (double *)R_alloc(lshsgc, sizeof(double));
    F77_CALL(shsgci)(&nlat, &nlon, wshsgc, &lshsgc, dwork, &ldwork, &ierror);
    check_ier(ierror, "shsgci");
    F77_CALL(igradgc)(&nlat, &nlon, &isym, &nt, zm, &nlat, &nlon,
                      br, bi, &mmax, &nlat, wshsgc, &lshsgc, work, &lwork, &ierror);
    check_ier(ierror, "igradgc");
  }

  SEXP s_z = PROTECT(Rf_allocVector(REALSXP, n_grid));
  math2geo_scalar_batch(nlat, nlon, nt, zm, REAL(s_z), gwork, radius);
  UNPROTECT(1);
  return s_z;
}

/* 12. lapsf / lapsg */
SEXP c_laps(SEXP s_z, SEXP s_nlon, SEXP s_nlat, SEXP s_nt,
            SEXP s_grid_type, SEXP s_radius) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  double radius = Rf_asReal(s_radius);
  int isym = 0;

  double *a, *b;
  int mmax;
  analyze_scalar(nlon, nlat, nt, grid_type, REAL(s_z), &a, &b, &mmax);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *zlapm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int ierror = 0;

  memset(zlapm, 0, n_grid * sizeof(double));
  if (grid_type == 0) {
    int lshsec = calc_lshaec(nlat, nlon);
    double *wshsec = (double *)R_alloc(lshsec, sizeof(double));
    F77_CALL(shseci)(&nlat, &nlon, wshsec, &lshsec, dwork, &ldwork, &ierror);
    check_ier(ierror, "shseci");
    F77_CALL(slapec)(&nlat, &nlon, &isym, &nt, zlapm, &nlat, &nlon,
                     a, b, &mmax, &nlat, wshsec, &lshsec, work, &lwork, &ierror);
    check_ier(ierror, "slapec");
  } else {
    int lshsgc = calc_lshagc(nlat, nlon);
    double *wshsgc = (double *)R_alloc(lshsgc, sizeof(double));
    F77_CALL(shsgci)(&nlat, &nlon, wshsgc, &lshsgc, dwork, &ldwork, &ierror);
    check_ier(ierror, "shsgci");
    F77_CALL(slapgc)(&nlat, &nlon, &isym, &nt, zlapm, &nlat, &nlon,
                     a, b, &mmax, &nlat, wshsgc, &lshsgc, work, &lwork, &ierror);
    check_ier(ierror, "slapgc");
  }

  SEXP s_zlap = PROTECT(Rf_allocVector(REALSXP, n_grid));
  math2geo_scalar_batch(nlat, nlon, nt, zlapm, REAL(s_zlap), gwork, 1.0 / (radius * radius));
  UNPROTECT(1);
  return s_zlap;
}

/* 13. ilapsf / ilapsg */
SEXP c_ilaps(SEXP s_zlap, SEXP s_zlmbda, SEXP s_nlon, SEXP s_nlat, SEXP s_nt,
             SEXP s_grid_type, SEXP s_radius) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  double radius = Rf_asReal(s_radius);
  int isym = 0;

  double *a, *b;
  int mmax;
  analyze_scalar(nlon, nlat, nt, grid_type, REAL(s_zlap), &a, &b, &mmax);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *zm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));
  double *pertrb = (double *)R_alloc(nt, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    int lshsec = calc_lshaec(nlat, nlon);
    double *wshsec = (double *)R_alloc(lshsec, sizeof(double));
    F77_CALL(shseci)(&nlat, &nlon, wshsec, &lshsec, dwork, &ldwork, &ierror);
    check_ier(ierror, "shseci");
    F77_CALL(islapec)(&nlat, &nlon, &isym, &nt, REAL(s_zlmbda), zm, &nlat, &nlon,
                      a, b, &mmax, &nlat, wshsec, &lshsec, work, &lwork, pertrb, &ierror);
    check_ier(ierror, "islapec");
  } else {
    int lshsgc = calc_lshagc(nlat, nlon);
    double *wshsgc = (double *)R_alloc(lshsgc, sizeof(double));
    F77_CALL(shsgci)(&nlat, &nlon, wshsgc, &lshsgc, dwork, &ldwork, &ierror);
    check_ier(ierror, "shsgci");
    F77_CALL(islapgc)(&nlat, &nlon, &isym, &nt, REAL(s_zlmbda), zm, &nlat, &nlon,
                      a, b, &mmax, &nlat, wshsgc, &lshsgc, work, &lwork, pertrb, &ierror);
    check_ier(ierror, "islapgc");
  }

  SEXP s_z = PROTECT(Rf_allocVector(REALSXP, n_grid));
  math2geo_scalar_batch(nlat, nlon, nt, zm, REAL(s_z), gwork, radius * radius);
  UNPROTECT(1);
  return s_z;
}

/* 14. lapvf / lapvg (mode = 0) and ilapvf / ilapvg (mode = 1) */
SEXP c_vlap_op(SEXP s_u, SEXP s_v, SEXP s_nlon, SEXP s_nlat, SEXP s_nt,
               SEXP s_grid_type, SEXP s_inverse, SEXP s_radius) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  int inverse = Rf_asInteger(s_inverse);
  double radius = Rf_asReal(s_radius);
  int ityp = 0;

  double *br, *bi, *cr, *ci;
  int mmax;
  analyze_vector(nlon, nlat, nt, grid_type, REAL(s_u), REAL(s_v),
                 &br, &bi, &cr, &ci, &mmax);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *vm = (double *)R_alloc(n_grid, sizeof(double));
  double *wm = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int lvhs = calc_lvhaec(nlat, nlon);
  double *wvhs = (double *)R_alloc(lvhs, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    F77_CALL(vhseci)(&nlat, &nlon, wvhs, &lvhs, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhseci");
    if (!inverse) {
      F77_CALL(vlapec)(&nlat, &nlon, &ityp, &nt, vm, wm, &nlat, &nlon,
                       br, bi, cr, ci, &mmax, &nlat, wvhs, &lvhs, work, &lwork, &ierror);
      check_ier(ierror, "vlapec");
    } else {
      F77_CALL(ivlapec)(&nlat, &nlon, &ityp, &nt, vm, wm, &nlat, &nlon,
                        br, bi, cr, ci, &mmax, &nlat, wvhs, &lvhs, work, &lwork, &ierror);
      check_ier(ierror, "ivlapec");
    }
  } else {
    F77_CALL(vhsgci)(&nlat, &nlon, wvhs, &lvhs, dwork, &ldwork, &ierror);
    check_ier(ierror, "vhsgci");
    if (!inverse) {
      F77_CALL(vlapgc)(&nlat, &nlon, &ityp, &nt, vm, wm, &nlat, &nlon,
                       br, bi, cr, ci, &mmax, &nlat, wvhs, &lvhs, work, &lwork, &ierror);
      check_ier(ierror, "vlapgc");
    } else {
      F77_CALL(ivlapgc)(&nlat, &nlon, &ityp, &nt, vm, wm, &nlat, &nlon,
                        br, bi, cr, ci, &mmax, &nlat, wvhs, &lvhs, work, &lwork, &ierror);
      check_ier(ierror, "ivlapgc");
    }
  }

  double scale = inverse ? (radius * radius) : (1.0 / (radius * radius));
  SEXP s_out1 = PROTECT(Rf_allocVector(REALSXP, n_grid));
  SEXP s_out2 = PROTECT(Rf_allocVector(REALSXP, n_grid));
  math2geo_vector_batch(nlat, nlon, nt, vm, wm, REAL(s_out1), REAL(s_out2), gwork, scale);

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar(inverse ? "u" : "ulap"));
  SET_STRING_ELT(nms, 1, Rf_mkChar(inverse ? "v" : "vlap"));
  SET_VECTOR_ELT(res, 0, s_out1);
  SET_VECTOR_ELT(res, 1, s_out2);
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(4);
  return res;
}

/* 15. lderuvf / lderuvg */
SEXP c_lderuv(SEXP s_u, SEXP s_v, SEXP s_nlon, SEXP s_nlat, SEXP s_nt,
              SEXP s_grid_type, SEXP s_radius) {
  int nlon = Rf_asInteger(s_nlon);
  int nlat = Rf_asInteger(s_nlat);
  int nt = Rf_asInteger(s_nt);
  int grid_type = Rf_asInteger(s_grid_type);
  double radius = Rf_asReal(s_radius);
  int ityp = 0;

  double *br, *bi, *cr, *ci;
  int mmax;
  analyze_vector(nlon, nlat, nt, grid_type, REAL(s_u), REAL(s_v),
                 &br, &bi, &cr, &ci, &mmax);

  size_t n_grid = (size_t)nlon * (size_t)nlat * (size_t)nt;
  double *vt = (double *)R_alloc(n_grid, sizeof(double));
  double *wt = (double *)R_alloc(n_grid, sizeof(double));
  double *gwork = (double *)R_alloc((size_t)nlon * (size_t)nlat, sizeof(double));

  int ldwork = calc_ldwork(nlat);
  double *dwork = (double *)R_alloc(ldwork, sizeof(double));
  int lwork = calc_lwork(nlat, nlon, nt);
  double *work = (double *)R_alloc(lwork, sizeof(double));
  int lwvts = calc_lvhaec(nlat, nlon);
  double *wvts = (double *)R_alloc(lwvts, sizeof(double));
  int ierror = 0;

  if (grid_type == 0) {
    F77_CALL(vtseci)(&nlat, &nlon, wvts, &lwvts, dwork, &ldwork, &ierror);
    check_ier(ierror, "vtseci");
    F77_CALL(vtsec)(&nlat, &nlon, &ityp, &nt, vt, wt, &nlat, &nlon,
                    br, bi, cr, ci, &mmax, &nlat, wvts, &lwvts, work, &lwork, &ierror);
    check_ier(ierror, "vtsec");
  } else {
    F77_CALL(vtsgci)(&nlat, &nlon, wvts, &lwvts, dwork, &ldwork, &ierror);
    check_ier(ierror, "vtsgci");
    F77_CALL(vtsgc)(&nlat, &nlon, &ityp, &nt, vt, wt, &nlat, &nlon,
                    br, bi, cr, ci, &mmax, &nlat, wvts, &lwvts, work, &lwork, &ierror);
    check_ier(ierror, "vtsgc");
  }

  SEXP s_uy = PROTECT(Rf_allocVector(REALSXP, n_grid));
  SEXP s_vy = PROTECT(Rf_allocVector(REALSXP, n_grid));
  /* Note: d/d(lat) = -d/d(theta), so math2geov output is negated and scaled by 1/radius */
  math2geo_vector_batch(nlat, nlon, nt, vt, wt, REAL(s_uy), REAL(s_vy), gwork, -1.0 / radius);

  SEXP res = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("uy"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("vy"));
  SET_VECTOR_ELT(res, 0, s_uy);
  SET_VECTOR_ELT(res, 1, s_vy);
  Rf_setAttrib(res, R_NamesSymbol, nms);
  UNPROTECT(4);
  return res;
}

static const R_CallMethodDef CallEntries[] = {
  {"c_gaqd",        (DL_FUNC) &c_gaqd,        1},
  {"c_sha",         (DL_FUNC) &c_sha,         5},
  {"c_shs",         (DL_FUNC) &c_shs,         7},
  {"c_vha",         (DL_FUNC) &c_vha,         6},
  {"c_vhs",         (DL_FUNC) &c_vhs,         9},
  {"c_uv2vrdv_op",  (DL_FUNC) &c_uv2vrdv_op,  8},
  {"c_vrdv2uv_op",  (DL_FUNC) &c_vrdv2uv_op,  8},
  {"c_uv2sfvp",     (DL_FUNC) &c_uv2sfvp,     7},
  {"c_sfvp2uv",     (DL_FUNC) &c_sfvp2uv,     7},
  {"c_grads",       (DL_FUNC) &c_grads,       6},
  {"c_igrads",      (DL_FUNC) &c_igrads,      7},
  {"c_laps",        (DL_FUNC) &c_laps,        6},
  {"c_ilaps",       (DL_FUNC) &c_ilaps,       7},
  {"c_vlap_op",     (DL_FUNC) &c_vlap_op,     8},
  {"c_lderuv",      (DL_FUNC) &c_lderuv,      7},
  {NULL, NULL, 0}
};

void R_init_spherepack(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
