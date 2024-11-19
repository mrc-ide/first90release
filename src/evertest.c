#include <R.h>
#include <Rinternals.h>

SEXP evertestC(SEXP s_mod,
               SEXP s_testnegpop,
               SEXP s_diagnpop,
               SEXP s_artpop,
               SEXP s_haidx,
               SEXP s_hagspan,
               SEXP s_sidx,
               SEXP s_hvidx,
               SEXP s_yidx,
               SEXP s_agfirst_idx,
               SEXP s_h_ag_span){

  int nval = Rf_length(s_haidx);
  SEXP s_val = PROTECT(Rf_allocVector(REALSXP, nval));

  int *dm = INTEGER(Rf_getAttrib(s_mod, R_DimSymbol));;
  double *mod = REAL(s_mod);

  int *dtn = INTEGER(Rf_getAttrib(s_testnegpop, R_DimSymbol));;
  double *testnegpop = REAL(s_testnegpop);

  int *ddg = INTEGER(Rf_getAttrib(s_diagnpop, R_DimSymbol));;
  double *diagnpop = REAL(s_diagnpop);

  int *da = INTEGER(Rf_getAttrib(s_artpop, R_DimSymbol));;
  double *artpop = REAL(s_artpop);

  for(int i = 0; i < nval; i++){

    int sidx = INTEGER(s_sidx)[i];
    int s1 = (sidx == 2) ? 1 : 0;
    int s2 = (sidx == 1) ? 1 : 2;

    int hvidx = INTEGER(s_hvidx)[i];
    int yidx = INTEGER(s_yidx)[i] - 1;

    double tested = 0;
    double pop = 0;

    for(int sx = s1; sx < s2; sx++){

      for(int ha = INTEGER(s_haidx)[i] - 1;
          ha < INTEGER(s_haidx)[i] - 1 + INTEGER(s_hagspan)[i];
          ha++) {

        if(hvidx != 2)
          tested += testnegpop[ha + dtn[0] * (sx + dtn[1] * (0 + dtn[2] * yidx))];

        if(hvidx != 1){
          tested += testnegpop[ha + dtn[0] * (sx + dtn[1] * (1 + dtn[2] * yidx))];

          for(int m = 0; m < ddg[0]; m++){
            tested += diagnpop[m + ddg[0] * (ha + ddg[1] * (sx + ddg[2] * yidx))];

            for(int u = 0; u < da[0]; u++)
              tested += artpop[u + da[0] * (m + da[1] * (ha + da[2] * (sx + da[3] * yidx)))];
          }

        } // if(hvidx != 1)

        int pag = INTEGER(s_agfirst_idx)[ha] - 1;
        for(int a = 0; a < INTEGER(s_h_ag_span)[ha]; a++){

          if(hvidx != 2)
            pop += mod[pag + a + dm[0] * (sx + dm[1] * (0 + dm[2] * yidx))];

          if(hvidx != 1)
            pop += mod[pag + a + dm[0] * (sx + dm[1] * (1 + dm[2] * yidx))];
        }
      }
    }
    REAL(s_val)[i] = tested / pop;
  }

  UNPROTECT(1);
  return(s_val);
}


SEXP diagnosedC(SEXP s_mod,
		SEXP s_diagnpop,
		SEXP s_artpop,
		SEXP s_haidx,
		SEXP s_hagspan,
		SEXP s_sidx,
		SEXP s_yidx,
		SEXP s_agfirst_idx,
		SEXP s_h_ag_span){

  int nval = Rf_length(s_haidx);
  SEXP s_val = PROTECT(Rf_allocVector(REALSXP, nval));

  int *dm = INTEGER(Rf_getAttrib(s_mod, R_DimSymbol));;
  double *mod = REAL(s_mod);

  int *ddg = INTEGER(Rf_getAttrib(s_diagnpop, R_DimSymbol));;
  double *diagnpop = REAL(s_diagnpop);

  int *da = INTEGER(Rf_getAttrib(s_artpop, R_DimSymbol));;
  double *artpop = REAL(s_artpop);

  for(int i = 0; i < nval; i++){

    int sidx = INTEGER(s_sidx)[i];
    int s1 = (sidx == 2) ? 1 : 0;
    int s2 = (sidx == 1) ? 1 : 2;

    int yidx = INTEGER(s_yidx)[i] - 1;

    double diagnosed = 0;
    double hivpop = 0;

    for(int sx = s1; sx < s2; sx++){

      for(int ha = INTEGER(s_haidx)[i] - 1;
          ha < INTEGER(s_haidx)[i] - 1 + INTEGER(s_hagspan)[i];
          ha++) {

          for(int m = 0; m < ddg[0]; m++){
            diagnosed += diagnpop[m + ddg[0] * (ha + ddg[1] * (sx + ddg[2] * yidx))];

            for(int u = 0; u < da[0]; u++)
              diagnosed += artpop[u + da[0] * (m + da[1] * (ha + da[2] * (sx + da[3] * yidx)))];
          }

	  int pag = INTEGER(s_agfirst_idx)[ha] - 1;
	  for(int a = 0; a < INTEGER(s_h_ag_span)[ha]; a++)
	    hivpop += mod[pag + a + dm[0] * (sx + dm[1] * (1 + dm[2] * yidx))];
      }
    }
    REAL(s_val)[i] = diagnosed / hivpop;
  }

  UNPROTECT(1);
  return(s_val);
}
