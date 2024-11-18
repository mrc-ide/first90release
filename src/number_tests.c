#include <R.h>
#include <Rinternals.h>

SEXP number_testsC(SEXP s_mod,
		   SEXP s_testnegpop,
		   SEXP s_hivpop,
		   SEXP s_diagnpop,
		   SEXP s_artpop,
		   SEXP s_hts_rate,
		   SEXP s_diagn_rate,
		   SEXP s_haidx,
		   SEXP s_hagspan,
		   SEXP s_sidx,
		   SEXP s_hvidx,
		   SEXP s_yidx,
		   SEXP s_agfirst_idx,
		   SEXP s_h_ag_span,
		   //Mathieu added
		   SEXP s_late_diagnoses){

  int nval = Rf_length(s_haidx);
  SEXP s_pop = PROTECT(Rf_allocVector(REALSXP, nval));
  SEXP s_tests = PROTECT(Rf_allocVector(REALSXP, nval));

  SEXP s_out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(s_out, 0, s_pop);
  SET_VECTOR_ELT(s_out, 1, s_tests);

  int *dm = INTEGER(Rf_getAttrib(s_mod, R_DimSymbol));;
  double *mod = REAL(s_mod);

  int *dtn = INTEGER(Rf_getAttrib(s_testnegpop, R_DimSymbol));;
  double *testnegpop = REAL(s_testnegpop);

  int *dh = INTEGER(Rf_getAttrib(s_hivpop, R_DimSymbol));;
  double *hivpop = REAL(s_hivpop);

  int *ddg = INTEGER(Rf_getAttrib(s_diagnpop, R_DimSymbol));;
  double *diagnpop = REAL(s_diagnpop);

  int *da = INTEGER(Rf_getAttrib(s_artpop, R_DimSymbol));;
  double *artpop = REAL(s_artpop);

  int *dhr = INTEGER(Rf_getAttrib(s_hts_rate, R_DimSymbol));;
  double *hts_rate = REAL(s_hts_rate);

  int *ddr = INTEGER(Rf_getAttrib(s_diagn_rate, R_DimSymbol));;
  double *diagn_rate = REAL(s_diagn_rate);

  //Mathieu added
  int *ldr = INTEGER(Rf_getAttrib(s_late_diagnoses, R_DimSymbol));;
  double *late_diagnoses = REAL(s_late_diagnoses);

  for(int i = 0; i < nval; i++){

    int sidx = INTEGER(s_sidx)[i];
    int s1 = (sidx == 2) ? 1 : 0;
    int s2 = (sidx == 1) ? 1 : 2;

    int hvidx = INTEGER(s_hvidx)[i];
    int yidx = INTEGER(s_yidx)[i] - 1;

    double tests = 0;
    double pop = 0;

    for(int sx = s1; sx < s2; sx++){

      for(int ha = INTEGER(s_haidx)[i] - 1;
	  ha < INTEGER(s_haidx)[i] - 1 + INTEGER(s_hagspan)[i];
	  ha++) {

	// testing among HIV-
	if(hvidx != 2) {

	  double pop_ha = 0.0, tested_ha;
	  int pag = INTEGER(s_agfirst_idx)[ha] - 1;
	  for(int a = 0; a < INTEGER(s_h_ag_span)[ha]; a++)
	    pop_ha += mod[pag + a + dm[0] * (sx + dm[1] * (0 + dm[2] * yidx))];

	  tested_ha = testnegpop[ha + dtn[0] * (sx + dtn[1] * (0 + dtn[2] * yidx))];

	  tests += (pop_ha - tested_ha) * hts_rate[ha + dhr[0] * (sx + dhr[1] * (0 + dtn[2] * yidx))];
	  tests += tested_ha * hts_rate[ha + dhr[0] * (sx + dhr[1] * (1 + dtn[2] * yidx))];
	  pop += pop_ha;
	}

	if(hvidx != 1){

	  double testneg_ha, undiagnosed_ha = 0.0;
	  testneg_ha = testnegpop[ha + dtn[0] * (sx + dtn[1] * (1 + dtn[2] * yidx))];
	  for(int m = 0; m < dh[0]; m++){
	    int idx = m + dh[0] * (ha + dh[1] * (sx + dh[2] * yidx));
	    undiagnosed_ha += hivpop[idx] - diagnpop[idx];
	  }
	  double prop_testneg = testneg_ha / undiagnosed_ha;


	  for(int m = 0; m < dh[0]; m++){
	    int idx = m + dh[0] * (ha + dh[1] * (sx + dh[2] * yidx));
	    double undiagnosed_ha_hm = hivpop[idx] - diagnpop[idx];
	    double artpop_ha_hm = 0.0;
	    for(int u = 0; u < da[0]; u++)
	      artpop_ha_hm += artpop[u + da[0] * (m + da[1] * (ha + da[2] * (sx + da[3] * yidx)))];

	    double naive_ha_hm = undiagnosed_ha_hm * (1 - prop_testneg);
	    double testneg_ha_hm = undiagnosed_ha_hm * prop_testneg;

	    tests += naive_ha_hm * diagn_rate[m + ddr[0] * (ha + ddr[1] * (sx + ddr[2] * (0 + ddr[3] * yidx)))];
	    tests += testneg_ha_hm * diagn_rate[m + ddr[0] * (ha + ddr[1] * (sx + ddr[2] * (1 + ddr[3] * yidx)))];
	    tests += diagnpop[idx] * diagn_rate[m + ddr[0] * (ha + ddr[1] * (sx + ddr[2] * (2 + ddr[3] * yidx)))];
	    tests += artpop_ha_hm * diagn_rate[m + ddr[0] * (ha + ddr[1] * (sx + ddr[2] * (3 + ddr[3] * yidx)))];
      // Mathieu added
      tests += late_diagnoses[m + ddr[0] * (ha + ddr[1] * (sx + ddr[2] * yidx))];

	    pop += hivpop[idx] + artpop_ha_hm;
	  }
	}
      }
    }
    REAL(s_pop)[i] = pop;
    REAL(s_tests)[i] = tests;
  }
  UNPROTECT(3);
  return(s_out);
}
