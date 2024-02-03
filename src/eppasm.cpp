#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>

#include <R.h>
#include <Rinternals.h>


#define AGE_START 15

#define NG 2
#define pAG 66
#define pDS 2

#define pIDX_FERT 0
#define pAG_FERT 35
#define pIDX_15TO49 0
#define pAG_15TO49  35
#define pIDX_15PLUS 0
#define pAG_15PLUS  66

#define hAG 9
#define hDS 7
#define hTS 3

#define hIDX_FERT 0
#define hAG_FERT 8
#define hIDX_15TO49 0
#define hAG_15TO49  8
#define hIDX_15PLUS 0
#define hAG_15PLUS  9

#define hIDX_CD4_350 2

#define MALE 0
#define FEMALE 1

#define HIVN 0
#define HIVP 1

#define ART0MOS 0
#define ART6MOS 1
#define ART1YR 2

#define ART_STAGE_PROG_RATE 2.0 // HARD CODED: ART stage progression rate

#define EPP_DIRECTINCID 2  // annual direct incidence inputs (as Spectrum)
#define EPP_DIRECTINFECTIONS 3  // annual number of new infections by sex and age
#define EPP_DIRECTINFECTIONS_HTS 4  // annual number of new infections by sex and age but intercalating infections by DT

#define INCIDPOP_15TO49 0 // age range corresponding to incidence input
#define INCIDPOP_15PLUS 1

#define PROJPERIOD_MIDYEAR 0   // mid-year projection period 
#define PROJPERIOD_CALENDAR 1  // calendar-year projection (Spectrum 6.2 update; December 2022)

using namespace boost;


// Function declarations
SEXP getListElement(SEXP list, const char *str);
int checkListElement(SEXP list, const char *str);

extern "C" {

  SEXP checkBoostAsserts(){
#ifndef BOOST_DISABLE_ASSERTS
    Rprintf("BOOST ASSERTS ENABLED\n");
#endif
    return R_NilValue;
  }

  SEXP eppasmC(SEXP s_fp){

    ////////////////////////////////
    ////  set parameter values  ////
    ////////////////////////////////

    using namespace boost;

    // state space dimensions
    SEXP s_ss = getListElement(s_fp, "ss");
    int PROJ_YEARS = *INTEGER(getListElement(s_ss, "PROJ_YEARS"));
    int HIVSTEPS_PER_YEAR = *INTEGER(getListElement(s_ss, "hiv_steps_per_year"));
    double DT = 1.0/HIVSTEPS_PER_YEAR;
    int *hAG_SPAN = INTEGER(getListElement(s_ss, "h.ag.span"));

    int hAG_START[hAG];
    hAG_START[0] = 0;
    for(int ha = 1; ha < hAG; ha++)
      hAG_START[ha] = hAG_START[ha-1] + hAG_SPAN[ha-1];

    int SIM_YEARS = *INTEGER(getListElement(s_fp, "SIM_YEARS"));

    int projection_period_int = *INTEGER(getListElement(s_fp, "projection_period_int"));

    // demographic projection
    multi_array_ref<double, 2> basepop(REAL(getListElement(s_fp, "basepop")), extents[NG][pAG]);
    multi_array_ref<double, 3> Sx(REAL(getListElement(s_fp, "Sx")), extents[PROJ_YEARS][NG][pAG]);
    multi_array_ref<double, 3> netmigr(REAL(getListElement(s_fp, "netmigr")), extents[PROJ_YEARS][NG][pAG]);
    multi_array_ref<double, 2> asfr(REAL(getListElement(s_fp, "asfr")), extents[PROJ_YEARS][pAG_FERT]);

    multi_array_ref<double, 2> entrantpop(REAL(getListElement(s_fp, "entrantpop")), extents[PROJ_YEARS][NG]);

    int bin_popadjust = *INTEGER(getListElement(s_fp, "popadjust"));
    multi_array_ref<double, 3> target_hivn_pop(REAL(getListElement(s_fp, "target_hivn_pop")), extents[PROJ_YEARS][NG][pAG]);
    multi_array_ref<double, 3> target_hivp_pop(REAL(getListElement(s_fp, "target_hivp_pop")), extents[PROJ_YEARS][NG][pAG]);
    multi_array_ref<double, 3> target_hivn_ha(REAL(getListElement(s_fp, "target_hivn_ha")), extents[PROJ_YEARS][NG][hAG]);
    multi_array_ref<double, 3> target_hivpop_ha(REAL(getListElement(s_fp, "target_hivpop_ha")), extents[PROJ_YEARS][NG][hAG]);
    multi_array_ref<double, 3> target_artpop_ha(REAL(getListElement(s_fp, "target_artpop_ha")), extents[PROJ_YEARS][NG][hAG]);

    // disease progression
    multi_array_ref<double, 3> cd4_initdist(REAL(getListElement(s_fp, "cd4_initdist")), extents[NG][hAG][hDS]);
    multi_array_ref<double, 3> cd4_prog(REAL(getListElement(s_fp, "cd4_prog")), extents[NG][hAG][hDS-1]);
    multi_array_ref<double, 3> cd4_mort(REAL(getListElement(s_fp, "cd4_mort")), extents[NG][hAG][hDS]);
    multi_array_ref<double, 4> art_mort(REAL(getListElement(s_fp, "art_mort")), extents[NG][hAG][hDS][hTS]);
    multi_array_ref<double, 2> artmx_timerr(REAL(getListElement(s_fp, "artmx_timerr")), extents[PROJ_YEARS][hTS]);

    // sub-fertility
    multi_array_ref<double, 3> frr_cd4(REAL(getListElement(s_fp, "frr_cd4")), extents[PROJ_YEARS][hAG_FERT][hDS]);
    multi_array_ref<double, 4> frr_art(REAL(getListElement(s_fp, "frr_art")), extents[PROJ_YEARS][hAG_FERT][hDS][hTS]);

    // ART inputs
    int t_ART_start = *INTEGER(getListElement(s_fp, "tARTstart")) - 1; // -1 for 0-based indexing in C vs. 1-based in R
    multi_array_ref<double, 2> artnum15plus(REAL(getListElement(s_fp, "art15plus_num")), extents[PROJ_YEARS][NG]);
    multi_array_ref<int, 2> art15plus_isperc(LOGICAL(getListElement(s_fp, "art15plus_isperc")), extents[PROJ_YEARS][NG]);

    int *artcd4elig_idx = INTEGER(getListElement(s_fp, "artcd4elig_idx"));  // NOTE: 1-based indexing
    double *specpop_percelig = REAL(getListElement(s_fp, "specpop_percelig"));
    double *pw_artelig = REAL(getListElement(s_fp, "pw_artelig"));
    double who34percelig = *REAL(getListElement(s_fp, "who34percelig"));

    int bin_art_dropout_recover_cd4 = *INTEGER(getListElement(s_fp, "art_dropout_recover_cd4"));
    double *art_dropout = REAL(getListElement(s_fp, "art_dropout"));
    double *median_cd4init = REAL(getListElement(s_fp, "median_cd4init"));

    int *med_cd4init_cat = INTEGER(getListElement(s_fp, "med_cd4init_cat"));
    int *med_cd4init_input = INTEGER(getListElement(s_fp, "med_cd4init_input"));

    int art_alloc_method = *INTEGER(getListElement(s_fp, "art_alloc_method"));
    double art_alloc_mxweight = *REAL(getListElement(s_fp, "art_alloc_mxweight"));
    int scale_cd4_mort = *INTEGER(getListElement(s_fp, "scale_cd4_mort"));


    // diagnosis model
    int t_hts_start = *INTEGER(getListElement(s_fp, "t_hts_start")) - 1; // -1 for 0-based indexing in C vs. 1-based in R
    multi_array_ref<double, 4> hts_rate(REAL(getListElement(s_fp, "hts_rate")), extents[PROJ_YEARS][pDS][NG][hAG]);
    multi_array_ref<double, 5> diagn_rate(REAL(getListElement(s_fp, "diagn_rate")), extents[PROJ_YEARS][4][NG][hAG][hDS]);

    // incidence model
    int eppmod = *INTEGER(getListElement(s_fp, "eppmodInt"));
    double *incrr_sex = nullptr;
    double *a_incrr_age = nullptr;
    double *incidinput = nullptr;
    int pIDX_INCIDPOP, pAG_INCIDPOP;
    double *a_infections = nullptr;
    if(eppmod == EPP_DIRECTINCID){
      incrr_sex = REAL(getListElement(s_fp, "incrr_sex"));
      a_incrr_age = REAL(getListElement(s_fp, "incrr_age"));
      incidinput = REAL(getListElement(s_fp, "incidinput"));
      pIDX_INCIDPOP = 0;
      if(*INTEGER(getListElement(s_fp, "incidpopage")) == INCIDPOP_15TO49) {
        pAG_INCIDPOP = pAG_15TO49;
      } else {
        pAG_INCIDPOP = pAG_15PLUS;
      }
    } else if (eppmod == EPP_DIRECTINFECTIONS ||
               eppmod == EPP_DIRECTINFECTIONS_HTS) {
      a_infections = REAL(getListElement(s_fp, "infections"));
    }
    multi_array_ref<double, 3> incrr_age(a_incrr_age, extents[PROJ_YEARS][NG][pAG]);
    multi_array_ref<double, 3> input_infections(a_infections, extents[PROJ_YEARS][NG][pAG]);


    multi_array_ref<double, 2> entrantprev(REAL(getListElement(s_fp, "entrantprev")), extents[PROJ_YEARS][NG]);
    multi_array_ref<double, 2> entrantartcov(REAL(getListElement(s_fp, "entrantartcov")), extents[PROJ_YEARS][NG]);

    multi_array_ref<double, 3> paedsurv_cd4dist(REAL(getListElement(s_fp, "paedsurv_cd4dist")), extents[PROJ_YEARS][NG][hDS]);
    multi_array_ref<double, 4> paedsurv_artcd4dist(REAL(getListElement(s_fp, "paedsurv_artcd4dist")), extents[PROJ_YEARS][NG][hDS][hTS]);

    // initialize output
    SEXP s_pop = PROTECT(allocVector(REALSXP, pAG * NG * pDS * PROJ_YEARS));
    SEXP s_pop_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_pop_dim)[0] = pAG;
    INTEGER(s_pop_dim)[1] = NG;
    INTEGER(s_pop_dim)[2] = pDS;
    INTEGER(s_pop_dim)[3] = PROJ_YEARS;
    setAttrib(s_pop, R_DimSymbol, s_pop_dim);
    memset(REAL(s_pop), 0, length(s_pop)*sizeof(double));

    SEXP s_hivpop = PROTECT(allocVector(REALSXP, hDS * hAG * NG * PROJ_YEARS));
    SEXP s_hivpop_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_hivpop_dim)[0] = hDS;
    INTEGER(s_hivpop_dim)[1] = hAG;
    INTEGER(s_hivpop_dim)[2] = NG;
    INTEGER(s_hivpop_dim)[3] = PROJ_YEARS;
    setAttrib(s_hivpop, R_DimSymbol, s_hivpop_dim);
    setAttrib(s_pop, install("hivpop"), s_hivpop);
    multi_array_ref<double, 4> hivpop(REAL(s_hivpop), extents[PROJ_YEARS][NG][hAG][hDS]);
    memset(REAL(s_hivpop), 0, length(s_hivpop)*sizeof(double));

    SEXP s_testnegpop = PROTECT(allocVector(REALSXP, hAG * NG * pDS * PROJ_YEARS));
    SEXP s_testnegpop_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_testnegpop_dim)[0] = hAG;
    INTEGER(s_testnegpop_dim)[1] = NG;
    INTEGER(s_testnegpop_dim)[2] = pDS;
    INTEGER(s_testnegpop_dim)[3] = PROJ_YEARS;
    setAttrib(s_testnegpop, R_DimSymbol, s_testnegpop_dim);
    setAttrib(s_pop, install("testnegpop"), s_testnegpop);
    multi_array_ref<double, 4> testnegpop(REAL(s_testnegpop), extents[PROJ_YEARS][pDS][NG][hAG]);
    memset(REAL(s_testnegpop), 0, length(s_testnegpop)*sizeof(double));

    SEXP s_diagnpop = PROTECT(allocVector(REALSXP, hDS * hAG * NG * PROJ_YEARS));
    SEXP s_diagnpop_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_diagnpop_dim)[0] = hDS;
    INTEGER(s_diagnpop_dim)[1] = hAG;
    INTEGER(s_diagnpop_dim)[2] = NG;
    INTEGER(s_diagnpop_dim)[3] = PROJ_YEARS;
    setAttrib(s_diagnpop, R_DimSymbol, s_diagnpop_dim);
    setAttrib(s_pop, install("diagnpop"), s_diagnpop);
    multi_array_ref<double, 4> diagnpop(REAL(s_diagnpop), extents[PROJ_YEARS][NG][hAG][hDS]);
    memset(REAL(s_diagnpop), 0, length(s_diagnpop)*sizeof(double));

    SEXP s_artpop = PROTECT(allocVector(REALSXP, hTS * hDS * hAG * NG * PROJ_YEARS));
    SEXP s_artpop_dim = PROTECT(allocVector(INTSXP, 5));
    INTEGER(s_artpop_dim)[0] = hTS;
    INTEGER(s_artpop_dim)[1] = hDS;
    INTEGER(s_artpop_dim)[2] = hAG;
    INTEGER(s_artpop_dim)[3] = NG;
    INTEGER(s_artpop_dim)[4] = PROJ_YEARS;
    setAttrib(s_artpop, R_DimSymbol, s_artpop_dim);
    setAttrib(s_pop, install("artpop"), s_artpop);
    multi_array_ref<double, 5> artpop(REAL(s_artpop), extents[PROJ_YEARS][NG][hAG][hDS][hTS]);
    memset(REAL(s_artpop), 0, length(s_artpop)*sizeof(double));

    SEXP s_infections = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_infections_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_infections_dim)[0] = pAG;
    INTEGER(s_infections_dim)[1] = NG;
    INTEGER(s_infections_dim)[2] = PROJ_YEARS;
    setAttrib(s_infections, R_DimSymbol, s_infections_dim);
    setAttrib(s_pop, install("infections"), s_infections);
    multi_array_ref<double, 3> infections(REAL(s_infections), extents[PROJ_YEARS][NG][pAG]);
    memset(REAL(s_infections), 0, length(s_infections)*sizeof(double));

    SEXP s_hivdeaths = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_hivdeaths_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_hivdeaths_dim)[0] = pAG;
    INTEGER(s_hivdeaths_dim)[1] = NG;
    INTEGER(s_hivdeaths_dim)[2] = PROJ_YEARS;
    setAttrib(s_hivdeaths, R_DimSymbol, s_hivdeaths_dim);
    setAttrib(s_pop, install("hivdeaths"), s_hivdeaths);
    multi_array_ref<double, 3> hivdeaths(REAL(s_hivdeaths), extents[PROJ_YEARS][NG][pAG]);
    memset(REAL(s_hivdeaths), 0, length(s_hivdeaths)*sizeof(double));

    SEXP s_hivpopdeaths = PROTECT(allocVector(REALSXP, hDS * hAG * NG * PROJ_YEARS));
    SEXP s_hivpopdeaths_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_hivpopdeaths_dim)[0] = hDS;
    INTEGER(s_hivpopdeaths_dim)[1] = hAG;
    INTEGER(s_hivpopdeaths_dim)[2] = NG;
    INTEGER(s_hivpopdeaths_dim)[3] = PROJ_YEARS;
    setAttrib(s_hivpopdeaths, R_DimSymbol, s_hivpopdeaths_dim);
    setAttrib(s_pop, install("hivpopdeaths"), s_hivpopdeaths);
    multi_array_ref<double, 4> hivpopdeaths(REAL(s_hivpopdeaths), extents[PROJ_YEARS][NG][hAG][hDS]);
    memset(REAL(s_hivpopdeaths), 0, length(s_hivpopdeaths)*sizeof(double));

    SEXP s_artpopdeaths = PROTECT(allocVector(REALSXP, hTS * hDS * hAG * NG * PROJ_YEARS));
    SEXP s_artpopdeaths_dim = PROTECT(allocVector(INTSXP, 5));
    INTEGER(s_artpopdeaths_dim)[0] = hTS;
    INTEGER(s_artpopdeaths_dim)[1] = hDS;
    INTEGER(s_artpopdeaths_dim)[2] = hAG;
    INTEGER(s_artpopdeaths_dim)[3] = NG;
    INTEGER(s_artpopdeaths_dim)[4] = PROJ_YEARS;
    setAttrib(s_artpopdeaths, R_DimSymbol, s_artpopdeaths_dim);
    setAttrib(s_pop, install("artpopdeaths"), s_artpopdeaths);
    multi_array_ref<double, 5> artpopdeaths(REAL(s_artpopdeaths), extents[PROJ_YEARS][NG][hAG][hDS][hTS]);
    memset(REAL(s_artpopdeaths), 0, length(s_artpopdeaths)*sizeof(double));

    SEXP s_natdeaths = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_natdeaths_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_natdeaths_dim)[0] = pAG;
    INTEGER(s_natdeaths_dim)[1] = NG;
    INTEGER(s_natdeaths_dim)[2] = PROJ_YEARS;
    setAttrib(s_natdeaths, R_DimSymbol, s_natdeaths_dim);
    setAttrib(s_pop, install("natdeaths"), s_natdeaths);
    multi_array_ref<double, 3> natdeaths(REAL(s_natdeaths), extents[PROJ_YEARS][NG][pAG]);
    memset(REAL(s_natdeaths), 0, length(s_natdeaths)*sizeof(double));

    // 0: negative, never tested
    // 1: negative, previously tested
    // 2: positive, never tested
    // 3: positive, previously tested negative
    // 4: positive, diagnosed and untreated
    // 5: positive, on ART
    SEXP s_hivtests = PROTECT(allocVector(REALSXP, hAG * NG * 6 * PROJ_YEARS));
    SEXP s_hivtests_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_hivtests_dim)[0] = hAG;
    INTEGER(s_hivtests_dim)[1] = NG;
    INTEGER(s_hivtests_dim)[2] = 6;
    INTEGER(s_hivtests_dim)[3] = PROJ_YEARS;
    setAttrib(s_hivtests, R_DimSymbol, s_hivtests_dim);
    setAttrib(s_pop, install("hivtests"), s_hivtests);
    multi_array_ref<double, 4> hivtests(REAL(s_hivtests), extents[PROJ_YEARS][6][NG][hAG]);
    memset(REAL(s_hivtests), 0, length(s_hivtests)*sizeof(double));

    SEXP s_diagnoses = PROTECT(allocVector(REALSXP, hDS * hAG * NG * PROJ_YEARS));
    SEXP s_diagnoses_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_diagnoses_dim)[0] = hDS;
    INTEGER(s_diagnoses_dim)[1] = hAG;
    INTEGER(s_diagnoses_dim)[2] = NG;
    INTEGER(s_diagnoses_dim)[3] = PROJ_YEARS;
    setAttrib(s_diagnoses, R_DimSymbol, s_diagnoses_dim);
    setAttrib(s_pop, install("diagnoses"), s_diagnoses);
    multi_array_ref<double, 4> diagnoses(REAL(s_diagnoses), extents[PROJ_YEARS][NG][hAG][hDS]);
    memset(REAL(s_diagnoses), 0, length(s_diagnoses)*sizeof(double));

    SEXP s_late_diagnoses = PROTECT(allocVector(REALSXP, hDS * hAG * NG * PROJ_YEARS));
    SEXP s_late_diagnoses_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_late_diagnoses_dim)[0] = hDS;
    INTEGER(s_late_diagnoses_dim)[1] = hAG;
    INTEGER(s_late_diagnoses_dim)[2] = NG;
    INTEGER(s_late_diagnoses_dim)[3] = PROJ_YEARS;
    setAttrib(s_late_diagnoses, R_DimSymbol, s_late_diagnoses_dim);
    setAttrib(s_pop, install("late_diagnoses"), s_late_diagnoses);
    multi_array_ref<double, 4> late_diagnoses(REAL(s_late_diagnoses), extents[PROJ_YEARS][NG][hAG][hDS]);
    memset(REAL(s_late_diagnoses), 0, length(s_late_diagnoses)*sizeof(double));

    SEXP s_artinits = PROTECT(allocVector(REALSXP, hDS * hAG * NG * PROJ_YEARS));
    SEXP s_artinits_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(s_artinits_dim)[0] = hDS;
    INTEGER(s_artinits_dim)[1] = hAG;
    INTEGER(s_artinits_dim)[2] = NG;
    INTEGER(s_artinits_dim)[3] = PROJ_YEARS;
    setAttrib(s_artinits, R_DimSymbol, s_artinits_dim);
    setAttrib(s_pop, install("artinits"), s_artinits);
    multi_array_ref<double, 4> artinits(REAL(s_artinits), extents[PROJ_YEARS][NG][hAG][hDS]);
    memset(REAL(s_artinits), 0, length(s_artinits)*sizeof(double));

    SEXP s_popadjust = PROTECT(allocVector(REALSXP, pAG * NG * PROJ_YEARS));
    SEXP s_popadjust_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(s_popadjust_dim)[0] = pAG;
    INTEGER(s_popadjust_dim)[1] = NG;
    INTEGER(s_popadjust_dim)[2] = PROJ_YEARS;
    setAttrib(s_popadjust, R_DimSymbol, s_popadjust_dim);
    setAttrib(s_pop, install("popadjust"), s_popadjust);
    multi_array_ref<double, 3> popadjust(REAL(s_popadjust), extents[PROJ_YEARS][NG][pAG]);
    memset(REAL(s_popadjust), 0, length(s_popadjust)*sizeof(double));

    SEXP s_prev15to49 = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("prev15to49"), s_prev15to49);
    double *prev15to49 = REAL(s_prev15to49);
    prev15to49[0] = 0.0;

    SEXP s_incid15to49 = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("incid15to49"), s_incid15to49);
    double *incid15to49 = REAL(s_incid15to49);
    memset(incid15to49, 0, length(s_incid15to49)*sizeof(double));

    SEXP s_entrantprev_out = PROTECT(allocVector(REALSXP, PROJ_YEARS));
    setAttrib(s_pop, install("entrantprev"), s_entrantprev_out);
    double *entrantprev_out = REAL(s_entrantprev_out);
    memset(entrantprev_out, 0, length(s_entrantprev_out)*sizeof(double));

    double *hivn15to49 = (double*) R_alloc(PROJ_YEARS, sizeof(double));
    double *hivp15to49 = (double*) R_alloc(PROJ_YEARS, sizeof(double));
    memset(hivn15to49, 0, PROJ_YEARS*sizeof(double));
    memset(hivp15to49, 0, PROJ_YEARS*sizeof(double));


    // initialize population

    // population by single-year age
    // double pop[PROJ_YEARS][pDS][NG][pAG];
    multi_array_ref<double, 4> pop(REAL(s_pop), extents[PROJ_YEARS][pDS][NG][pAG]);
    for(int g = 0; g < NG; g++)
      for(int a = 0; a < pAG; a++){
        pop[0][HIVN][g][a] = basepop[g][a];
        pop[0][HIVP][g][a] = 0.0;
        if(a >= pIDX_15TO49 && a < pIDX_15TO49+pAG_15TO49)
          hivn15to49[0] += basepop[g][a];
      }

    // HIV population with stage stratification
    // double hivpop[PROJ_YEARS][NG][hAG][hDS];
    for(int g = 0; g < NG; g++)
      for(int ha = 0; ha < hAG; ha++)
        for(int hm = 0; hm < hDS; hm++)
          hivpop[0][g][ha][hm] = 0.0;

    int everARTelig_idx = hDS;

    ////////////////////////////////////
    ////  do population projection  ////
    ////////////////////////////////////

    for(int t = 1; t < SIM_YEARS; t++){

      // age the population one year
      for(int m = 0; m < pDS; m++)
        for(int g = 0; g < NG; g++){
          for(int a = 1; a < pAG; a++)
            pop[t][m][g][a] = pop[t-1][m][g][a-1];
          pop[t][m][g][pAG-1] += pop[t-1][m][g][pAG-1]; // open age group
        }

      double hiv_ag_prob[NG][hAG];
      for(int g = 0; g < NG; g++){
        int a = 0;
        for(int ha = 0; ha < (hAG-1); ha++){
          hiv_ag_prob[g][ha] = 0;
          for(int i = 0; i < hAG_SPAN[ha]; i++){
            hiv_ag_prob[g][ha] += pop[t-1][HIVP][g][a];
            a++;
          }
          hiv_ag_prob[g][ha] = (hiv_ag_prob[g][ha] > 0) ? pop[t-1][HIVP][g][a-1] / hiv_ag_prob[g][ha] : 0;
        }
        hiv_ag_prob[g][hAG-1] = 0.0; // no one ages out of the open-ended age group
      }

      double hivn_ag_prob[NG][hAG];
      if(t > t_hts_start) {
        for(int g = 0; g < NG; g++){
          int a = 0;
          for(int ha = 0; ha < (hAG-1); ha++){
            hivn_ag_prob[g][ha] = 0;
            for(int i = 0; i < hAG_SPAN[ha]; i++){
              hivn_ag_prob[g][ha] += pop[t-1][HIVN][g][a];
              a++;
            }
            hivn_ag_prob[g][ha] = (hivn_ag_prob[g][ha] > 0) ? pop[t-1][HIVN][g][a-1] / hivn_ag_prob[g][ha] : 0;
          }
          hivn_ag_prob[g][hAG-1] = 0.0; // no one ages out of the open-ended age group
        }

        for(int g = 0; g < NG; g++)
          for(int ha = 1; ha < hAG; ha++) {
            testnegpop[t][HIVN][g][ha] = (1-hivn_ag_prob[g][ha]) * testnegpop[t-1][HIVN][g][ha] + hivn_ag_prob[g][ha-1]*testnegpop[t-1][HIVN][g][ha-1];
            testnegpop[t][HIVP][g][ha] = (1-hiv_ag_prob[g][ha]) * testnegpop[t-1][HIVP][g][ha] + hiv_ag_prob[g][ha-1]*testnegpop[t-1][HIVP][g][ha-1];
          }
      }

      for(int g = 0; g < NG; g++)
        for(int ha = 1; ha < hAG; ha++)
          for(int hm = 0; hm < hDS; hm++){
            hivpop[t][g][ha][hm] = (1-hiv_ag_prob[g][ha]) * hivpop[t-1][g][ha][hm] + hiv_ag_prob[g][ha-1]*hivpop[t-1][g][ha-1][hm];
            if(t > t_hts_start)
              diagnpop[t][g][ha][hm] = (1-hiv_ag_prob[g][ha]) * diagnpop[t-1][g][ha][hm] + hiv_ag_prob[g][ha-1]*diagnpop[t-1][g][ha-1][hm];
            if(t > t_ART_start)
              for(int hu = 0; hu < hTS; hu++)
                artpop[t][g][ha][hm][hu] = (1-hiv_ag_prob[g][ha]) * artpop[t-1][g][ha][hm][hu] + hiv_ag_prob[g][ha-1]*artpop[t-1][g][ha-1][hm][hu];
          }

      // add lagged births to youngest age group
      for(int g = 0; g < NG; g++){

        double paedsurv_g;
        double entrant_prev;

        entrant_prev = entrantprev[t][g];

        pop[t][HIVN][g][0] =  entrantpop[t-1][g] * (1.0-entrant_prev);
        paedsurv_g = entrantpop[t-1][g] * entrant_prev;

        pop[t][HIVP][g][0] = paedsurv_g;

        entrantprev_out[t] = (pop[t][HIVP][MALE][0] + pop[t][HIVP][FEMALE][0]) / (pop[t][HIVN][MALE][0] + pop[t][HIVN][FEMALE][0] + pop[t][HIVP][MALE][0] + pop[t][HIVP][FEMALE][0]);

        if(t > t_hts_start){
          testnegpop[t][HIVN][g][0] = (1-hivn_ag_prob[g][0]) * testnegpop[t-1][HIVN][g][0];
          testnegpop[t][HIVP][g][0] = (1-hiv_ag_prob[g][0]) * testnegpop[t-1][HIVP][g][0];
        }

        for(int hm = 0; hm < hDS; hm++){
          hivpop[t][g][0][hm] = (1-hiv_ag_prob[g][0]) * hivpop[t-1][g][0][hm] + paedsurv_g * paedsurv_cd4dist[t][g][hm] * (1.0 - entrantartcov[t][g]);
          if(t > t_hts_start)
            diagnpop[t][g][0][hm] = (1-hiv_ag_prob[g][0]) * diagnpop[t-1][g][0][hm];
          if(t > t_ART_start){
            for(int hu = 0; hu < hTS; hu++){
              artpop[t][g][0][hm][hu] = (1-hiv_ag_prob[g][0]) * artpop[t-1][g][0][hm][hu];
              artpop[t][g][0][hm][hu] += paedsurv_g * paedsurv_artcd4dist[t][g][hm][hu] * entrantartcov[t][g];
            }
          }
        }
      }

      // non-HIV mortality and netmigration
      for(int g = 0; g < NG; g++){
        int a = 0;
        for(int ha = 0; ha < hAG; ha++){
          double deathsmig_ha = 0, hivpop_ha = 0;
          double deathsmig_hivn_ha = 0.0, hivnpop_ha = 0.0;
          for(int i = 0; i < hAG_SPAN[ha]; i++){

            hivpop_ha += pop[t][HIVP][g][a];
            hivnpop_ha += pop[t][HIVN][g][a];

            // non-HIV mortality
            double qx = 1.0 - Sx[t][g][a];
            double ndeaths_a = pop[t][HIVN][g][a] * qx;
            pop[t][HIVN][g][a] -= ndeaths_a; // survival HIV- population
            double hdeaths_a = pop[t][HIVP][g][a] * qx;
            pop[t][HIVP][g][a] -= hdeaths_a;   // survival HIV+ population
            natdeaths[t][g][a] = ndeaths_a + hdeaths_a;

	    deathsmig_ha -= hdeaths_a;
	    if (t > t_hts_start) {
              deathsmig_hivn_ha -= ndeaths_a;
	    }

	    if (projection_period_int == PROJPERIOD_MIDYEAR) {
	      // net migration
	      double migrate_a = netmigr[t][g][a] * (1+Sx[t][g][a])/2.0 / (pop[t][HIVN][g][a] + pop[t][HIVP][g][a]);
	      double nmig_a = migrate_a * pop[t][HIVN][g][a];
	      pop[t][HIVN][g][a] += nmig_a;
	      double hmig_a = migrate_a * pop[t][HIVP][g][a];
	      pop[t][HIVP][g][a] += hmig_a;
	      
	      deathsmig_ha += hmig_a;
	      if(t > t_hts_start) {
		deathsmig_hivn_ha += nmig_a;
	      }
	    }
            a++;
	  }

          // migration and deaths for hivpop
          double deathmigrate_ha = hivpop_ha > 0 ? deathsmig_ha / hivpop_ha : 0.0;
	  
          if(t > t_hts_start) {
            double deathmigrate_hivn_ha = hivnpop_ha > 0 ? deathsmig_hivn_ha / hivnpop_ha : 0.0;
            testnegpop[t][HIVN][g][ha] *= 1+deathmigrate_hivn_ha;
            testnegpop[t][HIVP][g][ha] *= 1+deathmigrate_ha;
          }

          for(int hm = 0; hm < hDS; hm++){
            hivpop[t][g][ha][hm] *= 1+deathmigrate_ha;
            if(t > t_hts_start)
              diagnpop[t][g][ha][hm] *= 1+deathmigrate_ha;
            if(t > t_ART_start)
              for(int hu = 0; hu < hTS; hu++)
                artpop[t][g][ha][hm][hu] *= 1+deathmigrate_ha;
          } // loop over hm
        } // loop over ha
      } // loop over g

      // fertility
      double births_by_ha[hAG_FERT];
      memset(births_by_ha, 0, hAG_FERT*sizeof(double));
      for(int m = 0; m < pDS; m++){
        int a = pIDX_FERT;
        for(int ha = hIDX_FERT; ha < hIDX_FERT+hAG_FERT; ha++){
          for(int i = 0; i < hAG_SPAN[ha]; i++){
            births_by_ha[ha-hIDX_FERT] += (pop[t-1][m][FEMALE][a] + pop[t][m][FEMALE][a])/2 * asfr[t][a];
            a++;
          }
        }
      }

      ////////////////////////////////
      ////  HIV model simulation  ////
      ////////////////////////////////

      int cd4elig_idx = artcd4elig_idx[t] - 1; // -1 for 0-based indexing vs. 1-based in R
      int anyelig_idx = (specpop_percelig[t] > 0 || pw_artelig[t] > 0) ? 0 : (who34percelig > 0) ? hIDX_CD4_350 : cd4elig_idx;
      everARTelig_idx = anyelig_idx < everARTelig_idx ? anyelig_idx : everARTelig_idx;

      for(int hts = 0; hts < HIVSTEPS_PER_YEAR; hts++){

        double hivdeaths_ha[NG][hAG];
        memset(hivdeaths_ha, 0, sizeof(double)*NG*hAG);

        // untreated population

        // disease progression and mortality
        double grad[NG][hAG][hDS];
        for(int g = 0; g < NG; g++)
          for(int ha = 0; ha < hAG; ha++){
            for(int hm = 0; hm < hDS; hm++){

              double cd4mx_scale = 1.0;
              if(scale_cd4_mort && t >= t_ART_start && hm >= everARTelig_idx){
                double artpop_hahm = 0.0;
                for(int hu = 0; hu < hTS; hu++)
                  artpop_hahm += artpop[t][g][ha][hm][hu];
                cd4mx_scale = hivpop[t][g][ha][hm] / (hivpop[t][g][ha][hm] + artpop_hahm);
              }

              double deaths = cd4mx_scale * cd4_mort[g][ha][hm] * hivpop[t][g][ha][hm];
              hivdeaths_ha[g][ha] += DT*deaths;
              hivpopdeaths[t][g][ha][hm] += DT*deaths;
              grad[g][ha][hm] = -deaths;
            }
            for(int hm = 1; hm < hDS; hm++){
              grad[g][ha][hm-1] -= cd4_prog[g][ha][hm-1] * hivpop[t][g][ha][hm-1];
              grad[g][ha][hm] += cd4_prog[g][ha][hm-1] * hivpop[t][g][ha][hm-1];
            }
          }

        // Add new infections per time step
        if (eppmod == EPP_DIRECTINFECTIONS_HTS) {
          for (int g = 0; g < NG; g++) {
            int a = 0;
            for (int ha = 0; ha < hAG; ha++) {
              double infections_a, infections_ha = 0.0;
              for (int i = 0; i < hAG_SPAN[ha]; i++) {
                infections_a = input_infections[t][g][a];
                infections_ha += input_infections[t][g][a];
                infections[t][g][a] += DT * infections_a;
		// if(t == 1) {
		//   Rprintf("%d %d %d: %f %s\n", hts, g, a, DT * infections_a, infections_a == 0 ? "true" : "false");
		// }
                pop[t][HIVN][g][a] -= DT * infections_a;
                pop[t][HIVP][g][a] += DT * infections_a;
                a++;
              }

              if (ha < (hIDX_15TO49 + hAG_15TO49)) {
                incid15to49[t] += DT * infections_ha;
              }

              for (int hm = 0; hm < hDS; hm++) {
                grad[g][ha][hm] += infections_ha * cd4_initdist[g][ha][hm];
              }
            }
          }
        } // if (eppmod == EPP_DIRECTINFECTIONS_HTS)

        // HIV testing and diagnosis
        if(t >= t_hts_start) {

          for(int g = 0; g < NG; g++) {
            int a = 0;
            for(int ha = 0; ha < hAG; ha++){

              double hivn_pop_ha = 0.0, infections_ha = 0.0;
              for(int i = 0; i < hAG_SPAN[ha]; i++){
                hivn_pop_ha += pop[t][HIVN][g][a];
                infections_ha += input_infections[t][g][a];
                a++;
              }

              // number of new infections among never tested (proportional to the HIV-neg pop)
              double testneg_infections_ha = infections_ha * (testnegpop[t][0][g][ha] / hivn_pop_ha);

              // Tests among HIV negative population
              hivtests[t][0][g][ha] += DT * hts_rate[t][0][g][ha] * (hivn_pop_ha - testnegpop[t][0][g][ha]);
              hivtests[t][1][g][ha] += DT * hts_rate[t][1][g][ha] * testnegpop[t][0][g][ha];
              testnegpop[t][HIVN][g][ha] += DT * hts_rate[t][HIVN][g][ha] * (hivn_pop_ha - testnegpop[t][HIVN][g][ha]);

              // Remove infections from HIV- testnegpop
              testnegpop[t][HIVN][g][ha] -= DT * testneg_infections_ha ;


              // New diagnoses among HIV positive population

              double grad_diagn[hDS], grad_testneg_hivp = 0.0;

              // Add infections among testnegpop
              grad_testneg_hivp += testneg_infections_ha;

              double hivpop_ha = 0.0, undiagnosed_ha = 0.0;
              for(int hm = 0; hm < hDS; hm++) {
                hivpop_ha += hivpop[t][g][ha][hm];
                undiagnosed_ha += hivpop[t][g][ha][hm] - diagnpop[t][g][ha][hm];
              }

              // Remove share of HIV deaths among testneg population
              grad_testneg_hivp -= hivdeaths_ha[g][ha] / DT * (hivpop_ha > 0 ? testnegpop[t][HIVP][g][ha] / hivpop_ha : 0.0);

              double prop_testneg = undiagnosed_ha > 0 ? testnegpop[t][HIVP][g][ha] / undiagnosed_ha : 0.0;

              for(int hm = 0; hm < hDS; hm++){
                double diagn_naive = diagn_rate[t][0][g][ha][hm] * (1.0 - prop_testneg) * (hivpop[t][g][ha][hm] - diagnpop[t][g][ha][hm]);
                double diagn_testneg = diagn_rate[t][1][g][ha][hm] * prop_testneg * (hivpop[t][g][ha][hm] - diagnpop[t][g][ha][hm]);

                hivtests[t][2][g][ha] += DT * diagn_naive;
                hivtests[t][3][g][ha] += DT * diagn_testneg;
                hivtests[t][4][g][ha] += DT * diagn_rate[t][2][g][ha][hm] * diagnpop[t][g][ha][hm];
                if(t >= t_ART_start){
                  for(int hu = 0; hu < hTS; hu++)
                    hivtests[t][5][g][ha] += DT * diagn_rate[t][3][g][ha][hm] * artpop[t][g][ha][hm][hu];
                }

                grad_testneg_hivp -= diagn_testneg;
                diagnoses[t][g][ha][hm] += DT * (diagn_naive + diagn_testneg);
                grad_diagn[hm] = (diagn_naive + diagn_testneg) - cd4_mort[g][ha][hm] * diagnpop[t][g][ha][hm];
              }
              for(int hm = 1; hm < hDS; hm++){
                grad_diagn[hm-1] -= cd4_prog[g][ha][hm-1] * diagnpop[t][g][ha][hm-1];
                grad_diagn[hm] += cd4_prog[g][ha][hm-1] * diagnpop[t][g][ha][hm-1];
              }

              testnegpop[t][HIVP][g][ha] += DT * grad_testneg_hivp;
              for(int hm = 0; hm < hDS; hm++)
                diagnpop[t][g][ha][hm] += DT * grad_diagn[hm];
            }
          }

        } // end HIV testing and diagnosis


        for(int g = 0; g < NG; g++)
          for(int ha = 0; ha < hAG; ha++)
            for(int hm = 0; hm < hDS; hm++)
              hivpop[t][g][ha][hm] += DT*grad[g][ha][hm];

        // ART progression, mortality, and initiation
        if(t >= t_ART_start){

          // progression and mortality
          for(int g = 0; g < NG; g++)
            for(int ha = 0; ha < hAG; ha++)
              for(int hm = everARTelig_idx; hm < hDS; hm++){
                double gradART[hTS];

                for(int hu = 0; hu < hTS; hu++){
                  double deaths = art_mort[g][ha][hm][hu] * artmx_timerr[t][hu] * artpop[t][g][ha][hm][hu];
                  hivdeaths_ha[g][ha] += DT*deaths;
                  artpopdeaths[t][g][ha][hm][hu] += DT*deaths;
                  gradART[hu] = -deaths;
                }

                gradART[ART0MOS] += -ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART0MOS];
                gradART[ART6MOS] += ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART0MOS] - ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART6MOS];
                gradART[ART1YR] += ART_STAGE_PROG_RATE * artpop[t][g][ha][hm][ART6MOS];

                for(int hu = 0; hu < hTS; hu++)
                  artpop[t][g][ha][hm][hu] += DT*gradART[hu];
              }

          // ART dropout
          if(art_dropout[t] > 0){
            for(int g = 0; g < NG; g++)
              for(int ha = 0; ha < hAG; ha++)
                for(int hm = everARTelig_idx; hm < hDS; hm++)
                  for(int hu = 0; hu < hTS; hu++){

		    double dropout_ts = DT * art_dropout[t] * artpop[t][g][ha][hm][hu];
		    
		    if (bin_art_dropout_recover_cd4 && hu >= 2 && hm >= 1) {
		      // recover people on ART >1 year to one higher CD4 category
		      hivpop[t][g][ha][hm-1] += dropout_ts;
		      if(t >= t_hts_start) {
			diagnpop[t][g][ha][hm-1] += dropout_ts;
		      }
		    } else {
		      hivpop[t][g][ha][hm] += dropout_ts;
		      if(t >= t_hts_start) {
			diagnpop[t][g][ha][hm] += dropout_ts;
		      }
		    }
		    
                    artpop[t][g][ha][hm][hu] -= dropout_ts;
                  }
          }

          // ART initiation
          for(int g = 0; g < NG; g++){

            double artelig_hahm[hAG_15PLUS][hDS], Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig15plus = 0.0;
            for(int ha = hIDX_15PLUS; ha < hAG; ha++){
              for(int hm = everARTelig_idx; hm < hDS; hm++){
                if(hm >= anyelig_idx){
                  double prop_elig = (hm >= cd4elig_idx) ? 1.0 : (hm >= hIDX_CD4_350) ? 1.0 - (1.0-specpop_percelig[t])*(1.0-who34percelig) : specpop_percelig[t];
                  Xartelig_15plus += artelig_hahm[ha-hIDX_15PLUS][hm] = prop_elig * hivpop[t][g][ha][hm] ;
                  expect_mort_artelig15plus += cd4_mort[g][ha][hm] * artelig_hahm[ha-hIDX_15PLUS][hm];
                }
                for(int hu = 0; hu < hTS; hu++)
                  Xart_15plus += artpop[t][g][ha][hm][hu];
              }

              // if pw_artelig, add pregnant women to artelig_hahm population
              if(g == FEMALE && pw_artelig[t] > 0 && ha < hAG_FERT){
                double frr_pop_ha = 0;
                for(int a =  hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
                  frr_pop_ha += pop[t][HIVN][g][a]; // add HIV- population
                for(int hm = 0; hm < hDS; hm++){
                  frr_pop_ha += frr_cd4[t][ha-hIDX_FERT][hm] * hivpop[t][g][ha][hm];
                  for(int hu = 0; hu < hTS; hu++)
                    frr_pop_ha += frr_art[t][ha-hIDX_FERT][hm][hu] * artpop[t][g][ha][hm][hu];
                }
                for(int hm = anyelig_idx; hm < cd4elig_idx; hm++){
                  double pw_elig_hahm = births_by_ha[ha-hIDX_FERT] * frr_cd4[t][ha-hIDX_FERT][hm] * hivpop[t][g][ha][hm] / frr_pop_ha;
                  artelig_hahm[ha-hIDX_15PLUS][hm] += pw_elig_hahm;
                  Xartelig_15plus += pw_elig_hahm;
                  expect_mort_artelig15plus += cd4_mort[g][ha][hm] * pw_elig_hahm;
                }
              }
            } // loop over ha

            // calculate number on ART at end of ts, based on number or percent
            double artnum_hts = 0.0;
            if (projection_period_int == PROJPERIOD_MIDYEAR && DT*(hts+1) < 0.5){
	      
              if(!art15plus_isperc[t-2][g] && !art15plus_isperc[t-1][g]){ // both numbers
                artnum_hts = (0.5-DT*(hts+1))*artnum15plus[t-2][g] + (DT*(hts+1)+0.5)*artnum15plus[t-1][g];
              } else if(art15plus_isperc[t-2][g] && art15plus_isperc[t-1][g]){ // both percentages
                double artcov_hts = (0.5-DT*(hts+1))*artnum15plus[t-2][g] + (DT*(hts+1)+0.5)*artnum15plus[t-1][g];
                artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
              } else if(!art15plus_isperc[t-2][g] && art15plus_isperc[t-1][g]){ // transition from number to percentage
                double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
                double artcov_hts = curr_coverage + (artnum15plus[t-1][g] - curr_coverage) * DT / (0.5-DT*hts);
                artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
              }
            } else {
	      
	      // If the projection period is calendar year (>= Spectrum v6.2), this condition is
 	      // always followed, and it interpolates between end of last year and current year (+ 1.0).
 	      // If projection period was mid-year (<= Spectrum v6.19), the second half of the projection
 	      // year interpolates the first half of the calendar year (e.g. hts 7/10 for 2019 interpolates
 	      // December 2018 to December 2019)

	      double art_interp_w = DT*(hts+1.0);
 	      if (projection_period_int == PROJPERIOD_MIDYEAR) {
 		art_interp_w -= 0.5;
 	      }
	      
              if(!art15plus_isperc[t-1][g] && !art15plus_isperc[t][g]){ // both numbers
		artnum_hts = (1.0 - art_interp_w)*artnum15plus[t-1][g] + art_interp_w*artnum15plus[t][g];
              } else if(art15plus_isperc[t-1][g] && art15plus_isperc[t][g]){ // both percentages
		double artcov_hts = (1.0 - art_interp_w)*artnum15plus[t-1][g] + art_interp_w*artnum15plus[t][g];		
                artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
              } else if(!art15plus_isperc[t-1][g] && art15plus_isperc[t][g]){ // transition from number to percentage
                double curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
		double artcov_hts = curr_coverage + (artnum15plus[t][g] - curr_coverage) * DT / (1.0 - art_interp_w + DT);
                artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
              }
            }

            double artinit_hts = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0.0;

            // median CD4 at initiation inputs
            if(med_cd4init_input[t]){

              const int CD4_LOW_LIM[hDS] = {500, 350, 250, 200, 100, 50, 0};
              const int CD4_UPP_LIM[hDS] = {1000, 500, 350, 250, 200, 100, 50};

              int medcd4_idx = med_cd4init_cat[t] - 1; // -1 for 0-based indexing vs. 1-based in R
              double medcat_propbelow = (median_cd4init[t] - CD4_LOW_LIM[medcd4_idx]) / (CD4_UPP_LIM[medcd4_idx] - CD4_LOW_LIM[medcd4_idx]);

              double elig_below = 0.0, elig_above = 0.0;
              for(int ha = hIDX_15PLUS; ha < hAG; ha++){
                for(int hm = anyelig_idx; hm < medcd4_idx; hm++)
                  elig_above += artelig_hahm[ha-hIDX_15PLUS][hm];
                elig_above += (1.0 - medcat_propbelow) * artelig_hahm[ha-hIDX_15PLUS][medcd4_idx];
                elig_below += medcat_propbelow * artelig_hahm[ha-hIDX_15PLUS][medcd4_idx];
                for(int hm = medcd4_idx+1; hm < hDS; hm++)
                  elig_below += artelig_hahm[ha-hIDX_15PLUS][hm];
              }

              double initprob_below = (elig_below > artinit_hts * 0.5) ? artinit_hts * 0.5 / elig_below : 1.0;
              double initprob_above = (elig_above > artinit_hts * 0.5) ? artinit_hts * 0.5 / elig_above : 1.0;
              double initprob_medcat = initprob_below * medcat_propbelow + initprob_above * (1.0-medcat_propbelow);

              for(int ha = hIDX_15PLUS; ha < hAG; ha++){

                double new_diagn_ha = 0.0, undiagnosed_ha = 0.0;

                if(t >= t_hts_start)
                  for(int hm = 0; hm < hDS; hm++)
                    undiagnosed_ha += hivpop[t][g][ha][hm] - diagnpop[t][g][ha][hm];

                for(int hm = anyelig_idx; hm < hDS; hm++){
                  double artinit_hahm;
                  if(hm < medcd4_idx)
                    artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm] * initprob_above;
                  else if(hm == medcd4_idx)
                    artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm] * initprob_medcat;
                  if(hm > medcd4_idx)
                    artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm] * initprob_below;
                  if(artinit_hahm > hivpop[t][g][ha][hm])
                    artinit_hahm = hivpop[t][g][ha][hm];

                  if(t >= t_hts_start){
                    double new_diagn = artinit_hahm > diagnpop[t][g][ha][hm] ? artinit_hahm - diagnpop[t][g][ha][hm] : 0;

                    new_diagn_ha += new_diagn;

                    diagnpop[t][g][ha][hm] -= artinit_hahm - new_diagn;
                    diagnoses[t][g][ha][hm] += new_diagn;
                    late_diagnoses[t][g][ha][hm] += new_diagn;
                  }

                  hivpop[t][g][ha][hm] -= artinit_hahm;
                  artpop[t][g][ha][hm][ART0MOS] += artinit_hahm;
                  artinits[t][g][ha][hm] += artinit_hahm;
                }

                // Remove share of excess ART initiations from either diagnpop or testnegpop
                if(t >= t_hts_start && new_diagn_ha > 0){

                  double diagn_surplus_ha = 0.0;
                  for(int hm = 0; hm < hDS; hm++)
                    diagn_surplus_ha += diagnpop[t][g][ha][hm];

                  double prop_put_back = new_diagn_ha < diagn_surplus_ha ? new_diagn_ha / diagn_surplus_ha : 1.0;
                  for(int hm = 0; hm < hDS; hm++) {
                    double put_back = prop_put_back * diagnpop[t][g][ha][hm];
                    diagnpop[t][g][ha][hm] -= put_back;
                    diagnoses[t][g][ha][hm] -= put_back;
                    late_diagnoses[t][g][ha][hm] -= put_back;
                  }

                  if(new_diagn_ha > diagn_surplus_ha){
                    new_diagn_ha -= diagn_surplus_ha;

                    double new_among_testneg = new_diagn_ha * (undiagnosed_ha > 0 ? testnegpop[t][HIVP][g][ha] / undiagnosed_ha : 0.0);
                    hivtests[t][2][g][ha] += new_diagn_ha - new_among_testneg; // new diagnoses among never tested
                    hivtests[t][3][g][ha] += new_among_testneg;
                    testnegpop[t][HIVP][g][ha] -= new_among_testneg;
                  }
                }

              } // end loop over ha

            } else if(art_alloc_method == 4) {  // lowest CD4 first

              for(int hm = hDS-1; hm >= anyelig_idx; hm--){
                double artelig_hm = 0;
                for(int ha = hIDX_15PLUS; ha < hAG; ha++)
                  artelig_hm += artelig_hahm[ha-hIDX_15PLUS][hm];
                double init_prop = (artelig_hm == 0 || artinit_hts > artelig_hm) ? 1.0 : artinit_hts / artelig_hm;

                for(int ha = hIDX_15PLUS; ha < hAG; ha++){
                  double artinit_hahm = init_prop * artelig_hahm[ha-hIDX_15PLUS][hm];

                  if(t >= t_hts_start){
                    Rprintf("Lowest CD4 first ART allocation option not yet implemented for HIV testing model\n. STOP WHAT YOU ARE DOING NOW");
                    break;
                  }

                  hivpop[t][g][ha][hm] -= artinit_hahm;
                  artpop[t][g][ha][hm][ART0MOS] += artinit_hahm;
                  artinits[t][g][ha][hm] += artinit_hahm;
                }
                if(init_prop < 1.0)
                  break;
                artinit_hts -= init_prop * artelig_hm;
              }

            } else { // Use mixture of eligibility and expected mortality for initiation distribution

              for(int ha = hIDX_15PLUS; ha < hAG; ha++) {

                double new_diagn_ha = 0.0, undiagnosed_ha = 0.0;

                if(t >= t_hts_start)
                  for(int hm = 0; hm < hDS; hm++)
                    undiagnosed_ha += hivpop[t][g][ha][hm] - diagnpop[t][g][ha][hm];

                for(int hm = anyelig_idx; hm < hDS; hm++){

                  double artinit_hahm = artinit_hts * artelig_hahm[ha-hIDX_15PLUS][hm] * ((1.0 - art_alloc_mxweight)/Xartelig_15plus + art_alloc_mxweight * cd4_mort[g][ha][hm] / expect_mort_artelig15plus);
                  if(Xartelig_15plus == 0)
                    artinit_hahm =  0;
                  if(artinit_hahm > artelig_hahm[ha-hIDX_15PLUS][hm])
                    artinit_hahm = artelig_hahm[ha-hIDX_15PLUS][hm];

                  if(t >= t_hts_start){
                    double new_diagn = artinit_hahm > diagnpop[t][g][ha][hm] ? artinit_hahm - diagnpop[t][g][ha][hm] : 0;

                    new_diagn_ha += new_diagn;
                    diagnpop[t][g][ha][hm] -= artinit_hahm - new_diagn;
                    diagnoses[t][g][ha][hm] += new_diagn;
                    late_diagnoses[t][g][ha][hm] += new_diagn;
                  }


                  hivpop[t][g][ha][hm] -= artinit_hahm;
                  artpop[t][g][ha][hm][ART0MOS] += artinit_hahm;
                  artinits[t][g][ha][hm] += artinit_hahm;
                }

                // Remove share of excess ART initiations from either diagnpop or testnegpop
                if(t >= t_hts_start && new_diagn_ha > 0){

                  double diagn_surplus_ha = 0.0;
                  for(int hm = 0; hm < hDS; hm++)
                    diagn_surplus_ha += diagnpop[t][g][ha][hm];

                  double prop_put_back = new_diagn_ha < diagn_surplus_ha ? new_diagn_ha / diagn_surplus_ha : 1.0;
                  for(int hm = 0; hm < hDS; hm++) {
                    double put_back = prop_put_back * diagnpop[t][g][ha][hm];
                    diagnpop[t][g][ha][hm] -= put_back;
                    diagnoses[t][g][ha][hm] -= put_back;
                    late_diagnoses[t][g][ha][hm] -= put_back;
                  }

                  if(new_diagn_ha > diagn_surplus_ha){
                    new_diagn_ha -= diagn_surplus_ha;

                    double new_among_testneg = new_diagn_ha * (undiagnosed_ha > 0 ? testnegpop[t][HIVP][g][ha] / undiagnosed_ha : 0.0);
                    hivtests[t][2][g][ha] += new_diagn_ha - new_among_testneg; // new diagnoses among never tested
                    hivtests[t][3][g][ha] += new_among_testneg;
                    testnegpop[t][HIVP][g][ha] -= new_among_testneg;
                  }
                }

              } // end loop over ha
            }

          }

        } // end ART initiation

        // remove hivdeaths from pop
        for(int g = 0; g < NG; g++){

          // sum HIV+ population size in each hivpop age group
          double hivpop_ha[hAG];
          int a = 0;
          for(int ha = 0; ha < hAG; ha++){
            hivpop_ha[ha] = 0.0;
            for(int i = 0; i < hAG_SPAN[ha]; i++){
              hivpop_ha[ha] += pop[t][HIVP][g][a];
              a++;
            }
          }

          // remove hivdeaths proportionally to age-distribution within each age group
          a = 0;
          for(int ha = 0; ha < hAG; ha++){
            if(hivpop_ha[ha] > 0){
              double hivqx_ha = hivdeaths_ha[g][ha] / hivpop_ha[ha];
              for(int i = 0; i < hAG_SPAN[ha]; i++){
                hivdeaths[t][g][a] += pop[t][HIVP][g][a] * hivqx_ha;
                pop[t][HIVP][g][a] *= (1.0-hivqx_ha);
                a++;
              }
            } else {
              a += hAG_SPAN[ha];
            }  // end if(pop_ha[ha] > 0)
          }
        }

      } // loop HIVSTEPS_PER_YEAR


      if (eppmod == EPP_DIRECTINCID || eppmod == EPP_DIRECTINFECTIONS) {
        // Calculating new infections once per year (like Spectrum)

        double Xhivn[NG], Xhivn_incagerr[NG];
        double incrate_g[NG];
        if(eppmod == EPP_DIRECTINCID) {

          for(int g = 0; g < NG; g++){
            Xhivn[g] = 0.0;
            Xhivn_incagerr[g] = 0.0;
            for(int a = pIDX_INCIDPOP; a < pIDX_INCIDPOP+pAG_INCIDPOP; a++){
              Xhivn[g] += pop[t-1][HIVN][g][a];
              Xhivn_incagerr[g] += incrr_age[t][g][a] * pop[t-1][HIVN][g][a];
            }
          }
          double incrate_i = incidinput[t];
          incrate_g[MALE] = incrate_i * (Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incrr_sex[t]*Xhivn[FEMALE]);
          incrate_g[FEMALE] = incrate_i * incrr_sex[t]*(Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incrr_sex[t]*Xhivn[FEMALE]);
        }

        for(int g = 0; g < NG; g++){
          int a = 0;
          for(int ha = 0; ha < hAG; ha++){
            double infections_a, infections_ha = 0.0;
            double hivn_pop_ha = 0.0;
            for(int i = 0; i < hAG_SPAN[ha]; i++){
              if(eppmod == EPP_DIRECTINCID)
                infections_ha += infections_a = pop[t-1][HIVN][g][a] * incrate_g[g] * incrr_age[t][g][a] * Xhivn[g] / Xhivn_incagerr[g];
              else if(eppmod == EPP_DIRECTINFECTIONS)
                infections_ha += infections_a = input_infections[t][g][a];
              hivn_pop_ha += pop[t][HIVN][g][a];
              infections[t][g][a] += infections_a;
              pop[t][HIVN][g][a] -= infections_a;
              pop[t][HIVP][g][a] += infections_a;
              a++;
            }
            if(ha < hIDX_15TO49+hAG_15TO49)
              incid15to49[t] += infections_ha;

            if(t >= t_hts_start) {
              double testneg_infections_ha = infections_ha * testnegpop[t][HIVN][g][ha] / hivn_pop_ha;
              testnegpop[t][HIVN][g][ha] -= testneg_infections_ha;
              testnegpop[t][HIVP][g][ha] += testneg_infections_ha;
            }

            // add infections to hivpop
            for(int hm = 0; hm < hDS; hm++)
              hivpop[t][g][ha][hm] += infections_ha * cd4_initdist[g][ha][hm];
          }
        }
      } // if (eppmod == EPP_DIRECTINCID || eppmod == EPP_DIRECTINFECTIONS)

      // Net migration for calendar-year projection option with end-year migration
      if (projection_period_int == PROJPERIOD_CALENDAR) {
	
	for(int g = 0; g < NG; g++){
	  int a = 0;
	  for(int ha = 0; ha < hAG; ha++){
	    double mig_ha = 0, hivpop_ha = 0;
	    double mig_hivn_ha = 0.0, hivnpop_ha = 0.0;
	    
	    for(int i = 0; i < hAG_SPAN[ha]; i++){
	      
	      hivpop_ha += pop[t][HIVP][g][a];
	      hivnpop_ha += pop[t][HIVN][g][a];
	      
	      // net migration
	      double migrate_a = netmigr[t][g][a] / (pop[t][HIVN][g][a] + pop[t][HIVP][g][a]);
	      double nmig_a = migrate_a * pop[t][HIVN][g][a];
	      pop[t][HIVN][g][a] += nmig_a;
	      double hmig_a = migrate_a * pop[t][HIVP][g][a];
	      pop[t][HIVP][g][a] += hmig_a;
	      
	      mig_ha += hmig_a;
	      if (t >= t_hts_start) {
		mig_hivn_ha += nmig_a;
	      }
	      a++;
	    }
	    
	    // migration for hivpop
	    double migrate_ha = hivpop_ha > 0 ? mig_ha / hivpop_ha : 0.0;
	    
	    if(t >= t_hts_start) {
	      double migrate_hivn_ha = hivnpop_ha > 0 ? mig_hivn_ha / hivnpop_ha : 0.0;
	      testnegpop[t][HIVN][g][ha] *= 1+migrate_hivn_ha;
	      testnegpop[t][HIVP][g][ha] *= 1+migrate_ha;
	    }

	    for(int hm = 0; hm < hDS; hm++){
	      hivpop[t][g][ha][hm] *= 1+migrate_ha;
	      if(t >= t_hts_start)
		diagnpop[t][g][ha][hm] *= 1+migrate_ha;
	      if(t >= t_ART_start)
		for(int hu = 0; hu < hTS; hu++)
		  artpop[t][g][ha][hm][hu] *= 1+migrate_ha;
	    } // loop over hm
	  } // loop over ha
	} // loop over g
      } // if (projection_period_int == PROJPERIOD_CALENDAR)      

      // adjust population to match target population
      if(bin_popadjust){
        for(int g = 0; g < NG; g++){
          int a = 0;
          for(int ha = 0; ha < hAG; ha++){
            double hivnpop_ha = 0.0;
            for(int i = 0; i < hAG_SPAN[ha]; i++){

              hivnpop_ha += pop[t][HIVN][g][a];

              pop[t][HIVN][g][a] = target_hivn_pop[t][g][a];
              pop[t][HIVP][g][a] = target_hivp_pop[t][g][a];

              a++;
            }

            double hivpop_ha = 0.0;
            for(int hm = 0; hm < hDS; hm++)
              hivpop_ha += hivpop[t][g][ha][hm];

            double hivpop_adj_ha = hivpop_ha > 0 ? target_hivpop_ha[t][g][ha] / hivpop_ha : 1.0;

            if(t >= t_hts_start) {
              double hivn_adj_ha = hivnpop_ha > 0 ? target_hivn_ha[t][g][ha] / hivnpop_ha : 1.0;
              testnegpop[t][HIVN][g][ha] *= hivn_adj_ha;
              testnegpop[t][HIVP][g][ha] *= hivpop_adj_ha;
            }

            for(int hm = 0; hm < hDS; hm++){
              hivpop[t][g][ha][hm] *= hivpop_adj_ha;
              if(t >= t_hts_start)
                diagnpop[t][g][ha][hm] *= hivpop_adj_ha;
            }

            if(t >= t_ART_start){
              double artpop_ha = 0.0;
              for(int hm = 0; hm < hDS; hm++)
                for(int hu = 0; hu < hTS; hu++)
                  artpop_ha += artpop[t][g][ha][hm][hu];

              double artpop_adj_ha = artpop_ha > 0 ? target_artpop_ha[t][g][ha] / artpop_ha : 1.0;

              for(int hm = 0; hm < hDS; hm++)
                for(int hu = 0; hu < hTS; hu++)
                  artpop[t][g][ha][hm][hu] *= artpop_adj_ha;
            }
          } // loop over ha
        } // loop over g
      } // if(bin_popadjust)

      // prevalence 15 to 49
      for(int g = 0; g < NG; g++)
        for(int a = pIDX_15TO49; a < pIDX_15TO49 + pAG_15TO49; a++){
          hivn15to49[t] += pop[t][HIVN][g][a];
          hivp15to49[t] += pop[t][HIVP][g][a];
        }
      prev15to49[t] = hivp15to49[t]/(hivn15to49[t] + hivp15to49[t]);
      incid15to49[t] /= hivn15to49[t-1];
    }

    UNPROTECT(33);
    return s_pop;
  }
}



SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for ( i = 0; i < length(list); i++ )
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) {
      elmt = VECTOR_ELT(list, i);
      break;
    }

  if ( elmt == R_NilValue )
    error("%s missing from list", str);

  return elmt;
}

int checkListElement(SEXP list, const char *str)
{
  SEXP names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++ )
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0 )
      return 1;

  return 0;
}
