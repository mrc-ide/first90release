
#' @export
Rcpp::cppFunction('List prb_dx_one_yr_cpp(List s_fp,
                  IntegerVector year,
                  std::string age,
                  std::string sex,
                  std::string test_ever,
                  double dt) {
                  
                  int n_year = year.size();
                  int ind_age;
                  int ind_sex;
                  int ind_th;
                  List ss = s_fp["ss"];
                  int proj_start = ss["proj_start"];
                  NumericVector infections = s_fp["infections"];
                  NumericVector diagn_rate = s_fp["diagn_rate"];
                  NumericVector cd4_init = s_fp["cd4_initdist"];
                  // age index
                  if (age == "15-24") {
                  ind_age = 0;
                  } else if (age == "25-34") {
                  ind_age = 3;
                  } else if (age == "35-49") {
                  ind_age = 5;
                  } else {
                  ind_age = 8;
                  }
                  // sex index
                  if (sex == "male") {
                  ind_sex = 0 ;
                  } else {
                  ind_sex = 1;
                  }
                  // ever tested index
                  if (test_ever == "never") {
                  ind_th = 0;
                  } else {
                  ind_th = 1;
                  }
                  
                  // average number of years in a CD4 category
                  NumericVector fpcd4 = s_fp["cd4_prog"];
                  IntegerVector ind_cd4 = IntegerVector::create(
                  (6 * 9 * ind_sex) + (6 * ind_age) + 0,
                  (6 * 9 * ind_sex) + (6 * ind_age) + 1,
                  (6 * 9 * ind_sex) + (6 * ind_age) + 2,
                  (6 * 9 * ind_sex) + (6 * ind_age) + 3,
                  (6 * 9 * ind_sex) + (6 * ind_age) + 4,
                  (6 * 9 * ind_sex) + (6 * ind_age) + 5);
                  NumericVector cd4_prg = fpcd4[ind_cd4];
                  
                  // death rates by CD4 category
                  NumericVector fpcd4mort = s_fp["cd4_mort"];
                  IntegerVector ind_nmd15 = IntegerVector::create(
                  (7 * 9 * ind_sex) + (7 * 0) + 0,
                  (7 * 9 * ind_sex) + (7 * 0) + 1,
                  (7 * 9 * ind_sex) + (7 * 0) + 2,
                  (7 * 9 * ind_sex) + (7 * 0) + 3,
                  (7 * 9 * ind_sex) + (7 * 0) + 4,
                  (7 * 9 * ind_sex) + (7 * 0) + 5,
                  (7 * 9 * ind_sex) + (7 * 0) + 6);
                  NumericVector cd4_mort15 = fpcd4mort[ind_nmd15];
                  IntegerVector ind_nmd25 = IntegerVector::create(
                  (7 * 9 * ind_sex) + (7 * 3) + 0,
                  (7 * 9 * ind_sex) + (7 * 3) + 1,
                  (7 * 9 * ind_sex) + (7 * 3) + 2,
                  (7 * 9 * ind_sex) + (7 * 3) + 3,
                  (7 * 9 * ind_sex) + (7 * 3) + 4,
                  (7 * 9 * ind_sex) + (7 * 3) + 5,
                  (7 * 9 * ind_sex) + (7 * 3) + 6);
                  NumericVector cd4_mort25 = fpcd4mort[ind_nmd25];
                  IntegerVector ind_nmd35 = IntegerVector::create(
                  (7 * 9 * ind_sex) + (7 * 5) + 0,
                  (7 * 9 * ind_sex) + (7 * 5) + 1,
                  (7 * 9 * ind_sex) + (7 * 5) + 2,
                  (7 * 9 * ind_sex) + (7 * 5) + 3,
                  (7 * 9 * ind_sex) + (7 * 5) + 4,
                  (7 * 9 * ind_sex) + (7 * 5) + 5,
                  (7 * 9 * ind_sex) + (7 * 5) + 6);
                  NumericVector cd4_mort35 = fpcd4mort[ind_nmd35];
                  IntegerVector ind_nmd50 = IntegerVector::create(
                  (7 * 9 * ind_sex) + (7 * 8) + 0,
                  (7 * 9 * ind_sex) + (7 * 8) + 1,
                  (7 * 9 * ind_sex) + (7 * 8) + 2,
                  (7 * 9 * ind_sex) + (7 * 8) + 3,
                  (7 * 9 * ind_sex) + (7 * 8) + 4,
                  (7 * 9 * ind_sex) + (7 * 8) + 5,
                  (7 * 9 * ind_sex) + (7 * 8) + 6);
                  NumericVector cd4_mort50 = fpcd4mort[ind_nmd50];
                  // initial cd4 cell count distribution
                  double cd4_init_1 = cd4_init[(7 * 9 * ind_sex) + (7 * ind_age) + 0];
                  double cd4_init_2 = cd4_init[(7 * 9 * ind_sex) + (7 * ind_age) + 1];
                  
                  // for probability of being Dx or dead before certain time (expressed in dt)
                  int tind1yr = nearbyint(1 / dt);

                  // to store the results
                  NumericVector prb1yr(n_year);

                  // we calculate for each year using a loop
                  for (int n = 0; n < n_year; ++n) {
                  int ind_yr = year[n] - proj_start;
                  
                  // what is the average age of incident HIV infection
                  NumericVector age_15 = NumericVector::create(1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
                  NumericVector inf_15 = infections[seq(66 * 2 * ind_yr + (ind_sex * 66), 66 * 2 * ind_yr + (ind_sex * 66) + 9)];
                  NumericVector nw15 = age_15 * inf_15;
                  double n15 = sum(nw15);
                  double avg_inf15 = 14 + n15 / sum(inf_15);
                  NumericVector age_25 = NumericVector::create(11, 12, 13, 14, 15, 16, 17, 18, 19, 20);
                  NumericVector inf_25 = infections[seq(66 * 2 * ind_yr + (ind_sex * 66) + 10, 66 * 2 * ind_yr + (ind_sex * 66) + 19)];
                  NumericVector nw25 = age_25 * inf_25;
                  double n25 = sum(nw25);
                  double avg_inf25 = 14 + n25 / sum(inf_25);
                  NumericVector age_35 = NumericVector::create(21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35);
                  NumericVector inf_35 = infections[seq(66 * 2 * ind_yr + (ind_sex * 66) + 20, 66 * 2 * ind_yr + (ind_sex * 66) + 34)];
                  NumericVector nw35 = age_35 * inf_35;
                  double n35 = sum(nw35);
                  double avg_inf35 = 14 + n35 / sum(inf_35);
                  
                  // testing rates by CD4 category
                  IntegerVector ind_nm15 = IntegerVector::create(
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + 0,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + 1,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + 2,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + 3,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + 4,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + 5,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + 6);
                  NumericVector nmt15 = diagn_rate[ind_nm15];
                  IntegerVector ind_nm25 = IntegerVector::create(
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 3) + 0,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 3) + 1,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 3) + 2,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 3) + 3,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 3) + 4,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 3) + 5,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 3) + 6);
                  NumericVector nmt25 = diagn_rate[ind_nm25];
                  IntegerVector ind_nm35 = IntegerVector::create(
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 5) + 0,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 5) + 1,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 5) + 2,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 5) + 3,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 5) + 4,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 5) + 5,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 5) + 6);
                  NumericVector nmt35 = diagn_rate[ind_nm35];
                  IntegerVector ind_nm50 = IntegerVector::create(
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 8) + 0,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 8) + 1,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 8) + 2,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 8) + 3,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 8) + 4,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 8) + 5,
                  (7 * 9 * 2 * 4 * ind_yr) + (7 * 9 * 2 * ind_th) + (7 * 9 * ind_sex) + (7 * 8) + 6);
                  NumericVector nmt50 = diagn_rate[ind_nm50];
                  
                  // Creating vector for time and testing rates
                  IntegerVector time_int_ = seq(1, 35 / dt);
                  NumericVector time_int = as<NumericVector>(time_int_) * dt;
                  int vec_l = time_int.size();
                  NumericVector time_t(vec_l);
                  NumericVector test_rate_t(vec_l * 7);
                  NumericVector mort_rate_t(vec_l * 7);
                  
                  if (age == "15-24") {
                  int ind15 = std::round((25 - avg_inf15) / dt);
                  int ind25 = std::round(ind15 + 10 / dt);
                  int ind35 = std::round(ind25 + 15 / dt);
                  for (int j = 0; j < vec_l; ++j) {
                  IntegerVector ind_j = seq(j * 7 + 0, j * 7 + 6);
                  if (j < ind15) {
                  test_rate_t[ind_j] = nmt15;
                  mort_rate_t[ind_j] = cd4_mort15;
                  } else if (j >= ind15 && j < ind25) {
                  test_rate_t[ind_j] = nmt25;
                  mort_rate_t[ind_j] = cd4_mort25;
                  } else if (j >= ind25 && j < ind35) {
                  test_rate_t[ind_j] = nmt35;
                  mort_rate_t[ind_j] = cd4_mort35;
                  } else {
                  test_rate_t[ind_j] = nmt50;
                  mort_rate_t[ind_j] = cd4_mort50;
                  }}
                  } else if (age == "25-34") {
                  int ind25 = std::round((35 - avg_inf25) / dt);
                  int ind35 = std::round(ind25 + 15 / dt);
                  for (int j = 0; j < vec_l; ++j) {
                  IntegerVector ind_j = seq(j * 7 + 0, j * 7 + 6);
                  if (j < ind25) {
                  test_rate_t[ind_j] = nmt25;
                  mort_rate_t[ind_j] = cd4_mort25;
                  } else if (j >= ind25 && j < ind35) {
                  test_rate_t[ind_j] = nmt35;
                  mort_rate_t[ind_j] = cd4_mort35;
                  } else {
                  test_rate_t[ind_j] = nmt50;
                  mort_rate_t[ind_j] = cd4_mort50;
                  }}
                  } else if (age == "35-49") {
                  int ind35 = std::round((50 - avg_inf35) / dt);
                  for (int j = 0; j < vec_l; ++j) {
                  IntegerVector ind_j = seq(j * 7 + 0, j * 7 + 6);
                  if (j < ind35) {
                  test_rate_t[ind_j] = nmt35;
                  mort_rate_t[ind_j] = cd4_mort35;
                  } else {
                  test_rate_t[ind_j] = nmt50;
                  mort_rate_t[ind_j] = cd4_mort50;
                  }}
                  } else {
                  for (int j = 0; j < vec_l; ++j) {
                  IntegerVector ind_j = seq(j * 7 + 0, j * 7 + 6);
                  test_rate_t[ind_j] = nmt50;
                  mort_rate_t[ind_j] = cd4_mort50;
                  }}

                  // To speed up, we only do 2 year.
                  IntegerVector time_int_prb_one_year_ = seq(1, 2 / dt);
                  NumericVector time_int_prb_one_year = as<NumericVector>(time_int_prb_one_year_) * dt;
                  int vec_l_prb_one_year = time_int_prb_one_year.size();

                  // we initialize the model for competing risk
                  NumericVector X(7 * vec_l_prb_one_year);
                  NumericVector nb_dx_tot(vec_l_prb_one_year);
                  X[0] = cd4_init_1;
                  X[1] = cd4_init_2;
                  // loop for competing risk of cd4 progression, HIV-death, and diagnosis
                  for (int i = 1; i < vec_l_prb_one_year; ++i) {
                  int t_i = i - 1;
                  IntegerVector ind_t = seq(i * 7 + 0, i * 7 + 6);
                  
                  X[i * 7 + 0] = X[t_i * 7 + 0] + dt * ( - (cd4_prg[0] + mort_rate_t[i * 7 + 0] + test_rate_t[i * 7 + 0]) * X[t_i * 7 + 0]);
                  X[i * 7 + 1] = X[t_i * 7 + 1] + dt * (cd4_prg[0] * X[t_i * 7 + 0] - 
                  (cd4_prg[1] + mort_rate_t[i * 7 + 1] + test_rate_t[i * 7 + 1]) * X[t_i * 7 + 1]);
                  X[i * 7 + 2] = X[t_i * 7 + 2] + dt * (cd4_prg[1] * X[t_i * 7 + 1] - 
                  (cd4_prg[2] + mort_rate_t[i * 7 + 2] + test_rate_t[i * 7 + 2]) * X[t_i * 7 + 2]);
                  X[i * 7 + 3] = X[t_i * 7 + 3] + dt * (cd4_prg[2] * X[t_i * 7 + 2] - 
                  (cd4_prg[3] + mort_rate_t[i * 7 + 3] + test_rate_t[i * 7 + 3]) * X[t_i * 7 + 3]);
                  X[i * 7 + 4] = X[t_i * 7 + 4] + dt * (cd4_prg[3] * X[t_i * 7 + 3] - 
                  (cd4_prg[4] + mort_rate_t[i * 7 + 4] + test_rate_t[i * 7 + 4]) * X[t_i * 7 + 4]);
                  X[i * 7 + 5] = X[t_i * 7 + 5] + dt * (cd4_prg[4] * X[t_i * 7 + 4] - 
                  (cd4_prg[5] + mort_rate_t[i * 7 + 5] + test_rate_t[i * 7 + 5]) * X[t_i * 7 + 5]);
                  X[i * 7 + 6] = X[t_i * 7 + 6] + dt * (cd4_prg[5] * X[t_i * 7 + 5] - 
                  (mort_rate_t[i * 7 + 6] + test_rate_t[i * 7 + 6]) * X[t_i * 7 + 6]);
                  
                  NumericVector nb_dx_z(7);
                  IntegerVector ind_t0 = seq(t_i * 7 + 0, t_i * 7 + 6);
                  for (int z = 0; z < 7; ++z) {
                  int ind_tz = ind_t0[z];
                  nb_dx_z[z] = (X[ind_tz] * test_rate_t[ind_tz + 7]) * dt;
                  }
                  nb_dx_tot[i] = sum(nb_dx_z);
                  }
                  
                  // outcomes are stored here for each year
                  NumericVector prb1yr_n = nb_dx_tot[seq(0, tind1yr - 1)];
                  prb1yr[n] = sum(prb1yr_n);
                  }
                  
                  DataFrame val = DataFrame::create(
                  Named("year") = year,
                  Named("age") = age,
                  Named("sex") = sex,
                  Named("prb1yr") = prb1yr);
                  
                  return(val);
                  }')
