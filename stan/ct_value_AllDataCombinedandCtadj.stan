
data {
  int num_data1;                          // number of data points
  real<lower=0> Ct_value1[num_data1];
  real<lower=0> time1[num_data1];
  int<lower=0> cow_number1[num_data1];
  
  int num_data2;                          // number of data points
  real<lower=0> Ct_value2[num_data2];
  real<lower=0> time2[num_data2];
  int<lower=0> cow_number2[num_data2];
  
  int num_data4;                          // number of data points
  real<lower=0> Ct_value4[num_data4];
  real<lower=0> time4[num_data4];
  int<lower=0> cow_number4[num_data4];
  
  int num_dataC1;                          // number of data points
  real<lower=0> timeC1[num_dataC1];
  int<lower=0> cow_numberC1[num_dataC1];
  
  int num_dataC2;                          // number of data points
  real<lower=0> timeC2[num_dataC2];
  int<lower=0> cow_numberC2[num_dataC2];
  
  int num_dataC4;                          // number of data points
  real<lower=0> timeC4[num_dataC4];
  int<lower=0> cow_numberC4[num_dataC4];
  
  int num_cows;
  
  
  // Extra data from other paper
  int num_dataED;
  real<lower=0> Ct_valueED[num_dataED];
  real<lower=0> timeED[num_dataED];
  int<lower=0> cow_numberED[num_dataED];
  
  int num_dataEDC;
  real<lower=0> timeEDC[num_dataEDC];
  int<lower=0> cow_numberEDC[num_dataEDC];
  
  int num_cowsED;
  
}

parameters {
  
  // censored data
  real<lower=40> Ct_valueC1[num_dataC1];
  real<lower=38> Ct_valueC2[num_dataC2];
  real<lower=38> Ct_valueC4[num_dataC4];
  
  real<lower=45> Ct_valueEDC[num_dataEDC];
  
  // model parameters
  real<lower=0, upper=10> t_peak_i[num_cows];
  real<lower=0> peak_i[num_cows];
  vector<lower=0.01, upper=100>[num_cows] rise_i;
  vector<lower=0.01, upper=100>[num_cows] decay_i;
  
  
  real<lower=0, upper=10> t_peak_mn;
  real<lower=0> t_peak_sd;
  
  real<lower=0.01, upper=100> rise_mn;
  real<lower=0> rise_sd;
  
  real<lower=0> peak_mn;
  real<lower=0> peak_sd;
  real<lower=0.01, upper=100> decay_mn;
  real<lower=0> decay_sd;
  
  // Single decay parameters
  real<lower=0> start_i[num_cowsED];
  vector<lower=0.01, upper=100>[num_cowsED] decayED_i;
  
  // CT adjustment parameter 
  real<lower=-10, upper=10> ct_adj;
  real<lower=-10, upper=10> ct_adj4;
  
  // noise parameters
  real<lower=0> sigma_Ct;
}


transformed parameters{
  
  vector<lower=0>[num_dataED] mod_Ct_ED;
  for(i in 1:num_dataED){
    mod_Ct_ED[i] = start_i[cow_numberED[i]] + timeED[i]*decayED_i[cow_numberED[i]];
  }
  
  vector<lower=0>[num_dataEDC] mod_Ct_EDC;
  for(i in 1:num_dataEDC){
    mod_Ct_EDC[i] = start_i[cow_numberEDC[i]] + timeEDC[i]*decayED_i[cow_numberEDC[i]];
  }
  
  vector<lower=0>[num_data1] mod_Ct1;
  
  for(i in 1:num_data1){
    if(time1[i]<t_peak_i[cow_number1[i]] ){
      mod_Ct1[i] = peak_i[cow_number1[i]] - (time1[i]-t_peak_i[cow_number1[i]])*rise_i[cow_number1[i]];
    } else{
      mod_Ct1[i] = peak_i[cow_number1[i]] + (time1[i]-t_peak_i[cow_number1[i]])*decay_i[cow_number1[i]];
    }
  }
  
  vector<lower=0>[num_dataC1] mod_CtC1;
  
  for( i in 1:num_dataC1){
    if(timeC1[i]<t_peak_i[cow_numberC1[i]] ){
      mod_CtC1[i] = peak_i[cow_numberC1[i]] - (timeC1[i]-t_peak_i[cow_numberC1[i]])*rise_i[cow_numberC1[i]];
    } else{
      mod_CtC1[i] = peak_i[cow_numberC1[i]] + (timeC1[i]-t_peak_i[cow_numberC1[i]])*decay_i[cow_numberC1[i]];
    }
  }
  
  vector<lower=0>[num_data2] mod_Ct2;
  
  for(i in 1:num_data2){
    if(time2[i]<t_peak_i[cow_number2[i]] ){
      mod_Ct2[i] = ct_adj + peak_i[cow_number2[i]] - (time2[i]-t_peak_i[cow_number2[i]])*rise_i[cow_number2[i]];
    } else{
      mod_Ct2[i] = ct_adj + peak_i[cow_number2[i]] + (time2[i]-t_peak_i[cow_number2[i]])*decay_i[cow_number2[i]];
    }
  }
  
  vector<lower=0>[num_dataC2] mod_CtC2;
  
  for( i in 1:num_dataC2){
    if(timeC2[i]<t_peak_i[cow_numberC2[i]] ){
      mod_CtC2[i] = ct_adj + peak_i[cow_numberC2[i]] - (timeC2[i]-t_peak_i[cow_numberC2[i]])*rise_i[cow_numberC2[i]];
    } else{
      mod_CtC2[i] = ct_adj + peak_i[cow_numberC2[i]] + (timeC2[i]-t_peak_i[cow_numberC2[i]])*decay_i[cow_numberC2[i]];
    }
  }
  
  vector<lower=0>[num_data4] mod_Ct4;
  
  for(i in 1:num_data4){
    if(time4[i]<t_peak_i[cow_number4[i]] ){
      mod_Ct4[i] = ct_adj4 + peak_i[cow_number4[i]] - (time4[i]-t_peak_i[cow_number4[i]])*rise_i[cow_number4[i]];
    } else{
      mod_Ct4[i] = ct_adj4 + peak_i[cow_number4[i]] + (time4[i]-t_peak_i[cow_number4[i]])*decay_i[cow_number4[i]];
    }
  }
  
  vector<lower=0>[num_dataC4] mod_CtC4;
  
  for( i in 1:num_dataC4){
    if(timeC4[i]<t_peak_i[cow_numberC4[i]] ){
      mod_CtC4[i] = ct_adj4 + peak_i[cow_numberC4[i]] - (timeC4[i]-t_peak_i[cow_numberC4[i]])*rise_i[cow_numberC4[i]];
    } else{
      mod_CtC4[i] = ct_adj4 + peak_i[cow_numberC4[i]] + (timeC4[i]-t_peak_i[cow_numberC4[i]])*decay_i[cow_numberC4[i]];
    }
  }
}

model {
  
  Ct_value1 ~ normal(  mod_Ct1   , sigma_Ct);
  Ct_value2 ~ normal(  mod_Ct2   , sigma_Ct);
  Ct_value4 ~ normal(  mod_Ct4   , sigma_Ct);
  Ct_valueC1 ~ normal(  mod_CtC1   , sigma_Ct);
  Ct_valueC2 ~ normal(  mod_CtC2   , sigma_Ct);
  Ct_valueC4 ~ normal(  mod_CtC4   , sigma_Ct);
  
  Ct_valueED ~ normal(  mod_Ct_ED   , sigma_Ct);
  Ct_valueEDC ~ normal(  mod_Ct_EDC   , sigma_Ct);
  
  t_peak_i ~ normal(t_peak_mn, t_peak_sd);
  peak_i ~ normal(peak_mn, peak_sd);
  log(rise_i) ~ normal(log(rise_mn), rise_sd);
  log(decay_i) ~ normal(log(decay_mn), decay_sd);
  log(decayED_i) ~ normal(log(decay_mn), decay_sd);
  
}
