
data {
  int num_data;                          // number of data points
  real<lower=0> Ct_value[num_data];
  real<lower=0> time[num_data];
  int<lower=0> cow_number[num_data];
  
  int num_dataC1;                          // number of data points
  real<lower=0> timeC1[num_dataC1];
  int<lower=0> cow_numberC1[num_dataC1];
  
  int num_dataC2;                          // number of data points
  real<lower=0> timeC2[num_dataC2];
  int<lower=0> cow_numberC2[num_dataC2];
  
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
  
  vector<lower=0>[num_data] mod_Ct;
  
  for(i in 1:num_data){
    if(time[i]<t_peak_i[cow_number[i]] ){
      mod_Ct[i] = peak_i[cow_number[i]] - (time[i]-t_peak_i[cow_number[i]])*rise_i[cow_number[i]];
    } else{
      mod_Ct[i] = peak_i[cow_number[i]] + (time[i]-t_peak_i[cow_number[i]])*decay_i[cow_number[i]];
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
  
  vector<lower=0>[num_dataC2] mod_CtC2;
  
  for( i in 1:num_dataC2){
    if(timeC2[i]<t_peak_i[cow_numberC2[i]] ){
      mod_CtC2[i] = peak_i[cow_numberC2[i]] - (timeC2[i]-t_peak_i[cow_numberC2[i]])*rise_i[cow_numberC2[i]];
    } else{
      mod_CtC2[i] = peak_i[cow_numberC2[i]] + (timeC2[i]-t_peak_i[cow_numberC2[i]])*decay_i[cow_numberC2[i]];
    }
  }
}

model {
  
  Ct_value ~ normal(  mod_Ct   , sigma_Ct);
  Ct_valueC1 ~ normal(  mod_CtC1   , sigma_Ct);
  Ct_valueC2 ~ normal(  mod_CtC2   , sigma_Ct);
  
  Ct_valueED ~ normal(  mod_Ct_ED   , sigma_Ct);
  Ct_valueEDC ~ normal(  mod_Ct_EDC   , sigma_Ct);
  
  t_peak_i ~ normal(t_peak_mn, t_peak_sd);
  peak_i ~ normal(peak_mn, peak_sd);
  log(rise_i) ~ normal(log(rise_mn), rise_sd);
  log(decay_i) ~ normal(log(decay_mn), decay_sd);
  log(decayED_i) ~ normal(log(decay_mn), decay_sd);

}
