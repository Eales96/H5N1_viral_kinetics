
data {
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
  real<lower=45> Ct_valueEDC[num_dataEDC];
  
  //Model parameters
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
}

model {
  
  Ct_valueED ~ normal(  mod_Ct_ED   , sigma_Ct);
  Ct_valueEDC ~ normal(  mod_Ct_EDC   , sigma_Ct);
  
  log(decayED_i) ~ normal(log(decay_mn), decay_sd);

}
