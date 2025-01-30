// Copyright 2024 Oliver Eales
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

data {
  int num_data_x_y;                          // number of data points
  real<lower=0> logTCID50_x_y[num_data_x_y];
  real<lower=0, upper=40> Ct_x_y[num_data_x_y];
  
  int num_data_x_ny;                         // number of data points
  real<lower=0, upper=40> Ct_x_ny[num_data_x_ny];
}

parameters {
  vector<lower=0>[num_data_x_y] mod_Ct_x_y;
  vector<lower=0>[num_data_x_ny] mod_Ct_x_ny;
  
  //real<lower=0, upper=1.05> logTCID50_x_ny[num_data_x_ny];
  real<upper=1.05> logTCID50_x_ny[num_data_x_ny];
  
  real<lower=0> beta0;
  real<lower=-100, upper=0> beta1;
  real<lower=0, upper=30> beta2;
  
  real<lower=0> sigma_logTCID50;
  real<lower=0> sigma_Ct;
}


transformed parameters {
  vector<lower=0>[num_data_x_y] mod_logTCID50_x_y;
  vector<lower=0>[num_data_x_ny] mod_logTCID50_x_ny;
  
  mod_logTCID50_x_y = beta0 / (1 + exp(beta1*(beta2-mod_Ct_x_y) ) );
  mod_logTCID50_x_ny = beta0 / (1 + exp(beta1*(beta2-mod_Ct_x_ny) ) );
  
  

}

model {
  logTCID50_x_y ~ normal(mod_logTCID50_x_y, sigma_logTCID50);
  logTCID50_x_ny ~ normal(mod_logTCID50_x_ny, sigma_logTCID50);
  
  Ct_x_y ~ normal(mod_Ct_x_y, sigma_Ct);
  Ct_x_ny ~ normal(mod_Ct_x_ny, sigma_Ct);

}
