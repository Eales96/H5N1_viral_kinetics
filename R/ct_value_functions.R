
pop_distributions <- function(ft,thresh, n_samples=1000, time=seq(0,30,1)){
  
  t_peak_mn <- median(ft$t_peak_mn)
  t_peak_sd <- mean(ft$t_peak_sd)
  
  peak_mn <- median(ft$peak_mn)
  peak_sd <- median(ft$peak_sd)
  
  rise_mn <- median(ft$rise_mn)
  rise_sd <- median(ft$rise_sd)
  
  decay_mn <- median(ft$decay_mn)
  decay_sd <- median(ft$decay_sd)
  
  t_peak <- rnorm(n_samples, t_peak_mn, t_peak_sd)
  peak <- rnorm(n_samples, peak_mn, peak_sd)
  rise <- exp(rnorm(n_samples, log(rise_mn), rise_sd))
  decay <- exp(rnorm(n_samples, log(decay_mn), decay_sd))
  
  
  df_mod <- data.frame(time = time,
                       y = NA,
                       y_lb95 = NA,
                       y_ub95 = NA,
                       y_lb50 = NA,
                       y_ub50 = NA,
                       prop=NA)
  for(i in 1:nrow(df_mod)){
    post <- c()
    for(j in 1:length(t_peak)){
      if(df_mod$time[i]<t_peak[j]){
        post <- c(post, peak[j] - (df_mod$time[i]-t_peak[j])*rise[j])
      } else{
        post <- c(post, peak[j] + (df_mod$time[i]-t_peak[j])*decay[j])
        
      }
      
    }
    quan <- quantile(post, c(0.5,0.025,0.25,0.75,0.975))
    df_mod$y[i] <- quan[[1]]
    df_mod$y_lb95[i] <- quan[[2]]
    df_mod$y_lb50[i] <- quan[[3]]
    df_mod$y_ub50[i] <- quan[[4]]
    df_mod$y_ub95[i] <- quan[[5]]
    df_mod$prop[i] <- sum(post<thresh)/length(post)
    
  }
  return(df_mod)
  
  
}



t_inf <- function(params, thresh, n_samples=1000){
  
  t_peak_mn <- params[1]
  t_peak_sd <- params[2]
  
  peak_mn <- params[3]
  peak_sd <- params[4]
  
  rise_mn <- params[5]
  rise_sd <- params[6]
  
  decay_mn <- params[7]
  decay_sd <- params[8]
  
  t_peak <- rnorm(n_samples, t_peak_mn, t_peak_sd)
  peak <- rnorm(n_samples, peak_mn, peak_sd)
  rise <- exp(rnorm(n_samples, log(rise_mn), rise_sd))
  decay <- exp(rnorm(n_samples, log(decay_mn), decay_sd))
  
  post <- c()
  for(i in 1:n_samples){
    
    if(peak[i] <thresh){
      time1 = t_peak[i] + ( (peak[i]-thresh)/rise[i])
      time2 = t_peak[i] + ( (thresh - peak[i])/decay[i])
      post <- c(post, time2-time1)
      
    } else{
      post <- c(post, 0)
    }
    
  }
  
  return(mean(post))
  
}

get_posterior_mean_tinf <- function(ft,thresh, n_samples=1000){
  
  post <- c()
  for(i in 1:nrow(ft$peak_mn)){
    print(i)
    params=c(ft$t_peak_mn[i],
             ft$t_peak_sd[i],
             ft$peak_mn[i],
             ft$peak_sd[i],
             ft$rise_mn[i],
             ft$rise_sd[i],
             ft$decay_mn[i],
             ft$decay_sd[i])
    if(length(thresh)==1){
      post <- c(post, t_inf(params, thresh, n_samples) )
    } else{
      post <- c(post, t_inf(params, thresh[i], n_samples) )
    }
    
  }
  
  return(post)
}


get_dist_peakCt <- function(ft, n_samples=1000){
  
  post <- c()
  for(i in 1:nrow(ft$peak_mn)){
    print(i)
    post <- c(post, rnorm(n_samples, ft$peak_mn[i],ft$peak_sd[i]))
    
    
  }
  
  return(post)
}
################################################################################

get_posteriorH_all <- function(ft, time = seq(0,30,0.1)){
  df_mod <- data.frame(time = time)
  for(i in 1:nrow(df_mod)){
    post <- c()
    for(j in 1:nrow(ft$t_peak_mn)){
      if(df_mod$time[i]<ft$t_peak_mn[j]){
        post <- c(post, ft$peak_mn[j] - (df_mod$time[i]-ft$t_peak_mn[j])*ft$rise_mn[j])
      } else{
        post <- c(post, ft$peak_mn[j] + (df_mod$time[i]-ft$t_peak_mn[j])*ft$decay_mn[j])
        
      }
      
    }
    quan <- quantile(post, c(0.5,0.025,0.25,0.75,0.975))
    df_mod$y[i] <- quan[[1]]
    df_mod$y_lb95[i] <- quan[[2]]
    df_mod$y_lb50[i] <- quan[[3]]
    df_mod$y_ub50[i] <- quan[[4]]
    df_mod$y_ub95[i] <- quan[[5]]
    
  }
  return(df_mod)
}





get_posteriorH_all_i <- function(ft, index, time = seq(0,30,0.1)){
  df_mod <- data.frame(time = time, index=index)
  for(i in 1:nrow(df_mod)){
    post <- c()
    for(j in 1:length(ft$t_peak_i[,index])){
      if(df_mod$time[i]<ft$t_peak_i[,index][j]){
        post <- c(post, ft$peak_i[,index][j] - (df_mod$time[i]-ft$t_peak_i[,index][j])*ft$rise_i[,index][j])
      } else{
        post <- c(post, ft$peak_i[,index][j] + (df_mod$time[i]-ft$t_peak_i[,index][j])*ft$decay_i[,index][j])
        
      }
      
    }
    quan <- quantile(post, c(0.5,0.025,0.25,0.75,0.975))
    df_mod$y[i] <- quan[[1]]
    df_mod$y_lb95[i] <- quan[[2]]
    df_mod$y_lb50[i] <- quan[[3]]
    df_mod$y_ub50[i] <- quan[[4]]
    df_mod$y_ub95[i] <- quan[[5]]
    
  }
  return(df_mod)
}




posterior_decline <- function(ft, time = seq(0,30,0.1)){
  
  df_mod <- data.frame(time = time)
  for(i in 1:nrow(df_mod)){
    post <- c()
    for(j in 1:nrow(ft$t_peak_mn)){
      post <- c(post, 0 + df_mod$time[i]*ft$decay_mn[j])
      
    }
    quan <- quantile(post, c(0.5,0.025,0.25,0.75,0.975))
    df_mod$y[i] <- quan[[1]]
    df_mod$y_lb95[i] <- quan[[2]]
    df_mod$y_lb50[i] <- quan[[3]]
    df_mod$y_ub50[i] <- quan[[4]]
    df_mod$y_ub95[i] <- quan[[5]]
    
  }
  return(df_mod)
  
  
  
}


posterior_decline_i <- function(ft,index, time = seq(0,30,0.1)){
  
  df_mod <- data.frame(time = time)
  df_mod$index <- index
  for(i in 1:nrow(df_mod)){
    post <- c()
    for(j in 1:nrow(ft$t_peak_mn)){
      post <- c(post, ft$start_i[j,index] + df_mod$time[i]*ft$decayED_i[j, index])
      
    }
    quan <- quantile(post, c(0.5,0.025,0.25,0.75,0.975))
    df_mod$y[i] <- quan[[1]]
    df_mod$y_lb95[i] <- quan[[2]]
    df_mod$y_lb50[i] <- quan[[3]]
    df_mod$y_ub50[i] <- quan[[4]]
    df_mod$y_ub95[i] <- quan[[5]]
    
    
  }
  return(df_mod)
  
  
  
}



inf_virus_posterior <- function(ft, Ct=seq(10,30,0.1)){
  
  df_mod <- data.frame(Ct= Ct)
  for(i in 1:nrow(df_mod)){
    
    post <- ft$beta0 / (1 + exp(ft$beta1*(ft$beta2-df_mod$Ct[i]) ) )
    quan <- quantile(post, c(0.5,0.025,0.25,0.75,0.975))
    df_mod$y[i] <- quan[[1]]
    df_mod$y_lb95[i] <- quan[[2]]
    df_mod$y_lb50[i] <- quan[[3]]
    df_mod$y_ub50[i] <- quan[[4]]
    df_mod$y_ub95[i] <- quan[[5]]
    
  }
  
  return(df_mod)
  
}








pop_distributions_all <- function(ft,thresh, n_samples=1000, time=seq(0,30,1)){
  
  mat_mod <- matrix(data=NA, nrow = length(ft$t_peak_mn)*length(time), ncol = 7)
  
  for(i in 1:length(ft$t_peak_mn)){
    
    if(i%%100==0){
      print(i)
    }
    
    t_peak <- rnorm(n_samples, ft$t_peak_mn[i], ft$t_peak_sd[i])
    peak <- rnorm(n_samples, ft$peak_mn[i], ft$peak_sd[i])
    rise <- exp(rnorm(n_samples, log(ft$rise_mn[i]), ft$rise_sd[i]))
    decay <- exp(rnorm(n_samples, log(ft$decay_mn[i]), ft$decay_sd[i]))
    
    for(j in 1:length(time)){
      
      vec1 <- peak - (time[j]-t_peak)*rise
      vec2 <- peak + (time[j]-t_peak)*decay
      index<- t_peak<j
      vec<- c(vec1[index==FALSE], vec2[index==TRUE])
      
      quan <- quantile(vec, c(0.025,0.25, 0.5, 0.75, 0.975))
      
      
      mat_mod[(i-1)*length(time)+j,] <- c(i,time[j],quan)
      
    }
    
    
  }
  
  return(mat_mod)
  
  
}


t_inf_distribution <- function(ft,thresh, n_samples=1000, time=seq(0,30,1)){
  
  mat_mod <- matrix(data=NA, nrow = length(ft$t_peak_mn)*length(time), ncol = 7)
  
  
  t_peak <- rnorm(n_samples, median(ft$t_peak_mn), median(ft$t_peak_sd) )
  peak <- rnorm(n_samples, median(ft$peak_mn), median(ft$peak_sd) )
  rise <- exp(rnorm(n_samples, log(median(ft$rise_mn)), median(ft$rise_sd) ))
  decay <- exp(rnorm(n_samples, log(median(ft$decay_mn) ), median(ft$decay_sd) ))
  
  time1 = t_peak + ( (peak-thresh)/rise)
  time2 = t_peak + ( (thresh - peak)/decay)
  
  index <- peak<thresh
  
  row1 <- data.frame(time1=time1[index],
                     time2=time2[index],
                     t_inf = time2[index]-time1[index])
  
  row2 <- data.frame(time1=rep(NA, sum(index==FALSE)),
                     time2=rep(NA, sum(index==FALSE)),
                     t_inf = rep(0, sum(index==FALSE)))
  
  
  
  return(rbind(row1, row2))
  
  
}


t_inf_distribution_i <- function(ft,thresh, n_samples=1000, time=seq(0,30,1), i){
  
  #mat_mod <- matrix(data=NA, nrow = length(ft$t_peak_mn)*length(time), ncol = 7)
  
  t_peak <- rnorm(n_samples, ft$t_peak_mn[i], ft$t_peak_sd[i] )
  peak <- rnorm(n_samples, ft$peak_mn[i], ft$peak_sd[i] )
  rise <- exp(rnorm(n_samples, log(ft$rise_mn[i]), ft$rise_sd[i] ))
  decay <- exp(rnorm(n_samples, log(ft$decay_mn[i] ), ft$decay_sd[i] ))
  
  time1 = t_peak + ( (peak-thresh)/rise)
  time2 = t_peak + ( (thresh - peak)/decay)
  
  index <- peak<thresh
  
  row1 <- data.frame(time1=time1[index],
                     time2=time2[index],
                     t_inf = time2[index]-time1[index])
  
  row2 <- data.frame(time1=rep(NA, sum(index==FALSE)),
                     time2=rep(NA, sum(index==FALSE)),
                     t_inf = rep(0, sum(index==FALSE)))
  
  
  
  return(rbind(row1, row2))
  
  
}



prop_inf_prev <- function(tinf1, t_i, p_i){
  
  tinf1$tp1 <- t_i-tinf1$time1
  tinf1$tp2 <- t_i-tinf1$time2
  in1 <- tinf1$tp1 < 0 
  in2 <- tinf1$tp1 >= 0 & tinf1$tp2<0
  in3 <- tinf1$tp2 >= 0
  
  tinf1$prev <- NA
  if(sum(in1)>0){
    tinf1[in1,]$prev <- tinf1[in1,]$t_inf * p_i
  }
  if(sum(in2)>0){
    tinf1[in2,]$prev <- -tinf1[in2,]$tp2 * p_i
  }
  if(sum(in3)>0){
    tinf1[in3,]$prev <- 0
  }
  
  return(sum(tinf1$prev))
  
  
  
}


get_posteriors <- function(post, params, model_name){
  df <- data.frame()
  for(i in 1:length(params)){
    vals <- post[names(post)==params[i]][[1]]
    quan <- quantile(vals, c(0.025,0.25,0.5,0.75,0.975))
    row <- data.frame(model_name = model_name,
                      param_name = params[i],
                      y= quan[[3]],
                      y_lb95 = quan[[1]],
                      y_lb50 = quan[[2]],
                      y_ub95 = quan[[5]],
                      y_ub50 = quan[[4]])
    
    df <- rbind(df, row)
  }
  return(df)
}


prop_inf_prev_eff <- function(tinf1, t_i, p_i){
  tinf1 <- tinf1[is.na(tinf1$time1)==FALSE,]
  sm <- c()
  for(i in 1:length(t_i)){
    sm <- c(sm, sum(tinf1[t_i[i]-tinf1$time1<0,]$t_inf) + sum( (tinf1[t_i[i]-tinf1$time1>=0 & t_i[i]-tinf1$time2 < 0,]$time2-t_i[i]) ))
  }
  sm <- as.matrix(sm)
  return(sm%*%p_i)
  
  
  
}
