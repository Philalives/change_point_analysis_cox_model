# Author: Philani Brian Mpofu
# Purpose: Develop code for identifying change-points for piecewise Cox Models

library(tidyverse)
library(dplyr)
library(survival)
set.seed(2985)


#-- Generate piecewise exponential data

# 1. Generate uniform random variable
U <- runif(1000)

generate_piece_exp <- function(u){
  n = length(u)
  id = 1:n
  # 1. Treatment indicator
  Z <- rbinom(n=n, size=1, prob=0.5)
  # hazards of events
  #-- before time 2
  h1 <- exp(-0.25 + 0.5*Z)
  #-- after time 2
  h2 <- exp(-0.25 - 0.5*Z)
  
  cut = exp(-2*h1)
  
  event_time = if_else(u > cut, -log(u)/h1, 2 - (log(u) + 2*h1)/h2)
  
  cens_time = rexp(n=n, rate = 0.1)
  
  t = pmin(event_time, cens_time)
  
  event = as.numeric(event_time == t)
  
  out = data.frame(id =id, Z = Z,
                   t = t,
                   event = event)
  
  return(out)
}

data_gen <- generate_piece_exp(U)

data_gen %>% 
  group_by(Z,event) %>% 
  summarise(n= n()) %>% 
  mutate(freq = round(100*n/sum(n), 2))

# Cox model
cox.mod <- coxph(Surv(t, event = event)~ Z, data= data_gen)
summary(cox.mod)

#-- 1. Scaled-Schoenfeld residual approach
mod_sch <-  cox.zph(cox.mod)
mod_sch
plot(mod_sch)

#-- Find all the event times
KM_times <-  survfit(Surv(t, event)~ Z, 
                     data= data_gen)
k <- summary(KM_times)
event_times <- KM_times$time
event_times <- event_times[c(-1, -length(event_times ))]

plot(k$time[k$strata=="Z=1"], 
     k$cumhaz[k$strata=="Z=1"], type="l")

#-- Let's try a piecewise linear model;((
find_cut = function(time_vec){
  time_vec_list = as.list(time_vec)
  extract_part_lik <-  function(t_thresh){
    # create counting process style dataset
    cp_data <- survSplit(Surv(t, event) ~.,
                         data= data_gen,
                         zero= 0, 
                         cut= c(t_thresh),
                         episode = "interval",
                         event="event")
    
    cox_mod_piece <- coxph(Surv(tstart, t, event) ~ Z:strata(interval) ,
                           cluster =id,
                           data = cp_data,
                           ties ="breslow")
    loglike = as.numeric(logLik(cox_mod_piece))
    return(loglike)
  }
  
  list_log_lik = lapply(time_vec_list, function(x) extract_part_lik(x))
  loglike = do.call(rbind,  list_log_lik)
  out = data.frame(time_vec, loglike)
  return(out)
}

times <- seq(0.5, 7, length.out = 1000)
partial_like_time = find_cut(event_times)

plot(partial_like_time$time_vec, partial_like_time$loglike, type="l")


cp_data <- survSplit(Surv(t, event) ~.,
                     data= data_gen,
                     zero= 0, 
                     cut= c(1),
                     episode = "interval",
                     event="event")

cox_mod_piece <- coxph(Surv(tstart, t, event) ~ Z:strata(interval),
                       cluster = id,
                       data = cp_data,
                       ties ="breslow")
summary(cox_mod_piece)
mod_sch1 <-  cox.zph(cox_mod_piece)
mod_sch1
plot(mod_sch1)
