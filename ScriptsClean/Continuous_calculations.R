########################################################################################
# Author: Houston F. Lester                                                                
# Description: This file is intended to perform the calculations necessary to obtain
# the desired R squared values for the data generation assessing the informative hypothesis
# testing approach. Some of this stuff is redundant with other code, but is repeated 
# to make each script reproducible. The end result of this code is to create the Design
# matrix that is included in the Data folder. 
#####################################################################################

#Setup ----
library(dplyr)


# Mean and covariance matrices -----

mean_vec3 <- matrix(data = c(0, 0, 0),  nrow = 3, ncol = 1)
cov_mat_int3_cont <- matrix(data = c(1, .3, 0, 
                                     .3, 1, 0, 0,
                                     0, 1.09), nrow = 3, ncol = 3)

## Corresponds to the 2 interaction terms conditions (i.e., 2 two way interaction model)
mean_vec5 <- matrix(data = c(0, 0, 0, 0, 0),  nrow = 5, ncol = 1)
cov_mat_int5_cont <- matrix(data = c(1, .3, .3, 0, 0, 
                                     .3, 1, .3, 0, 0, 
                                     .3, .3, 1, 0, 0,
                                     0, 0, 0, 1.09, .39,
                                     0, 0, 0, .39, 1.09), nrow = 5, ncol = 5) 


## Corresponds to the 4 interaction terms conditions (i.e., 3 way interaction model)
mean_vec7 <- matrix(data = c(0, 0, 0, 0, 0, 0, 0),  nrow = 7, ncol = 1)
cov_mat_int7_cont <- matrix(data = c(1, .3, .3, 0, 0, 0, .48,  
                                     .3, 1, .3, 0, 0, 0, .48, 
                                     .3, .3, 1, 0, 0, 0, .48,
                                     0, 0, 0, 1.09, .39, .39, 0,
                                     0, 0, 0, .39, 1.09, .39, 0,
                                     0, 0, 0, .39, .39, 1.09, 0,
                                     .48, .48, .48, 0, 0, 0, 1.756), nrow = 7, ncol = 7) 


mu_list <- list(vec3 = mean_vec3, vec5 = mean_vec5, mean_vec7 = mean_vec7)
cov_mat_list <- list(single_mod = cov_mat_int3_cont, two_mod = cov_mat_int5_cont, four_mod = cov_mat_int7_cont)


# Functions to calculate the quantities of interest -----
# Parameter explanations - fsqu is a measure of the interaction effect size
#                       - Rsqu_red is the amount of variability explained by the reduced model
#                       - npred_red is the number of predictors in the reduced model
#                       - npred_full is the number of predictors in the full model
#                       - full_cov_mat is the full model covariance matrix obtained from the Expected_value_calculation.R file



## Clunky but effective way to calculate the regression coefficients needed for the one and 2 interaction term conditions. 
Rsq_f_gamma2 <- function(fsqu, Rsqu_red = .104, npred_red, npred_full, full_cov_mat){
  Rsq_needed <- (fsqu + Rsqu_red)/(1 + fsqu)
  diffz <- Rsq_needed - Rsqu_red
  test <- data.frame(full_cov_mat, rnumz = 1:nrow(full_cov_mat))
  cov_sum <- test %>%
    mutate(cov_contrib = if_else(rnumz > npred_red, rowSums(test) - rnumz, 0)) %>%
    summarise(new = sum(cov_contrib))
  new_stuff <- cov_sum
  beta <- sqrt(diffz/new_stuff)
  toret <- list(desired_rsq = Rsq_needed, beta = beta)
  return(toret)
}

##Quadratic Equation used to calculate desired regression coefficients in four interaction term conditions using Rsq_f_gamma4. 
quadratic <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
    result = c(x_1,x_2)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
  }
  else {"There are no real roots."} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
  b^2-4*a*c
}

Rsq_f_gamma4 <- function(fsqu, Rsqu_red = .104, npred_red, npred_full, full_cov_mat, maineff){
  Rsq_needed <- (fsqu + Rsqu_red)/(1 + fsqu)
  diffz <- Rsq_needed - Rsqu_red
  if(diffz == 0){ 
    beta <- 0
  }
  else {
    new_var <- sum(diag(full_cov_mat)[(npred_red + 1):npred_full])
    new_cov <- 6*full_cov_mat[5,4]
    beta <- quadratic(a = new_var + new_cov, b = 6*maineff*full_cov_mat[3,7], c = - diffz)
  }
  toret <- list(desired_rsq = Rsq_needed, beta = beta)
  return(toret)
}

#Creating Design Matrix -----
## Manipulated factors ----
N <- c(500, 1000, 2000)
f_sq <- c(0, .002, .01, .02)
n_inter <- c(1, 2, 4)
continuous <- c(0, 1) 

three_pred_main <- 0.147196

stat_design <- expand.grid(N, f_sq, n_inter, continuous)
names(stat_design) <- c("N", "f_sq", "n_inter", "continuous")

cont <- as.data.frame(stat_design)
cont <- filter(cont, continuous == 1)

cont51 <- cont %>% 
  mutate(npred_red = if_else(n_inter == 1, 2, 3),
         npred_full = if_else(n_inter == 1, 3, if_else(n_inter == 2, 5, 7)),
         beta1 = if_else(npred_red == 2, .2, three_pred_main), beta2 = if_else(npred_red == 2, .2, three_pred_main), 
         beta3 = if_else(npred_red == 3, three_pred_main, 0))

continuous_Rsquared_betas <- vector(mode = "list", length = nrow(cont51))

for(i in 1:nrow(cont51)) {
  if(i < 13) {
    continuous_Rsquared_betas[[i]] <- Rsq_f_gamma2(fsqu = cont51$f_sq[i], 
                                                  npred_red = cont51$npred_red[i], 
                                                  npred_full = cont51$npred_full[i], 
                                                  full_cov_mat = cov_mat_int3_cont)
  }
  else if((i >=  13) & (i < 25) ){ # second case D=0
    continuous_Rsquared_betas[[i]] <- Rsq_f_gamma2(fsqu = cont51$f_sq[i], npred_red = cont51$npred_red[i], 
                                                    npred_full = cont51$npred_full[i], full_cov_mat = cov_mat_int5_cont)
  }
  else {
    continuous_Rsquared_betas[[i]] <- Rsq_f_gamma4(fsqu = cont51$f_sq[i], 
                                                    npred_red = cont51$npred_red[i], 
                                                    npred_full = cont51$npred_full[i], 
                                                    full_cov_mat = cov_mat_int7_cont, maineff =  cont51$beta1[i])
    }
  }




better_list <- lapply(continuous_Rsquared_betas, function(x) {
      temp1 <- unlist(x)
      r_sq <- temp1[1]
      new_beta <- temp1[2]
      new_dat <- data.frame(r_sq = r_sq, new_beta = new_beta)
  }
)

R_squared_betas <- do.call(rbind, better_list)
rownames(R_squared_betas) <- NULL

cont52 <- cbind(cont51, R_squared_betas)

cont53 <- cont52 %>%
  mutate(beta1 = if_else(npred_red == 2, .2, three_pred_main), beta2 = if_else(npred_red == 2, .2, three_pred_main), 
         beta3 = if_else(npred_red == 3, three_pred_main, 0), int1 = new_beta, int2 = if_else((npred_full - npred_red) >= 2, new_beta, 0),
         int3 = if_else((npred_full - npred_red) == 4, new_beta, 0 ), three_way = if_else((npred_full - npred_red) == 4, new_beta, 0)) 

Design_redo <- cont53 %>%
  select(1:6, 10, 11, 7:9, 12:15)

round(Design_redo, 4) == round(Design, 4)



round(Design_redo, 4) == round(Design, 4)
sum(round(Design_redo, 4) == round(Design, 4))

#Design_redo matches Design. 

