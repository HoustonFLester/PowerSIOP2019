#####################################################################################
# Author: Houston F. Lester                                                                
# Description: this is the expanded version of the Continuous predictors simulation
# using the informative hypothesis method
# If none of these packages are installed, then you can use the commented code below.                                                           
#####################################################################################


# install.packages("mvtnorm")
# install.packages("dplyr")
# install.packages("SimDesign")
# install.packages("restriktor")

library(mvtnorm)
library(dplyr)
library(SimDesign)
library(restriktor)

####Three different covariance matrices associated with the data generating populations. You don't actually use these to generate the data.
####These are theoretical covariances (calculated by using expected value analytical calculations instead of simulating a large sample 
####and calcuating it using the sample statistics. This is helpful for knowing exactly what the R squared values are, but practically 
####speaking there is not difference (differnces < .001) between the observed and theoretical. 
  

mean_vec3 <- matrix(data = c(0, 0, 0),  nrow = 3, ncol = 1)
cov_mat_int3_cont <- matrix(data = c(1, .3, 0, 
                                     .3, 1, 0, 0,
                                     0, 1.09), nrow = 3, ncol = 3)



##Corresponds to the 2 interaction terms conditions (i.e., 2 two way interaction model)
mean_vec5 <- matrix(data = c(0, 0, 0, 0, 0),  nrow = 5, ncol = 1)
cov_mat_int5_cont <- matrix(data = c(1, .3, .3, 0, 0, 
                                     .3, 1, .3, 0, 0, 
                                     .3, .3, 1, 0, 0,
                                     0, 0, 0, 1.09, .39,
                                     0, 0, 0, .39, 1.09), nrow = 5, ncol = 5) 


##Corresponds to the 4 interaction terms conditions (i.e., 3 way interaction model)
mean_vec7 <- matrix(data = c(0, 0, 0, 0, 0, 0, 0),  nrow = 7, ncol = 1)
cov_mat_int7_cont <- matrix(data = c(1, .3, .3, 0, 0, 0, .48,  
                                     .3, 1, .3, 0, 0, 0, .48, 
                                     .3, .3, 1, 0, 0, 0, .48,
                                     0, 0, 0, 1.09, .39, .39, 0,
                                     0, 0, 0, .39, 1.09, .39, 0,
                                     0, 0, 0, .39, .39, 1.09, 0,
                                     .48, .48, .48, 0, 0, 0, 1.756), nrow = 7, ncol = 7) 



#### Adding in the actual data generating stuff. 

mvec2 <- matrix(data = c(0,0), nrow = 2, ncol = 1)
mvec3 <- matrix(data = c(0, 0, 0), nrow = 3, ncol = 1)

cov2 <- matrix(data = c(1, .3, .3, 1), nrow = 2, ncol = 2)
cov3 <- matrix(data = c(1, .3, .3, 
                        .3, 1, .3,
                        .3, .3, 1), nrow = 3, ncol = 3)


data_gen_stuff <- list(mvec2 = mvec2, mvec3 = mvec3, cov2 = cov2, cov3 = cov3)

N <- c(500, 1000, 2000)
f_sq <- c(0, .002, .01, .02)
n_inter <- c(1, 2, 4)
continuous <- c(0, 1) # 1 indicates continuous moderator
Nreps <- 1000 # Change this to number of replications that you want. For the poster it was 1000. 

const_mod1 <- " int1 > 0"

const_mod2 <- " int1 > 0 
int2 > 0"

const_three_way <- "int1 > 0
int2 > 0
int3 > 0
three_way > 0"

# Design <- expand.grid(N, f_sq, n_inter, continuous)
# names(Design) <- c("N", "f_sq", "n_inter", "continous")
## Additional calculations were performed elsewhere to arrive at the design matrix that is imported below. 

# The path below is relative to the working directory which would be where the Rproject is located on your computer. I don't think you should
# have to adjust anything 
Design <- read.csv("Data/Power_ORM_Design.csv")




# I used this, but it is only necessary if you are using the parallel processing abilities built into the SimDesign package. It can dramatically speed things
# up if your computer has a lot of cores. I did not do parallel processing here, because it is pretty fast anyway. 
# fixed_objects <- list(data_gen_stuff = data_gen_stuff,
#                       const_mod1 = const_mod1,
#                       const_mod2 = const_mod2,
#                       const_three_way = const_three_way)

#### Generating the data that will be used analyzed and summarized. 
Generate <- function(condition, fixed_objects = NULL){
  npred_full <- condition$npred_full
  N <- condition$N
  beta1 <- condition$beta1
  beta2 <- condition$beta2
  beta3 <- condition$beta
  error_var <- (1 - condition$rsq_full)
  if(npred_full == 3) {
    X <- rmvnorm(n = N, mean = data_gen_stuff$mvec2, 
                 data_gen_stuff$cov2)
    X <- data.frame(X)
    X$int1 <- X[,1]*X[,2]
    X <- as.matrix(X)
    B <- matrix(data = c(condition$beta1, condition$beta2, condition$int1), nrow = 3, ncol = 1) 
    XB <- X %*% B
    e <- rnorm(n = N, mean = 0, sd = sqrt(error_var))
    y <- XB + e
    dat <- data.frame(id = 1:N, beta0 = 0, x1 =X[,1], beta1 = beta1, x2 = X[,2], beta2 = beta2, int1 = X[,3], y = y, residz = error_var)
  }
  else if (npred_full == 5) {
    X <- rmvnorm(n = N, mean = data_gen_stuff$mvec3, 
                 data_gen_stuff$cov3)
    X <- as.data.frame(X)
    B <- matrix(data = c(condition$beta1, condition$beta2, condition$beta3, condition$int1, condition$int2), nrow = 5, ncol = 1) 
    X$int1 <- X[,1]*X[,2]
    X$int2 <- X[,1]*X[,3]
    X <- as.matrix(X)
    XB <- X %*% B
    e <- rnorm(n = N, mean = 0, sd = sqrt(error_var))
    y <- XB + e
    dat <- data.frame(id = 1:N, beta0 = 0, x1 = X[,1], beta1 = beta1, x2 = X[,2], beta2 = beta2, x3 = X[,3], beta3 = beta3,
                      int1 = X[,4], int2 = X[,5], y = y, residz = error_var)
  }
  else {
    X <- rmvnorm(n = N, mean = data_gen_stuff$mvec3, 
                 data_gen_stuff$cov3)
    X <- as.data.frame(X)
    B <- matrix(data = c(condition$beta1, condition$beta2, condition$beta3, condition$int1, condition$int2, 
                         condition$int3, condition$three_way), nrow = 7, ncol = 1) 
    X$int1 <- X[,1]*X[,2]
    X$int2 <- X[,1]*X[,3]
    X$int3 <- X[,2]*X[,3]
    X$int4 <- X[,1]*X[,2]*X[,3]
    X <- as.matrix(X)
    XB <- X %*% B
    e <- rnorm(n = N, mean = 0, sd = sqrt(error_var))
    y <- XB + e
    dat <- data.frame(id = 1:N, beta0 = 0, x1 = X[,1], beta1 = beta1, x2 = X[,2], beta2 = beta2, x3 = X[,3], beta3 = beta3,
                      int1 = X[,4], int2 = X[,5], int3 = X[,6], three_way = X[,7], y = y, residz = error_var)
  }
  dat
}


#### Analyzing the generating data. 
Analyse <- function(condition, dat, fixed_objects = NULL){
  npred_full <- condition$npred_full
  if(npred_full == 3) {
    red_1mod <- lm(y ~  x1 + x2 , data = dat)
    full_1mod <- lm(y ~  x1 + x2 + int1, data = dat)
    anova1INT <- anova(red_1mod, full_1mod)
    full_v_red <- anova1INT$`Pr(>F)`[2]
    fixef_pvalues <- summary(full_1mod)
    fixef_pvalues <- coef(fixef_pvalues)
    ##### Inequality constraint analysis 
    inequal_const_out <- restriktor::conTest(full_1mod, constraints = const_mod1)
    A_pvalue <- inequal_const_out$A$pvalue[1]
    B_pvalue <- inequal_const_out$B$pvalue[1]
    J_pvalue <- dplyr::if_else(B_pvalue >= .05, A_pvalue, 1)
    J_pvalue_other_way <- (B_pvalue >=.05)*(A_pvalue < .05)
  
    ##### Typical null hypothesis significance testing (NHST)
    beta0 <- fixef_pvalues[1, 1]
    beta1 <- fixef_pvalues[2, 1]
    beta2 <- fixef_pvalues[3, 1]
    beta0_p <- fixef_pvalues[1, 4]
    beta1_p <- fixef_pvalues[2, 4]
    beta2_p <- fixef_pvalues[3, 4]
    int1 <- fixef_pvalues[4, 1]
    int1_p <- fixef_pvalues[4, 4]
    beta3 <- -999
    beta3_p <- -999
    int2 <- -999
    int3 <- -999
    int3way <- -999
    int2_p <- -999
    int3_p <- -999
    int3way_p <- -999
    ret <- c(full_v_red = full_v_red, int1 = int1, int1_p = int1_p, 
             int2 = int2, int2_p = int2_p,
             int3 = int3, int3_p = int3_p,
             int3way = int3way, int3way_p = int3way_p,
             beta0 = beta0, beta0_p = beta0_p, 
             beta1 = beta1, beta1_p = beta1_p,
             beta2 = beta2, beta2_p = beta2_p,
             beta3 = beta3, beta3_p = beta3_p,
             A_pvalue = A_pvalue, B_pvalue = B_pvalue,
             J_pvalue = J_pvalue, J_pvalue_other_way = J_pvalue_other_way)
  }
  else if (npred_full == 5) {
    red_2mod <- lm(y ~  x1 + x2 + x3, data = dat)
    full_2mod <- lm(y ~  x1 + x2 + x3 + int1 + int2, data = dat)
    anova2INT <- anova(red_2mod, full_2mod) 
    full_v_red <- anova2INT$`Pr(>F)`[2]
    fixef_pvalues <- summary(full_2mod)
    fixef_pvalues <- coef(fixef_pvalues)
    
    ##### Inequality constraint analysis 
    inequal_const_out <- restriktor::conTest(full_2mod, constraints = const_mod2)
    A_pvalue <- inequal_const_out$A$pvalue[1]
    B_pvalue <- inequal_const_out$B$pvalue[1]
    J_pvalue <- dplyr::if_else(B_pvalue >= .05, A_pvalue, 1)
    J_pvalue_other_way <- (B_pvalue >=.05)*(A_pvalue < .05)
    
    ##### Typical null hypothesis significance testing (NHST)
    beta0 <- fixef_pvalues[1, 1]
    beta1 <- fixef_pvalues[2, 1]
    beta2 <- fixef_pvalues[3, 1]
    beta3 <- fixef_pvalues[4, 1]
    beta0_p <- fixef_pvalues[1, 4]
    beta1_p <- fixef_pvalues[2, 4]
    beta2_p <- fixef_pvalues[3, 4]
    beta3_p <- fixef_pvalues[4, 4]
    int1 <- fixef_pvalues[5, 1]
    int1_p <- fixef_pvalues[5, 4]
    int2 <- fixef_pvalues[6, 1]
    int2_p <- fixef_pvalues[6, 4]
    int3 <- -999
    int3way <- -999
    int3_p <- -999
    int3way_p <- -999
    ret <- c(full_v_red = full_v_red, int1 = int1, int1_p = int1_p, 
             int2 = int2, int2_p = int2_p,
             int3 = int3, int3_p = int3_p,
             int3way = int3way, int3way_p = int3way_p,
             beta0 = beta0, beta0_p = beta0_p, 
             beta1 = beta1, beta1_p = beta1_p,
             beta2 = beta2, beta2_p = beta2_p,
             beta3 = beta3, beta3_p = beta3_p,
             A_pvalue = A_pvalue, B_pvalue = B_pvalue,
             J_pvalue = J_pvalue, J_pvalue_other_way = J_pvalue_other_way)
  }
  else {
    red_3mod <- lm(y ~  x1 + x2 + x3, data = dat)
    full_3mod <- lm(y ~  x1 + x2 + x3 + int1 + int2 + int3 + three_way, data = dat)
    anova3INT <- anova(red_3mod, full_3mod) 
    full_v_red <- anova3INT$`Pr(>F)`[2]
    fixef_pvalues <- summary(full_3mod)
    fixef_pvalues <- coef(fixef_pvalues)
    
    
    ##### Inequality constraint analysis 
    inequal_const_out <- restriktor::conTest(full_3mod, constraints = const_three_way)
    A_pvalue <- inequal_const_out$A$pvalue[1]
    B_pvalue <- inequal_const_out$B$pvalue[1]
    J_pvalue <- dplyr::if_else(B_pvalue >= .05, A_pvalue, 1)
    J_pvalue_other_way <- (B_pvalue >=.05)*(A_pvalue < .05)
    
    beta0 <- fixef_pvalues[1, 1]
    beta1 <- fixef_pvalues[2, 1]
    beta2 <- fixef_pvalues[3, 1]
    beta3 <- fixef_pvalues[4, 1]
    beta0_p <- fixef_pvalues[1, 4]
    beta1_p <- fixef_pvalues[2, 4]
    beta2_p <- fixef_pvalues[3, 4]
    beta3_p <- fixef_pvalues[4, 4]
    int1 <- fixef_pvalues[5, 1]
    int1_p <- fixef_pvalues[5, 4]
    int2 <- fixef_pvalues[6, 1]
    int2_p <- fixef_pvalues[6, 4]
    int3 <- fixef_pvalues[7, 1]
    int3way <- fixef_pvalues[8, 1]
    int3_p <- fixef_pvalues[7, 4]
    int3way_p <- fixef_pvalues[8, 4]
    ret <- c(full_v_red = full_v_red, 
             int1 = int1, 
             int1_p = int1_p, 
             int2 = int2,
             int2_p = int2_p,
             int3 = int3,
             int3_p = int3_p,
             int3way = int3way,
             int3way_p = int3way_p,
             beta0 = beta0,
             beta0_p = beta0_p,
             beta1 = beta1,
             beta1_p = beta1_p,
             beta2 = beta2,
             beta2_p = beta2_p,
             beta3 = beta3,
             beta3_p = beta3_p, 
             A_pvalue = A_pvalue,
             B_pvalue = B_pvalue,
             J_pvalue = J_pvalue, J_pvalue_other_way = J_pvalue_other_way)
  }
  ret
}

Summarise <- function(condition, results, fixed_objects = NULL){
  power_mod_comp <- EDR(results[, 1], alpha = .05)
  int1_mean <- mean(results[,2], na.rm = T)
  int1_bias <- bias(results[,2], parameter = condition$int1, type = 'relative')
  int1_power <- EDR(results[,3], alpha = .05)
  int2_mean <- mean(results[,4])
  int2_bias <- bias(results[,4], parameter = condition$int2, type = 'relative')
  TypeA_power <- EDR(results[,18], .05)
  TypeB_power <- EDR(results[,19], .05)
  TypeJ_power <- EDR(results[,20], .05)
  TypeJ_power_other <- mean(results[,21])
  
  
  ret <- c(power_mod_comp = power_mod_comp, int1_mean = int1_mean, int1_bias = int1_bias,
           int1_power = int1_power, int2_mean = int2_mean, int2_bias = int2_bias,
           TypeA_power = TypeA_power, TypeB_power = TypeB_power, 
           TypeJ_power = TypeJ_power, TypeJ_power_other = TypeJ_power_other)
  ret
}

#These results appear to work. I need to go back through and make sure that he did not do anything else weird. 

set.seed(123)
seeds <- round(runif(n = 36, min = 1, max = 1E3), 0)
seeds
#Using parallel computing. 
# results_1000_new <- runSimulation(Design, replications = 1000, seed = seeds, 
#                                   generate = Generate, analyse = Analyse, 
#                                   summarise = Summarise, parallel = T, 
#                                   ncores = 51)

# cl <- parallel::makeCluster(type='PSOCK', master=primary, spec=spec)


results_info_hypo_1000reps <- runSimulation(Design, replications = 1000, seed = seeds, 
                                generate = Generate, analyse = Analyse, 
                                summarise = Summarise)








