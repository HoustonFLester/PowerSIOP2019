#####################################################################################
# Author: Houston F. Lester                                                                
# Description: Bayesian simulation. This is only for 10 reps. This simulation can be
# pretty slow. Please see the Stan installation files for how to install it. 
                                                     
#####################################################################################


# setwd("C:/Rbatch")

library(mvtnorm)
library(tidyr)
library(brms)
library(dplyr)

npersons <- 400
mean_dif <- .2
N_reps <- 1000

beta_0 <- 0
beta_1 <- .2
beta_2 <- .2
beta_3 <- .1

mean_vec2 <- matrix(data = c(0, 0), nrow = 2, ncol = 1)
cov_mat2 <- matrix(data = c(1, 0, 0, 1), nrow = 2, ncol = 2)

####Generating single data set to precompile the model of interest 
dat <- rmvnorm(n = npersons, mean = mean_vec2, cov_mat2)
dat <- as.data.frame(dat)
names(dat) <- c("x1", "x2")
dat$int1 <- dat$x1*dat$x2
dat$resids <- rnorm(n = npersons, mean = 0, sd = sqrt(.91))
dat$y <- dat$x1*beta_1 + dat$x2*beta_2 + dat$int1*beta_3 +  dat$resids
nhst_ttest <- lm(y ~ x1 + x2 + int1, data = dat)

###Bayesian model 

bayes_form <- bf(y ~ x1 + x2 + int1)
slightly_informed_prior <- set_prior("normal(.12,.25)", class = "b", coef = "int1")
info_prior <- set_prior("normal(.12,.061)", class = "b", coef = "int1")

##Uninformed prior. 
like_lihood <- brm(bayes_form, data = dat, chains = 4, cores = 1, iter = 5000,
                   control = list(adapt_delta = .80, max_treedepth = 15), seed = 1235, save_all_pars = FALSE)

slightly_informed_prior_res <- brm(bayes_form, data = dat, chains = 4, cores = 1, iter = 5000,
                                   control = list(adapt_delta = .80, max_treedepth = 15), 
                                   prior = slightly_informed_prior, seed = 5189, save_all_pars = FALSE)

info_prior_res <- brm(bayes_form, data = dat, chains = 4, cores = 1, iter = 5000,
                      control = list(adapt_delta = .80, max_treedepth = 15), prior = info_prior, seed = 845, save_all_pars = FALSE)


Design_Bayesian <- expand.grid(effect_sizeRsq = c(0, .01), cont_v_cat = c(0, 1)) # cont is a 1. 

Design_Bayesian2 <- Design_Bayesian %>% 
  mutate(beta_0 = 0, beta_1 = if_else(cont_v_cat == 1, .2, 0.2529822), 
         beta_2 = if_else(cont_v_cat == 1, .2, 0.2529822), 
         beta_3 = if_else(effect_sizeRsq == 0, 0, if_else(cont_v_cat == 1, .1, 0.0632455)),
         residual_var = if_else(effect_sizeRsq == 0, .92, .91))


fixed_objects <- list(like_lihood = like_lihood, 
                      slightly_informed_prior_res = slightly_informed_prior_res,
                      info_prior_res = info_prior_res, mean_vec2 = mean_vec2, 
                      cov_mat2 = cov_mat2)

# Generate_bayes <- function(condition, fixed_objects = NULL){
condition <- Design_Bayesian2
N <- 400
R <- 10


out_hold_1000 <- vector(mode = "list", length = nrow(condition)*R)
seedz1 <- runif(1000*400, min = 1000, max = 1E8)
seedz2 <- runif(1000*400, min = 1000, max = 1E8)
seedz3 <- runif(1000*400, min = 1000, max = 1E8)

for(i in 1:nrow(condition)) {
  for(j in 1:R) {
    counter <- (i-1)*R + j
    ###Generate portion ----  
    # mean_vec2 <- fixed_objects$mean_vec2
    # cov_mat2 <- fixed_objects$cov_mat2
    effect_sizeRsq <- condition$effect_sizeRsq[i]
    cont_v_cat <- condition$cont_v_cat[i]
    beta_1 <- condition$beta_1[i]
    beta_2 <- condition$beta_2[i]
    beta_3 <- condition$beta_3[i]
    residual_var <- condition$residual_var[i]
    
    if(effect_sizeRsq == 0) {
      set.seed(seedz1[counter])
      dat1 <- rmvnorm(n = N, mean = mean_vec2, cov_mat2)
      dat1 <- as.data.frame(dat1)
      names(dat1) <- c("V1", "V2")
      dat1$int1 <- dat1$V1*dat1$V2
      set.seed(seedz2[counter])
      dat1$resids <- rnorm(n = N, mean = 0, sd = sqrt(residual_var))
      dat1$y <- dat1$V1*beta_1 + dat1$V2*beta_2 + dat1$int*beta_3 +  dat1$resids
      dat <- data.frame(id = 1:N, beta0 = 0, x1 = dat1[,1], beta1 = beta_1, x2 = dat1[,2], beta2 = beta_2, 
                        int1 = dat1$int1, beta3 = beta_3, y = dat1$y, 
                        residz = dat1$resids)
    }
    else {
      set.seed(seedz2[counter])
      x <- rnorm(n = N, mean = 0, sd = 1)
      g <- rbinom(n = N, size = 1, prob = .5)
      
      x1 <- x + 0*g 
      
      cat_1mod <- data.frame(x = x, g = g, x1 = x1)
      
      cat_1mod <- cat_1mod %>%
        mutate(intz = x1*g) %>%
        select(-x)
      
      X <- as.matrix(cat_1mod)
      B <- matrix(data = c(beta_1, beta_2, beta_3), nrow = 3, ncol = 1) 
      XB <- X %*% B
      
      set.seed(seedz3[counter])
      e <- rnorm(n = N, mean = 0, sd = sqrt(residual_var))
      y <- XB + e
      dat <- data.frame(id = 1:N, beta0 = 0, x1 =X[,1], beta1 = beta_1, x2 = X[,2], beta2 = beta_2, int1 = X[,3], 
                        beta3 = beta_3, y = y, 
                        residz = residual_var, inderrors = e)
    }
    
    ####Analyze portion ----- 
    
    nhst <- lm(y ~ x1 + x2 + int1, data = dat)
    nhst_ttest <- summary(nhst)
    b0_est_nhst <- nhst_ttest$coefficients[1 ,"Estimate"]
    b1_est_nhst <- nhst_ttest$coefficients[2 ,"Estimate"]
    b2_est_nhst <- nhst_ttest$coefficients[3 ,"Estimate"]
    b3_est_nhst <- nhst_ttest$coefficients[4 ,"Estimate"]
    
    b0_pval_nhst <- nhst_ttest$coefficients[1 ,"Pr(>|t|)"]
    b1_pval_nhst <- nhst_ttest$coefficients[2 ,"Pr(>|t|)"]
    b2_pval_nhst <- nhst_ttest$coefficients[3 ,"Pr(>|t|)"]
    b3_pval_nhst <- nhst_ttest$coefficients[4 ,"Pr(>|t|)"]
    
    not_informed <- update(like_lihood, newdata = dat, seed = round(seedz1[counter]/1000, 0))
    not_informed <- summary(not_informed)$fixed
    
    slightly_informed_sim <- update(slightly_informed_prior_res, newdata = dat, seed = round(seedz2[counter]/1000, 0))
    slightly_informed_sim <- summary(slightly_informed_sim)$fixed
    
    info_prior_sim <- update(info_prior_res, newdata = dat, seed = round(seedz3[counter]/1000,0))
    info_prior_sim <- summary(info_prior_sim)$fixed
    
    
    b0_est_noinfo <- not_informed[1,"Estimate"]
    b1_est_noinfo <- not_informed[2,"Estimate"]
    b2_est_noinfo <- not_informed[3,"Estimate"]
    b3_est_noinfo <- not_informed[4,"Estimate"]
    
    b0_l95_noinfo <-not_informed[1,"l-95% CI"]
    b1_l95_noinfo <-not_informed[2,"l-95% CI"]
    b2_l95_noinfo <-not_informed[3,"l-95% CI"]
    b3_l95_noinfo <-not_informed[4,"l-95% CI"]
    
    b0_u95_noinfo <-not_informed[1,"u-95% CI"]
    b1_u95_noinfo <-not_informed[2,"u-95% CI"]
    b2_u95_noinfo <-not_informed[3,"u-95% CI"]
    b3_u95_noinfo <-not_informed[4,"u-95% CI"]
    
    b0_eff_noinfo <-not_informed[1,"Eff.Sample"]
    b1_eff_noinfo <-not_informed[2,"Eff.Sample"]
    b2_eff_noinfo <-not_informed[3,"Eff.Sample"]
    b3_eff_noinfo <-not_informed[4,"Eff.Sample"]
    
    b0_Rhat_noinfo <-not_informed[1,"Rhat"]
    b1_Rhat_noinfo <-not_informed[2,"Rhat"]
    b2_Rhat_noinfo <-not_informed[3,"Rhat"]
    b3_Rhat_noinfo <-not_informed[4,"Rhat"]
    
    ####slightly informed
    b0_est_slight<- slightly_informed_sim[1,"Estimate"]
    b1_est_slight<- slightly_informed_sim[2,"Estimate"]
    b2_est_slight<- slightly_informed_sim[3,"Estimate"]
    b3_est_slight<- slightly_informed_sim[4,"Estimate"]
    
    b0_l95_slight<- slightly_informed_sim[1,"l-95% CI"]
    b1_l95_slight<- slightly_informed_sim[2,"l-95% CI"]
    b2_l95_slight<- slightly_informed_sim[3,"l-95% CI"]
    b3_l95_slight<- slightly_informed_sim[4,"l-95% CI"]
    
    b0_u95_slight<- slightly_informed_sim[1,"u-95% CI"]
    b1_u95_slight<- slightly_informed_sim[2,"u-95% CI"]
    b2_u95_slight<- slightly_informed_sim[3,"u-95% CI"]
    b3_u95_slight<- slightly_informed_sim[4,"u-95% CI"]
    
    b0_eff_slight<- slightly_informed_sim[1,"Eff.Sample"]
    b1_eff_slight<- slightly_informed_sim[2,"Eff.Sample"]
    b2_eff_slight<- slightly_informed_sim[3,"Eff.Sample"]
    b3_eff_slight<- slightly_informed_sim[4,"Eff.Sample"]
    
    b0_Rhatslight <- slightly_informed_sim[1,"Rhat"]
    b1_Rhatslight <- slightly_informed_sim[2,"Rhat"]
    b2_Rhatslight <- slightly_informed_sim[3,"Rhat"]
    b3_Rhatslight <- slightly_informed_sim[4,"Rhat"]
    
    ##informative 
    
    b0_est_info <- info_prior_sim[1,"Estimate"]
    b1_est_info <- info_prior_sim[2,"Estimate"]
    b2_est_info <- info_prior_sim[3,"Estimate"]
    b3_est_info <- info_prior_sim[4,"Estimate"]
    
    b0_l95_info <- info_prior_sim[1,"l-95% CI"]
    b1_l95_info <- info_prior_sim[2,"l-95% CI"]
    b2_l95_info <- info_prior_sim[3,"l-95% CI"]
    b3_l95_info <- info_prior_sim[4,"l-95% CI"]
    
    b0_u95_info <- info_prior_sim[1,"u-95% CI"]
    b1_u95_info <- info_prior_sim[2,"u-95% CI"]
    b2_u95_info <- info_prior_sim[3,"u-95% CI"]
    b3_u95_info <- info_prior_sim[4,"u-95% CI"]
    
    b0_eff_info <- info_prior_sim[1,"Eff.Sample"]
    b1_eff_info <- info_prior_sim[2,"Eff.Sample"]
    b2_eff_info <- info_prior_sim[3,"Eff.Sample"]
    b3_eff_info <- info_prior_sim[4,"Eff.Sample"]
    
    b0_Rhatinfo <-info_prior_sim[1,"Rhat"]
    b1_Rhatinfo <-info_prior_sim[2,"Rhat"]
    b2_Rhatinfo <-info_prior_sim[3,"Rhat"]
    b3_Rhatinfo <-info_prior_sim[4,"Rhat"]
    
    ret <- c(design_row = i, effect_sizeRsq = effect_sizeRsq, cont_v_cat = cont_v_cat, beta_1 = beta_1, beta_2 = beta_2,
             beta_3 = beta_3, residual_var = residual_var,
             replication = j, b0_est_nhst = b0_est_nhst, b0_pval_nhst = b0_pval_nhst, b1_est_nhst = b1_est_nhst, b2_est_nhst = b2_est_nhst, b3_est_nhst = b3_est_nhst, b3_pval_nhst = b3_pval_nhst, 
             b0_est_noinfo = b0_est_noinfo, b1_est_noinfo = b1_est_noinfo, b2_est_noinfo = b2_est_noinfo, b3_est_noinfo = b3_est_noinfo, 
             b0_l95_noinfo = b0_l95_noinfo, b1_l95_noinfo = b1_l95_noinfo, b2_l95_noinfo = b2_l95_noinfo, b3_l95_noinfo = b3_l95_noinfo, 
             b0_u95_noinfo = b0_u95_noinfo, b1_u95_noinfo = b1_u95_noinfo, b2_u95_noinfo = b2_u95_noinfo, b3_u95_noinfo = b3_u95_noinfo, 
             b0_eff_noinfo = b0_eff_noinfo, b1_eff_noinfo = b1_eff_noinfo, b2_eff_noinfo = b2_eff_noinfo, b3_eff_noinfo = b3_eff_noinfo, 
             b0_Rhat_noinfo = b0_Rhat_noinfo, b1_Rhat_noinfo = b1_Rhat_noinfo, b2_Rhat_noinfo = b2_Rhat_noinfo, b3_Rhat_noinfo = b3_Rhat_noinfo, 
             b0_est_slight =  b0_est_slight, b1_est_slight =  b1_est_slight, b2_est_slight =  b2_est_slight, b3_est_slight =  b3_est_slight, 
             b0_l95_slight =  b0_l95_slight, b1_l95_slight =  b1_l95_slight, b2_l95_slight =  b2_l95_slight, b3_l95_slight =  b3_l95_slight, 
             b0_u95_slight =  b0_u95_slight, b1_u95_slight =  b1_u95_slight, b2_u95_slight =  b2_u95_slight, b3_u95_slight =  b3_u95_slight, 
             b0_eff_slight =  b0_eff_slight, b1_eff_slight =  b1_eff_slight, b2_eff_slight =  b2_eff_slight, b3_eff_slight =  b3_eff_slight, 
             b0_Rhatslight =  b0_Rhatslight, b1_Rhatslight =  b1_Rhatslight, b2_Rhatslight =  b2_Rhatslight, b3_Rhatslight =  b3_Rhatslight, 
             b0_est_info = b0_est_info, b1_est_info = b1_est_info, b2_est_info = b2_est_info, b3_est_info = b3_est_info, b0_l95_info = b0_l95_info, 
             b1_l95_info = b1_l95_info, b2_l95_info = b2_l95_info, b3_l95_info = b3_l95_info, b0_u95_info = b0_u95_info, 
             b1_u95_info = b1_u95_info, b2_u95_info = b2_u95_info, b3_u95_info = b3_u95_info, b0_eff_info = b0_eff_info, 
             b1_eff_info = b1_eff_info, b2_eff_info = b2_eff_info, b3_eff_info = b3_eff_info, b0_Rhatinfo = b0_Rhatinfo, 
             b1_Rhatinfo = b1_Rhatinfo, b2_Rhatinfo = b2_Rhatinfo, b3_Rhatinfo = b3_Rhatinfo)
    out_hold_1000[[(i-1)*R + j]] <- ret
  }
}

out <- do.call(rbind, out_hold_1000)
out <- as.data.frame(out)


#####Creating the plot for a single analyzed replication. 

###You can generate your own data to see how playing with sample size affects
###the distance between the likelihood and the posterior. 
###This is just using the last dataset that was generated from the simulation. The
###figure displayed on the poster is just illustrative and uses a smaller sample 
###size than the figure displayed here. 



for_plot_med_post <- brm(bayes_form, data = dat, chains = 4, cores = 4, iter = 5000,
                         control = list(adapt_delta = .99, max_treedepth = 15), prior = slightly_informed_prior, 
                         seed = 875, save_all_pars = T)

like_lihood <- brm(bayes_form, data = dat, chains = 4, cores = 4, iter = 5000,
                         control = list(adapt_delta = .80, max_treedepth = 15), prior = slightly_informed_prior, 
                         seed = 147, save_all_pars = T)

for_plot_med_post_samples <- posterior_samples(for_plot_med_post)
plot_likelihood <- posterior_samples(like_lihood)
plot_prior_med <- rnorm(n = 10000, mean = .12, sd = .061)

plot_wide_med <- data.frame(Posterior = for_plot_med_post_samples$b_int, Likelihood = plot_likelihood$b_int,
                            Prior = plot_prior_med)

plot_med_long_data <- gather(plot_wide_med, var_label, samples)

ggplot(plot_med_long_data,aes(x=samples, fill=var_label)) + 
  geom_density(alpha=0.25) +
  theme_bw() +
  xlab("Sampled Parameter Values") +
  ylab("Density") +
  theme(legend.title=element_blank())












