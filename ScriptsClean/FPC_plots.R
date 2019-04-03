#####################################################################################
# Author: Houston F. Lester                                                                
# Description: This is the plot that was created for the finite population correction                                                           
#####################################################################################


library(ggplot2)
library(dplyr)
library(tidyr)

FPC <- function(k, K){
  mult <- 1-k/K
  mult
}

#assuming residual variance is one and the corresponding diagonal element of solve(t(X)%*%X) is also 1
mult_factor <- FPC(1:99, 100)
FPC_plot <- data.frame(mult_factor, sampled = 1:99, t_stat = 1/mult_factor)

ggplot(data = FPC_plot, aes(x = sampled, y = mult_factor)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylab("Squared Standard Error") +
  xlab("Percentage of the Population Sampled") +
  geom_vline(xintercept = 75) +
  ggtitle("Top Panel")

ggplot(data = FPC_plot, aes(x = sampled, y = t_stat)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylab("t/z Statistic") +
  xlab("Percentage of the Population Sampled") +
  geom_vline(xintercept = 75) +
  ggtitle("Bottom Panel")
