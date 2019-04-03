###Continuous simulation study Results section. 

library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)

###General setup ----
results_to_plot <- results_info_hypo_1000reps



# View(Design)

# N <- c(500, 1000, 2000)
# f_sq <- c(0, .002, .01, .02)
# n_inter <- c(1, 2, 4)
# continuous <- c(0, 1) # 1 indicates continuous moderator
# #The way that it is analyzed is the only variable that is not on a different row. So, I think that it
# #is the only one that needs to be flipped around. 

results_to_plot$N <- as.factor(results_to_plot$N)
results_to_plot$f_sq <- as.factor(results_to_plot$f_sq)
results_to_plot$n_inter <- as.factor(results_to_plot$n_inter)

to_plot_infoHYPO <- results_to_plot %>%
  select(N, f_sq, n_inter,  continuous, TypeJ_power, TypeJ_power_other, power_mod_comp) %>%
  mutate(orig_order = 1:36)


longto_plot_infoHYPO <- to_plot_infoHYPO %>%
  gather(key = analysis_method, value = Power, TypeJ_power, power_mod_comp)

longto_plot_infoHYPO <- arrange(longto_plot_infoHYPO, orig_order)


##Plot will likely look much better as Type I error rates and power. 

###Type I error plot -----

TypeI <- longto_plot_infoHYPO %>%
  filter(f_sq == 0) %>%
  mutate(analysis_method_pretty = if_else(analysis_method == "TypeJ_power", 
                                          "Informative Hypothesis", "NHST"))



TypeI_error_plot <- ggplot(TypeI, aes(x = N, y = Power, color = analysis_method_pretty)) +
  geom_point(aes(shape = analysis_method_pretty)) +
  ylim(c(.03, .07)) +
  ylab("Type I Error") +
  facet_grid(rows = vars(n_inter)) +
  theme_bw() +
  scale_color_grey() + 
  theme(legend.title=element_blank())

###Power plot -----

power <- longto_plot_infoHYPO %>%
  filter(f_sq != 0) %>%
  mutate(analysis_method_pretty = if_else(analysis_method == "TypeJ_power", 
                                          "Informative Hypothesis", "NHST"),
         f_sq_num = as.numeric(f_sq))


#Getting the column labels to look good is annoying and confusing. 
testxx <- paste0("f^2 ==", levels(power$f_sq)[2:4])

power$fsq_better <- as.factor(power$f_sq_num)
levels(power$fsq_better) <- testxx

power_plot_cont <- ggplot(power, aes(x = N, y = Power, color = analysis_method_pretty)) +
  geom_point(aes(shape = analysis_method_pretty)) +
  facet_grid(cols = vars(fsq_better), labeller = label_parsed, 
             rows = vars(n_inter)) +
  geom_hline(yintercept = .8, linetype = 2) +
  theme_bw() +
  scale_color_grey() +
  theme(legend.title=element_blank(), legend.position = "left") + 
  ggtitle("Continuous Results") +
  theme(plot.title = element_text(hjust = 0.5))
  







#You get some weird results here. Like it is misreading the effect size value. 
# ggplot(power, aes(x = N, y = Power, color = analysis_method_pretty)) +
#   geom_point(aes(shape = analysis_method_pretty)) +
#   facet_grid(cols = vars(f_sq), labeller = label_bquote(cols = f^2 == .(f_sq)), 
#              rows = vars(n_inter)) +
#   theme_bw() +
#   scale_color_grey() +
#   theme(legend.title=element_blank())
