rm(list = ls())
setwd("C:/Users/Jun Cai/Documents/C++/VOC_burden")

library(tidyverse)
library(readr)
library(glue)

files <- list.files("output/single/leaky", pattern = "*.txt", full.names = T)
file <- files[1]
newI.dose.sim <- read_delim(file, delim = " ")

newI.gr <- newI.dose.sim %>% 
  select(sim, t, starts_with("newI."))

dose2.gr <- newI.dose.sim %>% 
  select(sim, t, starts_with("dose2."))

# extract parameters
R0 <- as.numeric(str_extract(file, "(?<=R0-).*(?=_ve2)"))
ve2 <- as.numeric(str_extract(file, "(?<=ve2-).*(?=_Tonset)"))
Tonset <- as.integer(str_extract(file, "(?<=Tonset-).*(?=_strategy)"))
strategy <- as.integer(str_extract(file, "(?<=strategy-).*(?=_nsim)"))
nsim <- as.integer(str_extract(file, "(?<=nsim-).*(?=.txt)"))

# check epidemic curve
N <- read.table("data/population") %>%
  unlist(use.names = FALSE)
(N.gr1 <- sum(N[1:3]))
(N.gr2 <- sum(N[4:11]))
(N.gr3 <- sum(N[12:16]))

newI.gr1 <- newI.gr %>% 
  mutate(newI = rowSums(select(., starts_with("newI"))), 
         newI.gr1 = rowSums(select(., ends_with(glue("newI.{1:3}")))), 
         newI.gr2 = rowSums(select(., ends_with(glue("newI.{4:11}")))), 
         newI.gr3 = rowSums(select(., ends_with(glue("newI.{12:16}"))))) %>% 
  group_by(sim) %>% 
  mutate(cumI.prop = cumsum(newI) / sum(N), 
         cumI.prop.gr1 = cumsum(newI.gr1) / sum(N.gr1), 
         cumI.prop.gr2 = cumsum(newI.gr2) / sum(N.gr2), 
         cumI.prop.gr3 = cumsum(newI.gr3) / sum(N.gr3)) %>% 
  ungroup() %>% 
  group_by(t) %>% 
  summarise(newI.m = round(mean(newI)), 
            newI.lwr = round(quantile(newI, probs = 0.025)), 
            newI.upr = round(quantile(newI, probs = 0.975)), 
            newI.gr1.m = round(mean(newI.gr1)), 
            newI.gr1.lwr = round(quantile(newI.gr1, probs = 0.025)), 
            newI.gr1.upr = round(quantile(newI.gr1, probs = 0.975)), 
            newI.gr2.m = round(mean(newI.gr2)), 
            newI.gr2.lwr = round(quantile(newI.gr2, probs = 0.025)), 
            newI.gr2.upr = round(quantile(newI.gr2, probs = 0.975)), 
            newI.gr3.m = round(mean(newI.gr3)), 
            newI.gr3.lwr = round(quantile(newI.gr3, probs = 0.025)), 
            newI.gr3.upr = round(quantile(newI.gr3, probs = 0.975)),
            cumI.prop.m = mean(cumI.prop), 
            cumI.prop.lwr = quantile(cumI.prop, probs = 0.025), 
            cumI.prop.upr = quantile(cumI.prop, probs = 0.975), 
            cumI.prop.gr1.m = mean(cumI.prop.gr1), 
            cumI.prop.gr1.lwr = quantile(cumI.prop.gr1, probs = 0.025), 
            cumI.prop.gr1.upr = quantile(cumI.prop.gr1, probs = 0.975), 
            cumI.prop.gr2.m = mean(cumI.prop.gr2), 
            cumI.prop.gr2.lwr = quantile(cumI.prop.gr2, probs = 0.025), 
            cumI.prop.gr2.upr = quantile(cumI.prop.gr2, probs = 0.975), 
            cumI.prop.gr3.m = mean(cumI.prop.gr3), 
            cumI.prop.gr3.lwr = quantile(cumI.prop.gr3, probs = 0.025), 
            cumI.prop.gr3.upr = quantile(cumI.prop.gr3, probs = 0.975))

# plot epidemic curve by age group and strategy ####
p1 <- ggplot(data = subset(newI.gr1, t >= Tonset)) + 
  geom_ribbon(aes(x = t, ymin = newI.gr1.lwr, ymax = newI.gr1.upr, fill = "gr1"), 
              alpha = 0.3) + 
  geom_line(aes(t, newI.gr1.m, color = "gr1")) + 
  geom_ribbon(aes(x = t, ymin = newI.gr2.lwr, ymax = newI.gr2.upr, fill = "gr2"), 
              alpha = 0.3) + 
  geom_line(aes(t, newI.gr2.m, color = "gr2")) + 
  geom_ribbon(aes(x = t, ymin = newI.gr3.lwr, ymax = newI.gr3.upr, fill = "gr3"), 
              alpha = 0.3) + 
  geom_line(aes(t, newI.gr3.m, color = "gr3")) + 
  geom_ribbon(aes(x = t, ymin = newI.lwr, ymax = newI.upr, fill = "All"), 
              alpha = 0.3) + 
  geom_line(aes(t, newI.m, color = "All")) + 
  geom_vline(xintercept = Tonset, linetype = "dashed", color = "gray") + 
  labs(x = "Days", y = "Number of new infections") + 
  scale_x_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) + 
  scale_color_manual(values = c(
    "All" = "#FC4E07",
    "gr1" = "#6698FF",
    "gr2" = "#00AFBB",
    "gr3" = "#E7B800"
  ), name = "Age group") + 
  scale_fill_manual(values = c(
    "All" = "#FC4E07",
    "gr1" = "#6698FF",
    "gr2" = "#00AFBB",
    "gr3" = "#E7B800"
  ), name = "Age group") + 
  theme_classic() + 
  theme(
    legend.position = c(0.9, 0.9), 
    legend.title = element_blank(), 
    legend.text = element_text(size = 11, color = "black"),
    strip.background = element_rect(color = NA), 
    strip.text = element_text(size = 13), 
    panel.grid.minor = element_line(colour = "gray", size = 0.1), 
    axis.text = element_text(size = 13, color = "black"),
    axis.ticks = element_line(color = "darkgray"), 
    axis.title = element_text(size = 13, color = "black", face = "bold")
  )
p1

outfile <- glue("figs/single/leaky/epi_curve_R0-{R0}_ve2-{ve2}_Tonset-{Tonset}_strategy-{strategy}_nsim-{nsim}.tif")
tiff(outfile, width = 12, height = 10, unit = "in", res = 300, compression = "lzw")
print(p1)
dev.off()


# plot infected proportion by age group ####
# final attack rate by age group
far <- newI.gr1 %>% 
  select(cumI.prop.m, cumI.prop.gr1.m, cumI.prop.gr2.m, cumI.prop.gr3.m) %>% 
  slice(n()) %>% 
  mutate(far.lab = format(cumI.prop.m, digits = 3, nsmall = 3), 
         far.gr1.lab = format(cumI.prop.gr1.m, digits = 3, nsmall = 3), 
         far.gr2.lab = format(cumI.prop.gr2.m, digits = 3, nsmall = 3), 
         far.gr3.lab = format(cumI.prop.gr3.m, digits = 3, nsmall = 3)) %>% 
  mutate(x = Tonset + 200, 
         gr1.x = Tonset + 250, 
         gr2.x = Tonset + 300, 
         gr3.x = Tonset + 350)

p2 <- ggplot(data = subset(newI.gr1, t >= Tonset)) + 
  geom_ribbon(aes(x = t, ymin = cumI.prop.gr1.lwr, ymax = cumI.prop.gr1.upr, fill = "gr1"), 
              alpha = 0.3) + 
  geom_line(aes(t, cumI.prop.gr1.m, color = "gr1")) + 
  geom_ribbon(aes(x = t, ymin = cumI.prop.gr2.lwr, ymax = cumI.prop.gr2.upr, fill = "gr2"), 
              alpha = 0.3) + 
  geom_line(aes(t, cumI.prop.gr2.m, color = "gr2")) + 
  geom_ribbon(aes(x = t, ymin = cumI.prop.gr3.lwr, ymax = cumI.prop.gr3.upr, fill = "gr3"), 
              alpha = 0.3) + 
  geom_line(aes(t, cumI.prop.gr3.m, color = "gr3")) + 
  geom_ribbon(aes(x = t, ymin = cumI.prop.lwr, ymax = cumI.prop.upr, fill = "All"), 
              alpha = 0.3) + 
  geom_line(aes(t, cumI.prop.m, color = "All")) + 
  geom_vline(xintercept = Tonset, linetype = "dashed", color = "gray") + 
  geom_text(data = far, aes(x = gr1.x, y = cumI.prop.gr1.m + 0.05, 
                            label = far.gr1.lab), color = "#6698FF") + 
  geom_text(data = far, aes(x = gr2.x, y = cumI.prop.gr2.m + 0.05, 
                            label = far.gr2.lab), color = "#00AFBB") + 
  geom_text(data = far, aes(x = gr3.x, y = cumI.prop.gr3.m + 0.05, 
                            label = far.gr3.lab), color = "#E7B800") + 
  geom_text(data = far, aes(x = x, y = cumI.prop.m + 0.05, 
                            label = far.lab), color = "#FC4E07") + 
  labs(x = "Days", y = "Proportion of cumulative infections") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 2, by = 0.2)) + 
  scale_color_manual(values = c(
    "All" = "#FC4E07",
    "gr1" = "#6698FF",
    "gr2" = "#00AFBB",
    "gr3" = "#E7B800"
  ), name = "Age group") + 
  scale_fill_manual(values = c(
    "All" = "#FC4E07",
    "gr1" = "#6698FF",
    "gr2" = "#00AFBB",
    "gr3" = "#E7B800"
  ), name = "Age group") + 
  theme_classic() + 
  theme(axis.line = element_blank()) +
  geom_segment(aes(x = 0.5, xend = Tonset + 400, y = -Inf, yend = -Inf), color = "darkgray") + 
  geom_segment(x = -Inf, xend = -Inf, y = 0, yend = 1, color = "darkgray") + 
  theme(
    legend.position = c(0.05, 0.9), 
    legend.title = element_blank(), 
    legend.text = element_text(size = 11, color = "black"),
    strip.background = element_rect(color = NA), 
    strip.text = element_text(size = 13), 
    panel.grid.minor = element_line(colour = "gray", size = 0.1), 
    axis.text = element_text(size = 13, color = "black"),
    axis.ticks = element_line(color = "darkgray"), 
    axis.title = element_text(size = 13, color = "black", face = "bold")
  )
p2

outfile <- glue("figs/single/leaky/cumI_prop_R0-{R0}_ve2-{ve2}_Tonset-{Tonset}_strategy-{strategy}_nsim-{nsim}.tif")
tiff(outfile, width = 12, height = 10, unit = "in", res = 300, compression = "lzw")
print(p2)
dev.off()


# check vaccine allocation curve
# daily coverage of people administered 2 doses ####
dose2.gr1 <- dose2.gr %>% 
  mutate(dose2 = rowSums(select(., ends_with(glue("dose2.{1:16}")))), 
         dose2.gr1 = rowSums(select(., ends_with(glue("dose2.{1:3}")))), 
         dose2.gr2 = rowSums(select(., ends_with(glue("dose2.{4:11}")))), 
         dose2.gr3 = rowSums(select(., ends_with(glue("dose2.{12:16}"))))) %>% 
  group_by(sim) %>% 
  mutate(vc = cumsum(dose2) / sum(N), 
         vc.gr1 = cumsum(dose2.gr1) / sum(N.gr1), 
         vc.gr2 = cumsum(dose2.gr2) / sum(N.gr2), 
         vc.gr3 = cumsum(dose2.gr3) / sum(N.gr3)) %>% 
  ungroup() %>% 
  group_by(t) %>% 
  summarise(vc.m = median(vc), 
            vc.lwr = quantile(vc, probs = 0.025), 
            vc.upr = quantile(vc, probs = 0.975), 
            vc.gr1.m = median(vc.gr1), 
            vc.gr1.lwr = quantile(vc.gr1, probs = 0.025), 
            vc.gr1.upr = quantile(vc.gr1, probs = 0.975), 
            vc.gr2.m = median(vc.gr2), 
            vc.gr2.lwr = quantile(vc.gr2, probs = 0.025), 
            vc.gr2.upr = quantile(vc.gr2, probs = 0.975), 
            vc.gr3.m = median(vc.gr3), 
            vc.gr3.lwr = quantile(vc.gr3, probs = 0.025), 
            vc.gr3.upr = quantile(vc.gr3, probs = 0.975))

# plot vaccine coverage by age group
vc.onset <- dose2.gr1 %>% 
  select(t, vc.m, vc.gr1.m, vc.gr2.m, vc.gr3.m) %>% 
  filter(t == Tonset) %>% 
  mutate(vc.lab = format(vc.m, digits = 3, nsmall = 3), 
         vc.gr1.lab = format(vc.gr1.m, digits = 3, nsmall = 3), 
         vc.gr2.lab = format(vc.gr2.m, digits = 3, nsmall = 3), 
         vc.gr3.lab = format(vc.gr3.m, digits = 3, nsmall = 3)) %>% 
  mutate(x = Tonset, 
         gr1.x = Tonset, 
         gr2.x = Tonset, 
         gr3.x = Tonset)

fvc <- dose2.gr1 %>% 
  select(vc.m, vc.gr1.m, vc.gr2.m, vc.gr3.m) %>% 
  slice(n()) %>% 
  mutate(vc.lab = format(vc.m, digits = 3, nsmall = 3), 
         vc.gr1.lab = format(vc.gr1.m, digits = 3, nsmall = 3), 
         vc.gr2.lab = format(vc.gr2.m, digits = 3, nsmall = 3), 
         vc.gr3.lab = format(vc.gr3.m, digits = 3, nsmall = 3)) %>% 
  mutate(x = Tonset + 350, 
         gr1.x = Tonset + 350, 
         gr2.x = Tonset + 350, 
         gr3.x = Tonset + 350)

p3 <- ggplot(dose2.gr1) + 
  geom_ribbon(aes(x = t, ymin = vc.gr1.lwr, ymax = vc.gr1.upr, fill = "gr1"), 
              alpha = 0.3) + 
  geom_line(aes(t, vc.gr1.m, color = "gr1")) + 
  geom_ribbon(aes(x = t, ymin = vc.gr2.lwr, ymax = vc.gr2.upr, fill = "gr2"), 
              alpha = 0.3) + 
  geom_line(aes(t, vc.gr2.m, color = "gr2")) + 
  geom_ribbon(aes(x = t, ymin = vc.gr3.lwr, ymax = vc.gr3.upr, fill = "gr3"), 
              alpha = 0.3) + 
  geom_line(aes(t, vc.gr3.m, color = "gr3")) + 
  geom_ribbon(aes(x = t, ymin = vc.lwr, ymax = vc.upr, fill = "All"), 
              alpha = 0.3) + 
  geom_line(aes(t, vc.m, color = "All")) + 
  geom_vline(xintercept = Tonset, linetype = "dashed", color = "gray") + 
  geom_text(data = vc.onset, aes(x = gr1.x, y = vc.gr1.m + 0.05, 
                                 label = vc.gr1.lab), color = "#6698FF") + 
  geom_text(data = vc.onset, aes(x = gr2.x, y = vc.gr2.m + 0.05, 
                                 label = vc.gr2.lab), color = "#00AFBB") + 
  geom_text(data = vc.onset, aes(x = gr3.x, y = vc.gr3.m + 0.05, 
                                 label = vc.gr3.lab), color = "#E7B800") + 
  geom_text(data = vc.onset, aes(x = x, y = vc.m + 0.05, 
                                 label = vc.lab), color = "#FC4E07") + 
  geom_text(data = fvc, aes(x = gr1.x, y = vc.gr1.m + 0.05, 
                            label = vc.gr1.lab), color = "#6698FF") + 
  geom_text(data = fvc, aes(x = gr2.x, y = vc.gr2.m + 0.02, 
                            label = vc.gr2.lab), color = "#00AFBB") + 
  geom_text(data = fvc, aes(x = gr3.x, y = vc.gr3.m + 0.07, 
                            label = vc.gr3.lab), color = "#E7B800") + 
  geom_text(data = fvc, aes(x = x, y = vc.m + 0.05, 
                            label = vc.lab), color = "#FC4E07") + 
  labs(x = "Days", y = "Coverage of 2 doses administration") + 
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 2, by = 0.2)) + 
  scale_color_manual(values = c(
    "All" = "#FC4E07",
    "gr1" = "#6698FF",
    "gr2" = "#00AFBB",
    "gr3" = "#E7B800"
  ), name = "Age group") + 
  scale_fill_manual(values = c(
    "All" = "#FC4E07",
    "gr1" = "#6698FF",
    "gr2" = "#00AFBB",
    "gr3" = "#E7B800"
  ), name = "Age group") + 
  theme_classic() + 
  theme(axis.line = element_blank()) +
  geom_segment(aes(x = 0.5, xend = Tonset + 400, y = -Inf, yend = -Inf), color = "darkgray") + 
  geom_segment(x = -Inf, xend = -Inf, y = 0, yend = 1.2, color = "darkgray") + 
  theme(
    legend.position = c(0.05, 0.9), 
    legend.title = element_blank(), 
    legend.text = element_text(size = 11, color = "black"),
    strip.background = element_rect(color = NA), 
    strip.text = element_text(size = 13), 
    panel.grid.minor = element_line(colour = "gray", size = 0.1), 
    axis.text = element_text(size = 13, color = "black"),
    axis.ticks = element_line(color = "darkgray"), 
    axis.title = element_text(size = 13, color = "black", face = "bold")
  )
p3

outfile <- glue("figs/single/leaky/2dose_coverage_R0-{R0}_ve2-{ve2}_Tonset-{Tonset}_strategy-{strategy}_nsim-{nsim}.tif")
tiff(outfile, width = 12, height = 10, unit = "in", res = 300, compression = "lzw")
print(p3)
dev.off()
