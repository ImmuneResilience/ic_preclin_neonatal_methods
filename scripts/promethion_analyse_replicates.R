##########################
# Adrien Eynaud
# 09 May 2026
# Promethion MS: analysis of technical replicates
##########################

# In this script, we will investigate MS Figure 1 and Supplementary Figure 1
# (technical replicates & coefficient of variation).


# load packages
library(plyr)
library(tidyverse)
library(lubridate)

# load data
dat.pub <- read.csv("Rdata/metab_replicates_for_analysis.csv")


# Figure 1 plots ----------------------------------------------------------

# Stratify by CV across the technical replicate readings.
high.cvs <- filter(dat.pub, cv.o2 >= 0 & cv.co2 >= 0)
high.cvs <- high.cvs[, which(colnames(high.cvs) %in% c("metab_identifier",
                                                       "reading",
                                                       "cv.o2",
                                                       "cv.co2",
                                                       "vo2",
                                                       "vco2"))]

long.high.cvs <- pivot_longer(
  high.cvs,
  cols = c(vo2, vco2),
  names_to = c("type"),
  values_to = c("value")
)

long.high.cvs <- pivot_longer(
  long.high.cvs,
  cols = c(cv.o2, cv.co2),
  values_to = c("cv_value")
)

ready.frame <- NULL
for (i in unique(long.high.cvs$metab_identifier)) {
  sub <- filter(long.high.cvs, metab_identifier == i)
  sub <- sub[c(5, 1, 8, 4), ]
  ready.frame <- rbind(ready.frame, sub)
}

ready.frame$cv_category <- NA
ready.frame$cv_category <- ifelse(ready.frame$cv_value < 25, "CV less than 25%", ready.frame$cv_category)
ready.frame$cv_category <- ifelse(ready.frame$cv_value >= 25 & ready.frame$cv_value < 50, "CV 25% to 50%", ready.frame$cv_category)
ready.frame$cv_category <- ifelse(ready.frame$cv_value >= 50, "CV 50%+", ready.frame$cv_category)

ready.frame$cv_category <- factor(ready.frame$cv_category,
                                  levels = c("CV less than 25%", "CV 25% to 50%", "CV 50%+"))
ready.frame$type <- factor(ready.frame$type, levels = c("vo2", "vco2"))

torm <- ready.frame[which(ready.frame$value > 0.75 | ready.frame$value < 0), ]$metab_identifier
`%notin%` <- Negate(`%in%`)
ready.frame <- ready.frame[which(ready.frame$metab_identifier %notin% torm), ]

# Figure 1B
p.reps <- ggplot(ready.frame, aes(x = as.factor(reading), y = value, group = metab_identifier)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_line(alpha = 0.5, linewidth = 0.15, na.rm = TRUE) +
  theme_bw() +
  facet_grid(type ~ cv_category, scales = "fixed", space = "free_x", switch = "y") +
  labs(x = "Reading", y = "Metabolic Output", color = "Metabolic Identifier") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.text.x = element_text(margin = margin(t = 10, b = 10)),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p.reps
ggsave("Figures/cv_replicates/cv_replicate_readings.png", plot = p.reps, device = "png", dpi = 300, width = 5.2, height = 4)

# Figure 1C
pct.cv <- data.frame(table(ready.frame$cv_category, ready.frame$name))
pct.n <- ddply(pct.cv,
               .(Var2),
               summarise,
               n = sum(Freq))

pct.cv$n <- pct.n$n[match(pct.cv$Var2, pct.n$Var2)]
pct.cv$percent <- pct.cv$Freq / pct.cv$n * 100
pct.cv$Var2 <- ifelse(pct.cv$Var2 == "cv.o2", "VO2", "VCO2")

p.pct <- ggplot(pct.cv, aes(x = Var2, y = percent, fill = Var1)) +
  geom_bar(color = "black", stat = "identity") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title = element_text(size = 12)) +
  labs(x = "Measurement", y = "Percent of observations")

p.pct
ggsave("Figures/cv_replicates/percent_of_observations.png", plot = p.pct, device = "png", dpi = 300, width = 2.5, height = 4)


# Supplementary Figure 1 --------------------------------------------------

# Fig. S1A-B: reading 1 versus reading 2 directionality
ready.frame2 <- ddply(ready.frame,
                      .(metab_identifier, name),
                      transform,
                      diff = value[reading == 1] - value[reading == 2])

ready.frame2$direction <- ifelse(ready.frame2$diff >= 0, "R1", "R2")

p.reps.r1r2 <- ggplot(ready.frame2, aes(x = as.factor(reading), y = value, group = metab_identifier, color = direction)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_line(alpha = 0.5, linewidth = 0.15, na.rm = TRUE) +
  theme_bw() +
  facet_grid(type ~ cv_category, scales = "fixed", space = "free_x", switch = "y") +
  labs(x = "Reading", y = "Metabolic Output", color = "Highest Reading") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 12, margin = margin(t = 10, b = 10)),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
p.reps.r1r2
ggsave("Figures/cv_replicates/cv_replicate_readings_R1R2.png", plot = p.reps.r1r2, device = "png", dpi = 300, width = 5.2, height = 4)

dat.pct <- ready.frame2 %>%
  select(cv_category, name, metab_identifier, direction) %>%
  unique()

pct.cv.r <- data.frame(table(dat.pct$cv_category, dat.pct$name, dat.pct$direction))
pct.n.r <- ddply(pct.cv.r,
                 .(Var2),
                 summarise,
                 n = sum(Freq))

pct.cv.r$n <- pct.n.r$n[match(pct.cv.r$Var2, pct.n.r$Var2)]
pct.cv.r$percent <- pct.cv.r$Freq / pct.cv.r$n * 100
pct.cv.r$Var2 <- ifelse(pct.cv.r$Var2 == "cv.o2", "VO2", "VCO2")

p.pct.r1r2 <- ggplot(pct.cv.r, aes(x = Var2, y = percent, fill = Var3)) +
  geom_bar(color = "black", stat = "identity") +
  facet_wrap(~Var1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Measurement", y = "Percent of observations")

p.pct.r1r2
ggsave("Figures/cv_replicates/percent_of_observations_R1R2.png", plot = p.pct.r1r2, device = "png", dpi = 300, width = 4.5, height = 3.5)

# Fig. S1C: proportion of high CV values across years
collapsed <- dat.pub[which(dat.pub$reading == 1), ]
full.tubes <- filter(collapsed, group == "pup", vo2 > 0, cv.o2 > 0)

full.tubes$year <- year(ymd(full.tubes$expDate))
full.tubes$cv_high <- ifelse(full.tubes$cv.o2 > 50, TRUE, FALSE)
full.tubes$year_group <- ifelse(full.tubes$year == 2021, "2021", "2022-2023")

contingency.table <- table(full.tubes$year_group, full.tubes$cv_high)
fisher.test(contingency.table)

pct.cv.year <- data.frame(table(full.tubes$year_group, full.tubes$cv_high))
pct.n.year <- ddply(pct.cv.year,
                    .(Var1),
                    summarise,
                    n = sum(Freq))

pct.cv.year$n <- pct.n.year$n[match(pct.cv.year$Var1, pct.n.year$Var1)]
pct.cv.year$percent <- pct.cv.year$Freq / pct.cv.year$n * 100
pct.cv.year$Var2 <- ifelse(pct.cv.year$Var2 == FALSE, "CV less than 50%", "CV 50%+")
pct.cv.year$Var2 <- factor(pct.cv.year$Var2, levels = c("CV less than 50%", "CV 50%+"))

p.year <- ggplot(pct.cv.year, aes(x = Var1, y = percent, fill = Var2)) +
  geom_bar(color = "black", stat = "identity") +
  scale_fill_manual(values = c("CV less than 50%" = "#00BFC4", "CV 50%+" = "#F8766D")) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title = element_text(size = 12)) +
  labs(x = "Experimental Year", y = "Percent of observations")

p.year

ggsave("Figures/cv_replicates/percent_of_obs_2021_vs_2022_2023.png", plot = p.year, device = "png", dpi = 300, width = 2.5, height = 4)
