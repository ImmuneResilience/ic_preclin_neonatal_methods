##########################
# Nelly Amenyogbe
# 07 May 2024
# Promethion MS: postQC analysis
##########################

# In this script, we will investigate MS figures 2-4 (2. Comparing sick & healthy pups, 3. Survival between promethion and no promethion pups, 4. Impact of chamber temperature)

# load packages
library(plyr)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(lubridate)
library(ggpubr)
library(survival)
library(survminer)

# load data
met <- read.csv("Rdata/metab_collapsed_for_analysis_pub.csv")
surv.data <- read.csv("Rdata/promethion_survival.csv")


# Pre-process data --------------------------------------------------------
# Eliminate values with a CV > 50
met <- filter(met, cv.o2 < 50, cv.co2 < 50)

# melt by metab measurement
metab.mes <- c("vo2", "vco2")
dat.m <- melt(met, measure.vars = metab.mes)

# Remove outliers for empty data
empties <- filter(dat.m, group == "empty")

# get lower and upper bounds for vo2

vo2.upper.bound <- quantile(empties$value[empties$variable == "vo2"], 0.975)
vco2.upper.bound <- quantile(empties$value[empties$variable == "vco2"], 0.975)

# set vco2 upper bound to 0.2 as well

empties$out <- ifelse(empties$value[empties$variable == "vo2"] > vo2.upper.bound, "Y", "N")
empties$out <- ifelse(empties$value[empties$variable == "vco2"] > vco2.upper.bound, "Y", empties$out)

table(empties$out, empties$variable)

id.empties <- empties$metab_identifier[empties$out == "Y"]
  
dat.m <- dat.m[-which(dat.m$metab_identifier %in% id.empties),]
  
# Define max VO2 and max VCO2 from blank cells
max.vals <- dat.m %>%
  filter(group == "empty") %>%
  ddply(.(variable),
        summarise,
        maxval = max(value))
max.vals

dat.m$blank.max[dat.m$variable == "vo2"] <- max.vals$maxval[max.vals$variable == "vo2"]
dat.m$blank.max[dat.m$variable == "vco2"] <- max.vals$maxval[max.vals$variable == "vco2"]

dat.m$above.blank <- ifelse(dat.m$value <= dat.m$blank.max, "Below blank",
                            ifelse(dat.m$value > dat.m$blank.max, "Above blank", dat.m$blank.max))


# split data by experiments
experiments <- split(dat.m, dat.m$expName)


# Analyse Temperature -----------------------------------------------------
temp <- experiments$Promethion30vsRT

# set additional variables
temp$temp.diff <- temp$temp - temp$temp2
unique(temp$group)

temp$challenge_stat <- gsub("No", "Unchallenged", temp$challenge_stat)
temp$challenge_stat <- gsub("Yes", "Challenged", temp$challenge_stat)

temp$visit.name <- factor(temp$visit.name, levels = c("baseline", "2HPC", "8HPC", "24HPC", "36HPC", "48HPC", "56HPC", "72HPC"))

# Temp differences ####
temp.p <- temp %>%
  select(temp, temp2, temp.diff, pup_ID, cage, group, hours_post_challenge, challenge_stat, visit.name, unique_pup_ID) %>%
  unique() %>%
  filter(temp.diff < 9) # There are two instances in these data where the temperature that was recorded was unusually low (20C and 30C for one pup, and 28C and 40C for another pup). These are being removed due to biological implausibility, especially since these pups did not display such different temperatures at later time points.

ggplot(temp.p, aes(x = hours_post_challenge, y = temp.diff)) +
  geom_point(shape = 21, size = 2, aes(fill = group)) +
  geom_smooth(aes(color = group)) +
  facet_wrap(~challenge_stat) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Hours post challenge", y = "Difference in temperature")

ggsave("Figures/temp/tempdiff_smooth.pdf", device = "pdf", dpi = 300, width = 6, height = 3)

# Temperature change across experiment
dat.p <- temp.p %>%
  melt(measure.vars = c("temp", "temp2"), variable.name = "temp.timepoint", value.name = "temperature")

ggplot(dat.p, aes(x = hours_post_challenge, y = temperature, fill = group)) +
  geom_smooth(aes(color = group)) +
  geom_point(shape = 21) +
  facet_grid(temp.timepoint~challenge_stat) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  labs(x = "Hours post challenge", y = "Temperature (C)")

ggsave("Figures/temp/30vsrt_temp1vs2_smooth.pdf", device = "pdf", dpi = 300, width = 6, height = 4.5)

# Plot paired boxplots for 8 HPC
dat.p$pup_visit <- paste0(dat.p$unique_pup_ID, dat.p$visit.name)
dat.p <- filter(dat.p, temperature < 38, visit.name %in% c("baseline","8HPC"))

ggplot(dat.p[dat.p$visit.name == "8HPC",], aes(x = temp.timepoint, y = temperature, fill = temp.timepoint)) +
  geom_boxplot(outlier.color = NA) +
  geom_point(shape = 21, size = 2) +
  geom_line(aes(group = unique_pup_ID)) +
  facet_grid(group ~ challenge_stat) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = "none",
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Temperature timepoint", y = "Temperature (C)") +
  geom_pwc(aes(group = temp.timepoint), method = "wilcox_test", label = "p.signif", y.position = 36) +
  ylim(c(24, 38))

ggsave("Figures/temp/temp1vs2_boxpots_wilcox.png", device = "pdf", dpi = 300, width = 3.8, height = 3)

# Temp and metab ####
ggplot(temp[temp$variable %in% c("vo2", "vco2"),], aes(x = hours_post_challenge, y = value, fill = group)) +
  geom_smooth(aes(color = group)) +
  geom_point(shape = 21) +
  facet_grid(variable ~ challenge_stat, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  ylim(c(0,0.5)) +
  labs(x = "Hours post challenge", y = "mL/min")

ggsave("Figures/temp/png/30vsrt_smooth.png", device = "png", dpi = 300, width = 6, height = 4.5)

# as boxplots  

visits <- unique(temp$visit.name)
visits <- gsub("HPC", "", visits)

temp$visit.name <- factor(temp$visit.name, levels = c("baseline", "2HPC", "8HPC", "24HPC", "36HPC", "48HPC", "56HPC", "72HPC"))

temp$visit.cat <- temp$visit.name
temp$visit.cat <- gsub("HPC", "", temp$visit.cat)
temp$visit.cat <- factor(temp$visit.cat, levels = c("baseline", "2", "8", "24", "36", "48", "56", "72"))

ggplot(temp[temp$variable %in% c("vo2", "vco2"),], aes(x = visit.cat, y = value, fill = group)) +
  geom_boxplot(outlier.color = NA) +
  geom_point(shape = 21, position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.3)) +
  facet_grid(variable ~ challenge_stat, scales = "free") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        legend.position = "bottom",
        strip.text = element_text(size = 11),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Timepoint (hours post challenge)", y = "mL/min") +
  ylim(c(0,0.6)) +
  geom_pwc(aes(group = group), method = "wilcox_test", label = "p.signif", y.position = 0.55,)

#ggsave("Figures/temp/30vsrt_vo2_boxplots.png", device = "png", dpi = 300, width = 8, height = 6)

# This shows us that the vh2o differed quite a bit between conditions.  Among challenged animals, the general kinetics following challenge are still observed but the response rate shows more heterogeneity at 30C compared to RT for both VO2 and VCO2.

# Investigate replicate values
# For the older mice, are the wider range of Vo2 values reflected by a larger drop from the first to the 2nd measurement?  We can use the raw, unaveraged data to find out.

dat.reps <- read.csv("Rdata/metab_replicates_fot_analysis.csv")
dat.reps$cv.o2.abs <- abs(dat.reps$cv.o2)
dat.reps$cv.co2.abs <- abs(dat.reps$cv.co2)
dat.reps <- filter(dat.reps, cv.co2.abs < 50, cv.o2.abs < 50)

unique(dat.reps$expName)

rep.temp <- filter(dat.reps, expName == "Promethion warming 30C vs RT")
# set aesthetics
rep.temp$visit.cat <- rep.temp$visit.name
rep.temp$visit.cat <- gsub("HPC", "", rep.temp$visit.cat)
rep.temp$visit.cat <- factor(rep.temp$visit.cat, levels = c("baseline", "2", "8", "24", "36", "48", "56", "72"))

# bin the VO2
rep.temp$vo2.cv <- cut(rep.temp$cv.o2, breaks = c(0,25,50,75,100,200))
rep.temp$vco2.cv <- cut(rep.temp$cv.co2, breaks = c(0,25,50,75,100,200))

# label each identifier as first reading higher or lower
rep.temp <- ddply(rep.temp,
                  .(metab_identifier),
                  transform,
                  vo2.diff = vo2[reading == 1] - vo2[reading == 2],
                  vco2.diff = vco2[reading == 1] - vco2[reading == 2])

rep.temp$vo2.change <- ifelse(rep.temp$vo2.diff >=0, "decrease", "increase")
rep.temp$vco2.change <- ifelse(rep.temp$vco2.diff >= 0, "decrease", "increase")

# plot increase or decrease
# VO2
ggplot(rep.temp[rep.temp$challenge_stat == "No",], aes(x = as.factor(reading), y = vo2, group = metab_identifier, color = vo2.change)) +
  geom_point() +
  geom_line() +
  facet_grid(group~visit.cat) +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Reading", y = "VO2 (mL/min")

ggsave("Figures/temp/reading1vs2_changetype.png", device = "png", dpi = 300, width = 6, height = 4)

# VCO2
ggplot(rep.temp[rep.temp$challenge_stat == "No",], aes(x = as.factor(reading), y = vco2, group = metab_identifier, color = vco2.change)) +
  geom_point() +
  geom_line() +
  facet_grid(group~visit.cat) +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Reading", y = "VCO2 (mL/min")

ggsave("../../Figures/temp/reading1vs2_vco2_changetype.png", device = "png", dpi = 300, width = 6, height = 4)


# plot CV
#VO2
ggplot(rep.temp[rep.temp$challenge_stat == "No",], aes(x = as.factor(reading), y = vo2, group = metab_identifier, color = vo2.cv)) +
  geom_point() +
  geom_line() +
  facet_grid(group~visit.cat) +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Reading", y = "VO2 (mL/min")

#ggsave("Figures/temp/reading1vs2_cvo2.png", device = "png", dpi = 300, width = 6, height = 4)

#VCO2
ggplot(rep.temp[rep.temp$challenge_stat == "No",], aes(x = as.factor(reading), y = vco2, group = metab_identifier, color = vco2.cv)) +
  geom_point() +
  geom_line() +
  facet_grid(group~visit.cat) +
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Reading", y = "VCO2 (mL/min")

#ggsave("Figures/temp/reading1vs2_cv_co2.png", device = "png", dpi = 300, width = 6, height = 4)

temp.samplesize <- temp %>% 
  select(unique_pup_ID, group, challenge_stat) %>%
  unique()

table(temp.samplesize$group, temp.samplesize$challenge_stat)


# Analyze age -------------------------------------------------------------
age <- experiments$Promethion_age

# calculate pup age
age$dob_pups <- dmy(age$DOB_pup.dd.mm.yyyy)

age$date_dd.mm.yyyy <-  gsub("-", "_", age$date_dd.mm.yyyy)
age$date_dd.mm.yyyy <- ymd(age$date_dd.mm.yyyy)

# calculate DOL
age$dol <- difftime(age$date_dd.mm.yyyy,age$dob_pups, units = 'day') %>%
  as.numeric()
# create factor
age$dol.f <- age$dol %>%
  as.character() 

age$dol.f <- paste0("D", age$dol.f)

# Age metab ####

# add empty cages
"%notin%" <- Negate("%in%")

empty <- experiments$Empty.Tubes
empty$dol <- 0

age.blank <- rbind.fill(age, empty)

# plot age line graph ####

ggplot(age.blank, aes(x = dol, y = value, fill = group)) +
  geom_smooth(color = "#636363") +
  geom_point(shape = 21) +
  facet_wrap(~variable, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  labs(x = "Age (days)", y = "mL/min")

#ggsave("Figures/age/png/vo2_age_linegraph.png", device = "png", dpi = 300, width = 5.5, height = 3.5)


# Age boxplots ####

# calculate number below blank.
age.blank$dol.f <- ifelse(age.blank$group == "empty", "empty", age.blank$dol.f)
age.blank$dol.f <- factor(age.blank$dol.f, levels = c("empty", "D1", "D3", "D5", "D7", "D8", "D11", "D14"))

age.p <- age.blank %>%
  select(unique_pup_ID, metab_identifier, dol.f, dol, variable, value, group, above.blank) %>%
  filter(dol < 8, variable %in% c("vo2", "vco2"))

age.p$above.blank <- factor(age.p$above.blank, levels = c("Below blank", "Above blank"))

# plot boxplot
ggplot(age.p, aes(x = dol.f, y = value)) +
  geom_boxplot(outlier.color = NA, fill = "#bdbdbd") +
  geom_point(shape = 21, size = 2, position = position_jitter(width = 0.3), aes(fill = above.blank)) +
  facet_grid(~variable) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = "bottom", 
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "", y = "mL/min")

ggsave("Figures/age/age_boxplot.png", device = "png", dpi = 300, width = 5, height = 4)

# get summary table

age.p$variable <- droplevels(age.p$variable)
pct <- ddply(age.p,
             .(variable, dol.f),
             summarise,
             n = length(metab_identifier),
             n_above = length(which(above.blank == "Above blank")),
             above.blank = 100 *length(which(above.blank == "Above blank")) / length(metab_identifier))
pct$below.blank = 100-pct$above.blank

write.csv(pct, "Tables/age_pctabove.csv", row.names = FALSE)

# plot threshold
pct.m <- melt(pct, measure.vars = c("above.blank", "below.blank"), variable.name = "threshold", value.name = "threshold.percentage")

pct.m$threshold <- ifelse(pct.m$threshold == "above.blank", "Above blank", "Below blank")
pct.m$threshold <- factor(pct.m$threshold, levels = c("Below blank", "Above blank"))

ggplot(pct.m[pct.m$dol.f != "empty",], aes(x = dol.f, y = threshold.percentage, fill = threshold)) +
  geom_bar(color = "black", stat = "identity") +
  facet_wrap(~variable) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Day of life", y = "Percent of observations")

ggsave("Figures/age/age_aboveblank_percentagebars.png", device = "png", dpi = 300, width = 4.5, height = 4)


# plot histograms
age.p$dol.f <- factor(age.p$dol.f, levels = c("empty", "D1", "D3", "D5", "D7"))

ggplot(age.p, aes(x = value, fill = dol.f)) +
  geom_density(alpha = 0.7) +
  geom_histogram(color = "black", alpha = 0.3, binwidth = 0.01) +
  facet_grid(~variable) +
  theme_bw() +
  scale_fill_manual(values = c("#737373", '#e66101','#fdb863','#b2abd2','#5e3c99')) +
  theme(axis.text = element_text(size = 12),
        legend.position = "bottom", 
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  xlab("mL/min")

#ggsave("Figures/age/png/age_vo2_histogram.png", device = "png", dpi = 300, width = 6, height = 3.5)

age.samplesize <- age %>%
  select(cage, pup_ID, group, dol.f) %>%
  unique()

table(age.samplesize$group, age.samplesize$dol.f)

# Analyse Promethion survival -------------------------------------------------------
df.surv <- surv.data %>%
  select(pup_ID, group, cage, sex, Exclude, Exclude_Notes, expName, outcome, sac_hpc) %>%
  unique()

table(df.surv$group)

# perform log-rank test May ####
surv <- with(df.surv, Surv(sac_hpc, outcome))

fit <- survfit(surv ~ group,
               data = df.surv,
               conf.int = 0.95,
               conf.type = "log")

test <- survdiff(surv ~ group, data  = df.surv, rho = 0)

pval.all <- broom::glance(test)$p.value
pval.all # p = 0.381

# plot kaplan-meier
km <- ggsurvplot(fit,
                 conf.int = FALSE,
                 fun = "pct",
                 linetype = c(2,1),
                 size = 1.0,
                 censor = FALSE,
                 legend = "bottom",
                 xlab = "Hours Post Challenge",
                 legend.labs = names(table(df.surv$group)),
                 palette = c("#969696", "#ef6548"))



# annotoate plot
anot <- paste0("p = ", round(pval.all, digits = 3))
anot

km$plot + annotate("text", size = 5, x = 10, y = 10, label = anot) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_blank())

ggsave("Figures/PPvsSS/PPvsSS_surv.pdf", device = "pdf", dpi = 300, width = 4, height = 4)

# Weight differences ####
surv.data <- ddply(surv.data,
                .(pup_ID, cage),
                transform,
                pct.weight = 100*(weight - weight[vis.type == "challenge"])/weight[vis.type == "challenge"])


ggplot(surv.data, aes(x = hours_post_challenge, y = pct.weight, fill = group)) +
  geom_smooth(aes(color = group)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = c("#969696", "#ef6548","blue")) +
  scale_color_manual(values = c("#969696", "#ef6548", "blue")) +
  theme_classic() +
  labs(y = "Percent Weight Change", x = "Hours post challenge") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title = element_text(size = 12))

ggsave("Figures/PPvsSS/PPvsSSweight_smooth.pdf", device = "pdf", dpi = 300, width = 4, height = 4)

# plot scores
ggplot(surv.data, aes(x = hours_post_challenge, y = avg_score, fill = group)) +
  geom_smooth(aes(color = group)) +
  scale_fill_manual(values = c("#969696", "#ef6548")) +
  scale_color_manual(values = c("#969696", "#ef6548")) +
  theme_classic() +
  labs(y = "Cliinical Score", x = "Hours post challenge") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title = element_text(size = 12))

# what was the visit schedule?
visits <- c(1:16)
visit <- paste0("v", visits)

surv.data$visit <- factor(surv.data$visit, levels = visit)
surv.data$Promethion <- ifelse(surv.data$group == "Promethion" & visit %in% c("v1", "v3", "v4", "v5", "v6", "v7", "v10", "v15"), "Yes", "No")


ggplot(surv.data, aes(x = visit, y = hours_post_challenge, fill = Promethion, shape = group)) +
  geom_point(size = 4) +
  scale_shape_manual(values = c("No promethion" = 21,"Promethion" =  22)) +
  theme_bw() +
  labs(x = "Visit", y = "Hours post challenge") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_fill_manual(values  = c("No" =  "#969696", "Yes" =  "#ef6548"))

pp <- filter(surv.data, group == "No promethion")
visit.names <- ddply(pp, 
                     .(visit),
                     summarise,
                     hpc = mean(hours_post_challenge))
visit.names


# get promethion time points
pr.surv <- experiments$Promethion_vs_NoPromethion

ggplot(pr.surv, aes(x = visit, y = hours_post_challenge)) +
  geom_point()

# For baseline, V1 for both
# for 12 HPC, V4 for No Promethion, V8 for Promethion
# for 20 HPC, V5 for NP, V9 for P
# for 48 HPC, V11 for NP, V15 for P

visit.names <- data.frame(group = c("No promethion", "No promethion", "No promethion", "No promethion", "Promethion", "Promethion", "Promethion", "Promethion"),
                          visit = c("v1", "v4", "v5", "v11", "v1", "v8", "v9", "v15"),
                          timepoint = c("0 HPC", "12 HPC", "20 HPC", "48 HPC"))

surv.data2 <- left_join(surv.data, visit.names, by = c("group" = "group", "visit" = "visit"))
surv.data2 <- surv.data2[-which(is.na(surv.data2$timepoint)),]

# plot score
surv.data2 <- surv.data2[-which(is.na(surv.data2$avg_score)),]
surv.data2$score.category <- as.character(surv.data2$avg_score)

# get proportions
ps.ag <- as.data.frame(table(surv.data2$score.category, surv.data2$timepoint, surv.data2$group))
ps.ag <- filter(ps.ag, Freq > 0)

colnames(ps.ag) <- c("score", "timepoint", "group", "Freq")

ggplot(ps.ag, aes(x = group, y = Freq, fill = score)) +
  geom_bar(color = "black", stat = "identity") +
  facet_grid(~timepoint) +
  scale_fill_manual(values = c('#b2182b','#ef8a62','#fddbc7','#d1e5f0','#67a9cf','#2166ac')) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1.0, vjust = 1.0),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "", y = "Number of observations")

ggsave("Figures/PPvsSS/ppvsss_scores_bargraph.png", device = "png", dpi = 300, width = 5, height = 4)

# Analyse sick vs. healthy ---------------------------------------------

# for this analysis, let's use the temp pilot and empties.

# prepare data
sh <- filter(dat.m, group %in% c("PP") | expName == "Promethion30vsRT") 

# subset to empties and 2 HPC only
table(sh$visit.name, sh$group) # all the empties have visit name as monitor. The PP mice do not have a visit name, but we can grab all the values under 3 to include the 1-2 hour mice.

sh.p <- filter(sh, visit.name == "monitor" | 
                 hours_post_challenge < 3 & hours_post_challenge > 1)

# add empty tubes
empty.add <- filter(dat.m, group == "empty")

sh.p <- rbind(sh.p, empty.add)

sh.p$challenge_stat <- ifelse(sh.p$group == "empty", "empty", sh.p$challenge_stat)

sh.p$group.p <- gsub("Yes", "sepsis", sh.p$challenge_stat)
sh.p$group.p <- gsub("No", "healthy", sh.p$group.p)

samplesize <- sh.p %>%
  select(metab_identifier, group.p) %>%
  unique()

table(samplesize$group.p)

# plot boxplot
sh.p$above.blank <- factor(sh.p$above.blank, levels  = c("Below blank", "Above blank"))

ggplot(sh.p, aes(x = group.p, y = value)) +
  geom_boxplot(outlier.color = NA, fill = "#bdbdbd") +
  geom_point(shape = 21, size = 2, position = position_jitter(width = 0.3), aes(fill = above.blank)) +
  facet_grid(~variable) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = "bottom", 
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "", y = "mL/min")

table(sh.p$group.p, sh.p$variable)
ggsave("Figures/sepsis/sickvshealthy_boxplot.png", device = "png", dpi = 300, width = 5, height = 4)

# get summary table

sh.p$variable <- droplevels(sh.p$variable)
pct <- ddply(sh.p,
             .(variable, group.p),
             summarise,
             n = length(metab_identifier),
             n_above = length(which(above.blank == "Above blank")),
             above.blank = 100 *length(which(above.blank == "Above blank")) / length(metab_identifier))
pct$below.blank = 100-pct$above.blank

write.csv(pct, "Tables/sickvshealthy_pctabove.csv", row.names = FALSE)

# plot threshold
pct.m <- melt(pct, measure.vars = c("above.blank", "below.blank"), variable.name = "threshold", value.name = "threshold.percentage")
pct.m$threshold <- ifelse(pct.m$threshold == "above.blank", "Above blank", "Below blank")
pct.m$threshold <- factor(pct.m$threshold, levels = c("Below blank", "Above blank"))

ggplot(pct.m[pct.m$group.p != "empty",], aes(x = group.p, y = threshold.percentage, fill = threshold)) +
  geom_bar(color = "black", stat = "identity") +
  facet_wrap(~variable) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "", y = "Percent of observations")

ggsave("Figures/sepsis/sickvshealthy_percentagebars.pdf", device = "pdf", dpi = 300, width = 3.8, height = 4)

# plot histograms
sh.p$group.p <- factor(sh.p$group.p, levels = c("empty", "sepsis", "healthy"))

ggplot(sh.p[sh.p$variable %in% c("vo2", "vco2"),], aes(x = value, fill = group.p)) +
  geom_density(alpha = 0.7) +
  geom_histogram(color = "black", alpha = 0.3, binwidth = 0.01) +
  facet_grid(~variable) +
  theme_bw() +
  scale_fill_manual(values = c("#737373", "#d6604d", "#4393c3")) +
  theme(axis.text = element_text(size = 12),
        legend.position = "bottom", 
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  xlab("mL/min")

ggsave("Figures/sepsis/sickvshealthy_histogram.pdf", device = "pdf", dpi = 300, width = 6, height = 3.5)

# END ####