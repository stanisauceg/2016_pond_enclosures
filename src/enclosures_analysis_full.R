# FIRST, load functions from "functions.R"
# SECOND, libraries below:

library(data.table)
library(lubridate)
library(car)
library(tidyverse)
library(vegan)
library(SimComp)
library(lattice)
library(corrplot)
library(lme4)
library(pbkrtest)
library(grid)
library(gridExtra)


path <- file.path("data", "data-keyed.csv")
data_keyed <- fread(path,header=TRUE)

# set column types
data_keyed <- data_keyed %>%
  as.tbl() %>%
  mutate(TimeCat = factor(Time, labels = c("May 5 2016", "May 15 2016", "May 27 2016", "June 9 2016")),
         TimeCat = mdy(TimeCat),
         block = as.factor(block),
         treatment = as.factor(treatment),
         sal = as.factor(sal),
         man = as.factor(man),
         ctrl = as.factor(ctrl)) %>%
  arrange(TimeCat, block)

# count days since salamander introduction
data_keyed$days <- data_keyed$Time
data_keyed$days[which(data_keyed$Time==1)]<-5
data_keyed$days[which(data_keyed$Time==2)]<-15
data_keyed$days[which(data_keyed$Time==3)]<-27
data_keyed$days[which(data_keyed$Time==4)]<-40

# path <- file.path("data", "mesozoop1234.csv")
# rawcounts <- fread(path, header=TRUE)
# 
# rawcounts <- drop_na(rawcounts)
# 
# countsmelt <- melt(rawcounts,measure.vars = 4,id.vars = 1:3)
# countscast <- dcast(countsmelt, Date + Block + Treatment ~ value)
# countscast$L <- rowSums(countscast[,4]!=0)
# countscast$N <- rowSums(countscast[,5]!=0)
# countscast$R <- rowSums(countscast[,6]!=0)

# load water_depths
path <- file.path("data", "water_depths.csv")
water_depths <- fread(path, header=TRUE)

colnames(water_depths) <- tolower(colnames(water_depths))

# make block a factor
water_depths$block <- as.factor(water_depths$block)

# extract initial depths from water_depths, append to data_keyed
depth_init<-unlist(c(subset(water_depths, select = "depth_may5"),
                      subset(water_depths, select = "depth_may15-pre"),
                      subset(water_depths, select = "depth_may27-pre"),
                      subset(water_depths, select = "depth_jun9-pre")))
data_keyed$depth_init <- depth_init

rm(water_depths, depth_init)

# use depths to calculate volumes of water in enclosures
data_keyed$vol_pre <- volume(data_keyed$depth_init)

# back-check new.depth fxn with example values: do I get the correct sample volume left over?
depth <- 10
sample_volume <- 500
old_vol <- volume(d = depth)
new_vol <- new.depth(d.start = depth, v.sample = sample_volume) %>% numeric.depth() %>% volume()
sample_est <- old_vol - new_vol
sample_est # compare to sample_volume = 500
# good, within 0.1%
rm(depth, sample_volume, old_vol, new_vol, sample_est)

# calculate sample volumes for all samples - initially INCORRECT for manual treatments
data_keyed$sample_vol_flush <- lapply(data_keyed$depth_init, FUN = sample.vol.flush)
data_keyed$sample_vol_hi <- lapply(data_keyed$depth_init, FUN = sample.vol.hi)

# identify the manual treatments
man_buckets <- which(data_keyed$treatment == "Manual")
# calculate correctly, 6 dips per sample
data_keyed$sample_vol_flush[man_buckets] <- lapply(X = data_keyed$depth_init[man_buckets], 
                                                   FUN = sample.vol.flush, ndips = 6)
data_keyed$sample_vol_hi[man_buckets] <- lapply(X = data_keyed$depth_init[man_buckets], 
                                                FUN = sample.vol.hi, ndips = 6)

rm(man_buckets)

data_keyed$sample_vol_flush <- unlist(data_keyed$sample_vol_flush)
data_keyed$sample_vol_hi <- unlist(data_keyed$sample_vol_hi)

data_keyed <- data_keyed %>% 
  mutate(sample_frac = sample_vol_flush / vol_pre,
         vol_post = vol_pre - sample_vol_flush)

depths <- numeric(length = 120)

# calculate enclosure depth following each sampling run
for (i in 1:120) {
  depths[i] <- new.depth(data_keyed$depth_init[i], data_keyed$sample_vol_flush[i]) %>% numeric.depth()
}

data_keyed$depth_post <- depths
rm(depths, path)

data_zoop <- data_keyed %>%
  select(-`adult coleopteran`,
         -Hebridae,
         -Isotomidae,
         -`larval coleopteran 1`,
         -`larval coleopteran 2`,
         -`larval dipteran`,
         -Sminthuridae,
         -`unID Hydracarina`,
         -`unID Nematoda`,
         -`unID Thysanoptera`,
         -`unID Turbellaria`,
         -`unID Tardigrada`,
         -`unknown Heteroptera 1`,
         -`unknown Heteroptera 2`,
         -`unknown Heteroptera 3`,
         -`unknown Heteroptera 4`,
         -`unknown Heteroptera 5`) %>%
  rename(Candona_sp_1 = 'Candona sp 1',
         Candona_sp_2 = 'Candona sp 2',
         Pleuroxus.denticulatus = 'Pleuroxus denticulatus',
         unknown_rotifer_1 = 'unknown Rotifer 1')

names(data_zoop)

# recalculate indices
data_zoop$nTot <- rowSums(data_zoop[,2:18]) # columns 2:18 are individual taxon counts after restriction
data_zoop$richness <- rowSums(data_zoop[,2:18] != 0)
data_zoop$InvSimpson <- 1 / rowSums((data_zoop[,2:18] / data_zoop$nTot) ^2)
data_zoop$evenness <- data_zoop$InvSimpson / data_zoop$richness

data_zoop$z_density <- data_zoop$nTot / data_zoop$sample_vol_flush

data_keyed$z_density <- data_keyed$nTot / data_keyed$sample_vol_flush

# restrict data to stocked spp only
data_seeded <- data_zoop %>%
  select(-Asplanchna,
         -Candona_sp_1,
         -Candona_sp_2,
         -Euchlanis,
         -Keratella,
         -Lecane,
         -Lepadella,
         -Pleuroxus.denticulatus,
         -Scapholeberis,
         -Trichocerca,
         -unknown_rotifer_1)
head(data_seeded)

# recalculate indices
data_seeded$nTot <- rowSums(data_seeded[,2:7]) # columns 2:7 are individual taxon counts after restriction
data_seeded$richness <- rowSums(data_seeded[,2:7] != 0)
data_seeded$InvSimpson <- 1 / rowSums((data_seeded[,2:7] / data_seeded$nTot)^2)
data_seeded <- data_seeded %>% mutate(evenness = InvSimpson / richness,
                                      z_density = nTot / sample_vol_flush)


### Rarefaction ----

data_zoop %>%
  group_by(TimeCat) %>%
  summarise(min = min(nTot),
            mean = mean(nTot),
            median = median(nTot),
            max = max(nTot))

# extract just the counts of animals
community_keyed <- data_keyed[,2:35]
community_zoop <- data_zoop[,2:18]
community_stock <- data_seeded[,2:7]

rarecurve(community_zoop, step = 20, sample = 100, col = "blue")

# rarefy richness to the minimum # of animals in each time step
rarefy1 <- rarefy(community_keyed, sample = min(data_keyed[data_keyed$Time == 1, 38]))
rarefy2 <- rarefy(community_keyed, sample = min(data_keyed[data_keyed$Time == 2, 38]))
rarefy3 <- rarefy(community_keyed, sample = min(data_keyed[data_keyed$Time == 3, 38]))
rarefy4 <- rarefy(community_keyed, sample = min(data_keyed[data_keyed$Time == 4, 38]))

# extract rarefied richness from each time step, combine into single vector
# using single lowest sample size
rar.S33 <- rarefy1
data_keyed$rarefied_S33 <- rar.S33

# using lowest sample size for each time step
rar.S <- unlist(c(rarefy1[1:30], rarefy2[31:60], rarefy3[61:90], rarefy4[91:120])) 
data_keyed$rarefied_S <- rar.S

# rarefy richness to the minimum # of animals in each time step
rarefy1 <- rarefy(community_zoop, sample = min(data_zoop[data_zoop$Time == 1, 21]))
rarefy2 <- rarefy(community_zoop, sample = min(data_zoop[data_zoop$Time == 2, 21]))
rarefy3 <- rarefy(community_zoop, sample = min(data_zoop[data_zoop$Time == 3, 21]))
rarefy4 <- rarefy(community_zoop, sample = min(data_zoop[data_zoop$Time == 4, 21]))

# extract rarefied richness from each time step, combine into single vector
# using single lowest sample size
rar.S33 <- rarefy1
data_zoop$rarefied_S33 <- rar.S33

# using lowest sample size for each time step
rar.S <- unlist(c(rarefy1[1:30], rarefy2[31:60], rarefy3[61:90], rarefy4[91:120])) 
data_zoop$rarefied_S <- rar.S

# rarefy richness to the minimum # of animals in each time step
rarefy1 <- rarefy(community_stock, sample = min(data_seeded[data_seeded$Time == 1, 10]))
rarefy2 <- rarefy(community_stock, sample = min(data_seeded[data_seeded$Time == 2, 10]))
rarefy3 <- rarefy(community_stock, sample = min(data_seeded[data_seeded$Time == 3, 10]))
rarefy4 <- rarefy(community_stock, sample = min(data_seeded[data_seeded$Time == 4, 10]))

# extract rarefied richness from each time step, combine into single vector
rar.S <- unlist(c(rarefy1[1:30], rarefy2[31:60], rarefy3[61:90], rarefy4[91:120])) 
data_seeded$rarefied_S <- rar.S

rm(rar.S, rar.S33, rarefy1, rarefy2, rarefy3, rarefy4)
rm(community_keyed, community_zoop, community_stock)

### CV Analysis ----

# create data frame of CV of each enclosure across the 4 time steps
cvs.enclosure <- data_zoop %>%
  select(Asplanchna:block) %>%
  group_by(treatment, block) %>%
  summarise_at(vars(Asplanchna:unknown_rotifer_1), ~ sd(.)/mean(.))

cvs.enclosure

# calculate treatment means of CVs
cvs.treatment <- cvs.enclosure %>% 
  group_by(treatment) %>%
  select(-block) %>%
  summarise_all(~mean(na.omit(.)))

cvs.treatment

# keep only species w/o NaN values
cvs.whole <- cvs.treatment %>%
  select_if(~ all(!is.nan(.)))

cvs.whole

# stocked species only
cvs.stocked <- cvs.treatment %>% # was cvs3.stocked
  select(treatment, Chydorus, Simocephalus, Cyclopoid, Harpacticoid, Cypricercus, Cyclocypris)

cvs.stocked

## Visualize CVs

# transform cvs.stocked for plotting
cvs.stocked.long <- cvs.stocked %>% gather(species, CV, Chydorus:Cyclocypris)

ggplot(cvs.stocked.long, aes(x = treatment, y = CV, color = species, group = species)) +
  geom_point() +
  geom_path() +
  ylim(c(0, 2.5))

# what are the error bars on those mean CV's?
cvs.enclosure.long <- cvs.enclosure %>%
  select(treatment, block, Chydorus, Simocephalus, Cyclopoid, Harpacticoid, Cypricercus, Cyclocypris) %>%
  gather(species, CV, Chydorus:Cyclocypris)

dodge = position_dodge(width = 0.1)

ggplot(cvs.enclosure.long, aes(x = treatment, y = CV, color = species, group = species)) +
  stat_summary(fun.data = mean_se, geom = "pointrange", position = dodge) +
  stat_summary(fun.y = mean,geom="line",aes(group=species),linetype=2,position= dodge) +
  scale_x_discrete(labels = c("Control", "Thinning", "Predation")) +
  ylim(c(0, 2.2))


# get CVs of diversity indices
cvs.diversity <- data_zoop %>%
  select(treatment:InvSimpson, evenness:rarefied_S) %>%
  group_by(treatment, block) %>%
  summarise_at(vars(richness:rarefied_S), ~ sd(.)/mean(.))

cvs.diversity

# test for differences in CVs
m1<-lm(cbind(Chydorus,Cyclocypris,Cyclopoid,Cypricercus,Harpacticoid,Simocephalus)~treatment,data=cvs.enclosure)
car::Manova(m1)

m2<-lm(cbind(Chydorus,Cyclocypris,Cyclopoid,Cypricercus,Harpacticoid,Simocephalus,
             Keratella)~treatment,data=cvs.enclosure)
car::Manova(m2)

m3<-aov(z_density~treatment,data=cvs.diversity)
m4<-aov(richness~treatment,data=cvs.diversity)
m5<-aov(InvSimpson~treatment,data=cvs.diversity)
m6<-aov(evenness~treatment,data=cvs.diversity)

anova(m3)
anova(m4)
anova(m5)
anova(m6)

car::Manova(lm(cbind(richness, InvSimpson, evenness) ~ treatment, data = cvs.diversity))
SimTestDiff(resp = c("z_density","richness","InvSimpson","evenness"),
            data = as.data.frame(cvs.diversity), grp = "treatment", type="Tukey")

TukeyHSD(m4)
TukeyHSD(m5)
TukeyHSD(m6)

bwplot(z_density ~ treatment, data = cvs.diversity, ylab = "CV zooplankton density (#/mL)")
bwplot(richness ~ treatment, data = cvs.diversity, ylab = "CV richness") # was rarefied_S
bwplot(InvSimpson ~ treatment, data = cvs.diversity, ylab = "CV inverse Simpson's D")
bwplot(evenness ~ treatment, data = cvs.diversity, ylab = "CV Simpson's E")

cvs.summarised <- cvs.diversity %>%
  group_by(treatment) %>%
  select(-block) %>%
  summarise_all(.funs = c(mean = ~mean(.),
                SE = ~1*sd(.)/sqrt(10)))

cvs.summarised

rm(m1, m2, m3, m4, m5, m6)

## Plot CV of Simpson's diversity & evenness ----

# base plot for CVs
cv.plot <- ggplot(cvs.summarised, aes(x = treatment)) +
  scale_x_discrete(labels = c("Control", "Thinning", "Predation")) +
  theme_classic(base_size = 18)

# Simpson's Diversity plot
cv.D <- cv.plot + geom_point(aes(y = InvSimpson_mean), size = 2) +
  geom_errorbar(aes(ymin = InvSimpson_mean - InvSimpson_SE, ymax = InvSimpson_mean + InvSimpson_SE),
                width = 0.05) +
  coord_cartesian(ylim = c(0.15, 0.45)) +
  labs(y = "CV of Simpson's Diversity", x = "Treatment")
cv.D

# Simpson's Evenness plot
cv.E <- cv.plot + geom_point(aes(y = evenness_mean), size = 2) +
  geom_errorbar(aes(ymin = evenness_mean - evenness_SE, ymax = evenness_mean + evenness_SE),
                width = 0.05) +
  coord_cartesian(ylim = c(0.15, 0.45)) +
  labs(y = "CV of Simpson's Evenness", x = "Treatment")
cv.E

# label each plot as a panel
cvD <- arrangeGrob(cv.D, bottom = textGrob("(a)", x = unit(0, "npc")
                                                , y   = unit(2, "npc"), just=c("left","top"),
                                                gp=gpar(col="black", fontsize=14)))
cvE <- arrangeGrob(cv.E, bottom = textGrob("(b)", x = unit(0, "npc")
                                           , y   = unit(2, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=14)))
# arrange panels into single figure
grid.arrange(cvD, cvE, nrow = 1)

# clean up
rm(cvs.diversity, cvs.enclosure, cvs.enclosure.long, cvs.stocked, 
   cvs.stocked.long, cvs.summarised, cvs.treatment, cvs.whole)
rm(cv.D, cv.E, cv.plot, cvD, cvE, dodge)

### Diversity Indices ----

# check significance of evenness differences using invsimpson/richness to get evenness alone
# check sensitivity of significance to individual taxa - e.g. try removing something like simo, retest
# opacum might be choosing simo
# evenness increases in control b/c ___ becomes abundant - but that's not the case w/ salamanders
# salamader preferences are having an effect that differs from manual, here it decreases evenness

# go back to original predictions - this structures the paper
# H1: did predators XXX? No. Yes Maybe. Effects on community but doesn't match diversity regulation suggested by feeding trials
# in wild it could go either way, I found YYY despite potential for regulation

# model: max predation coefficient on a given taxon, scale that by proportion of taxon in system

# MARS - multiple auto regressive state space
# https://journal.r-project.org/archive/2012-1/RJournal_2012-1_Holmes~et~al.pdf


# ## convert sample abundance to sample density - assuming sampler skirt was FLUSH not high
# community_keyed <- community_keyed / data_keyed$sample_vol_flush * data_keyed$vol_pre
# community_zoop <- community_zoop / data_zoop$sample_vol_flush * data_zoop$vol_pre
# community_stock <- community_stock / data_seeded$sample_vol_flush * data_seeded$vol_pre

### summary visualizations ----

## lattice plots of effects across blocks ----
# export xyplots in 1024*600

xyplot(nTot ~ TimeCat | block + treatment, data = data_keyed, type = c("b"))
xyplot(nTot ~ TimeCat | block + treatment, data = data_zoop, type = c("b"))
xyplot(nTot ~ TimeCat | block + treatment, data = data_seeded, type = c("b"))

xyplot(z_density ~ TimeCat | block + treatment, data = data_zoop, type = c("b"))

xyplot(richness ~ TimeCat | block + treatment, data = data_keyed, type=c("b")) # all taxa
xyplot(richness ~ TimeCat | block + treatment, data = data_zoop, type=c("b")) # zoop only
xyplot(richness ~ TimeCat | block + treatment, data = data_seeded, type=c("b")) # stocked spp only

xyplot(InvSimpson ~ TimeCat | block + treatment, data = data_keyed, type=c("b"))  # all taxa
xyplot(InvSimpson ~ TimeCat | block + treatment, data = data_zoop, type=c("b")) # zoop only
xyplot(InvSimpson ~ TimeCat | block + treatment, data = data_seeded, type=c("b")) # stocked spp only

## lattice box-and-whisker plots, summarizing by treatment ----

# stocked spp only
bwplot(nTot ~ as.factor(TimeCat) | treatment, data=data_seeded, aspect=2, main = "Stocked species") # abundance
bwplot(z_density ~ as.factor(TimeCat) | treatment, data=data_seeded, aspect=2, main = "Stocked species") # density
bwplot(log(z_density) ~ as.factor(TimeCat) | treatment, data=data_seeded, aspect=2, main = "Stocked species") # log density
bwplot(richness ~ as.factor(TimeCat) | treatment, data=data_seeded, aspect=2) # Richness
bwplot(rarefied_S ~ as.factor(TimeCat) | treatment, data=data_seeded, aspect=2, main = "Stocked species") # density
bwplot(InvSimpson ~ as.factor(TimeCat) | treatment, data=data_seeded, aspect=2) # Inverse Simpson's D
bwplot(InvSimpson / richness ~ as.factor(TimeCat) | treatment, data=data_seeded, aspect=2) # Simpson's E
bwplot(InvSimpson / rarefied_S ~ as.factor(TimeCat) |treatment, data=data_seeded, aspect=2) # rarefied Simpson's E

# zoop only
bwplot(nTot~as.factor(TimeCat)|treatment, data=data_zoop, aspect=2,main="zooplankton") # abundance
bwplot(z_density~as.factor(TimeCat)|treatment, data=data_zoop, aspect=2,main="zooplankton") # log density
bwplot(log(z_density)~as.factor(TimeCat)|treatment, data=data_zoop, aspect=2,main="zooplankton") # density
bwplot(richness~as.factor(TimeCat)|treatment, data=data_zoop, aspect=2) # Richness
bwplot(rarefied_S~as.factor(TimeCat)|treatment, data=data_zoop, aspect=2) # Richness
bwplot(InvSimpson~as.factor(TimeCat)|treatment, data=data_zoop, aspect=2) # Inverse Simpson's D
bwplot(InvSimpson/richness~as.factor(TimeCat)|treatment, data=data_zoop, aspect=2) # Simpson's E
bwplot(InvSimpson/rarefied_S~as.factor(TimeCat)|treatment, data=data_zoop, aspect=2) # rarefied Simpson's E

# all taxa
bwplot(nTot~as.factor(TimeCat)|treatment, data=data_keyed, aspect=2) # abundance
bwplot(z_density~as.factor(TimeCat)|treatment, data=data_keyed, aspect=2) # density
bwplot(log(z_density)~as.factor(TimeCat)|treatment, data=data_keyed, aspect=2) # log density
bwplot(richness~as.factor(TimeCat)|treatment, data=data_keyed, aspect=2) # Richness
bwplot(rarefied_S~as.factor(TimeCat)|treatment, data=data_keyed, aspect=2) # Richness
bwplot(InvSimpson~as.factor(TimeCat)|treatment, data=data_keyed, aspect=2) # Inverse Simpson's D
bwplot(InvSimpson/richness~as.factor(TimeCat)|treatment, data=data_keyed, aspect=2) # Simpson's E
bwplot(InvSimpson/rarefied_S~as.factor(TimeCat)|treatment, data=data_keyed, aspect=2) # rarefied Simpson's E

### view individual taxa ----

## view by block
xyplot(Asplanchna~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Candona_sp_1~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Candona_sp_2~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Chydorus~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Cyclocypris~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Cyclopoid~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Cypricercus~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Euchlanis~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Harpacticoid~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Isotomidae~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Keratella~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Lecane~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Lepadella~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Pleuroxus.denticulatus~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Scapholeberis~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Simocephalus~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Sminthuridae~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(Trichocerca~TimeCat|block+treatment,data=data_zoop,type="b")
xyplot(unknown_rotifer_1~TimeCat|block+treatment,data=data_zoop,type="b")

## summarize by treatment

# scale by nTot, or log transform the abundances (bounded by 0)
bwplot(Asplanchna/nTot~as.factor(TimeCat)|treatment,data=data_zoop)
bwplot(Candona_sp_1/nTot~as.factor(TimeCat)|treatment,data=data_zoop)
bwplot(Candona_sp_2/nTot~as.factor(TimeCat)|treatment,data=data_zoop)
bwplot(Chydorus/nTot~as.factor(TimeCat)|treatment,data=data_zoop) #this one
bwplot(Chydorus/nTot~as.factor(TimeCat)|treatment,data=data_seeded) #this one
bwplot(log(Chydorus)~as.factor(TimeCat)|treatment,data=data_zoop) #this one
bwplot(log(Chydorus)~as.factor(TimeCat)|treatment,data=data_seeded) #this one

bwplot(Cyclocypris/nTot~as.factor(TimeCat)|treatment,data=data_zoop)
bwplot(Cyclocypris/nTot~as.factor(TimeCat)|treatment,data=data_seeded)

bwplot(Cyclopoid/nTot~as.factor(TimeCat)|treatment,data=data_zoop) #this one
bwplot(Cyclopoid/nTot~TimeCat|treatment,data=data_seeded) #this one

bwplot(Cypricercus/nTot~as.factor(TimeCat)|treatment,data=data_zoop) #this one
bwplot(Cypricercus/nTot~as.factor(TimeCat)|treatment,data=data_seeded) #this one

bwplot(Euchlanis/nTot~as.factor(TimeCat)|treatment,data=data_zoop)
bwplot(Harpacticoid/nTot~as.factor(TimeCat)|treatment,data=data_zoop) #this one
bwplot(Harpacticoid/nTot~as.factor(TimeCat)|treatment,data=data_seeded) #this one

bwplot(Isotomidae/nTot~as.factor(TimeCat)|treatment,data=data_zoop)
bwplot(log(Keratella/sample_vol_flush)~as.factor(TimeCat)|treatment,data=data_zoop) #this one
bwplot(Keratella/nTot~as.factor(TimeCat)|treatment,data=data_zoop) #this one

bwplot(Lecane/nTot~as.factor(TimeCat)|treatment,data=data_zoop)
bwplot(Lepadella/nTot~as.factor(TimeCat)|treatment,data=data_zoop) #this one
bwplot(Pleuroxus.denticulatus/nTot~as.factor(TimeCat)|treatment,data=data_zoop) #sort-of?
bwplot(Scapholeberis/nTot~as.factor(TimeCat)|treatment,data=data_zoop) #sort-of?
bwplot(Simocephalus/sample_vol_flush~as.factor(TimeCat)|treatment,data=data_zoop) #this one
bwplot(Simocephalus/nTot~as.factor(TimeCat)|treatment,data=data_zoop) #this one
bwplot(Simocephalus/nTot~as.factor(TimeCat)|treatment,data=data_seeded) #this one

bwplot(Trichocerca/nTot~as.factor(TimeCat)|treatment,data=data_zoop)
bwplot(unknown_rotifer_1/nTot~as.factor(TimeCat)|treatment,data=data_zoop)

by(data_zoop, list(data_zoop$TimeCat, data_zoop$treatment), function(x){
  y <- subset(x, select = c(Simocephalus))
  z <- subset(x, select = c(nTot))
  w <- y/z
  apply(w, 2, mean)
}
)
# frequency: control: 0.2341503/0.0250856
# frequency: manual: 0.1914459/0.02653027
# density: control: 0.06349303/0.001300172
# density: manual: 0.05142766/0.001801125

# to compare chydorus, simocephalus, and keratella frequencies over time,
# species counts must be gathered into a single column (i.e., long form table)

data_zoop.long <- data_zoop %>% 
  gather(key = "species", value = "count", Asplanchna:unknown_rotifer_1) %>%
  select(V1:block, Time:TimeCat, sample_vol_flush, nTot, species, count) %>%
  filter(species %in% c("Chydorus", "Simocephalus", "Keratella")) %>%
  mutate(species.f = factor(species, levels = c("Simocephalus", "Chydorus", "Keratella")))

# now plot frequency as species counts divided by enclosure totals (nTot)
ggplot(data_zoop.long, aes(x = as.integer(TimeCat), y = count/nTot, color = treatment)) +
  facet_wrap(. ~ species.f, nrow = 1) +
  stat_summary(fun.y = mean, geom = "point", position = dodge) +
  stat_summary(fun.y = mean, geom = "line", position = dodge) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = dodge, width = 0.7) +
  scale_x_continuous(breaks = breaks,
                     labels = c("May 5", "May 15", "May 27", "June 9")) +
  scale_color_discrete(name = "Treatment") +
  labs(x = "Date", y = "Frequency") +
  theme_bw(base_size=size)+
  theme(axis.text.x = element_text(vjust=1, hjust = 1, angle = 45), 
        panel.spacing = unit(1, "lines"), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank())

### Investigate correlations among taxa - possible food web dynamics? ----
#
# are there correlations among taxa?
# food web model? repeatable pattern if A increases then B decreases? (correlation)
# this would take a matrix of interactions; bayesian techniques to search for interaction coefficients
#
## exclude the salamander treatments to remove influence of salamander predation (external forcing) on dynamics

# zooplankton only
# subset just the individuals that appear in 5+ samples
data.common <- subset(data_zoop, select = which(colSums(data_zoop[,1:18] != 0) > 4))
M <- cor(data.common[-seq(from=3, to=120, by=3),2:13]) # exclude salamander treatments
corrplot(M, method = "circle", type = "upper", order = "hclust")
corrplot(M, method = "circle", type = "upper", order = "AOE")

rm(data.common, M)

### Mixed Effects Modeling ----

### log(Density)
dm1<- lmer(z_density~treatment * days + (1|block) + (1|TimeCat), data=data_zoop)
dm2<- lmer(z_density~treatment + days + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(dm1, dm2)) # not significant, ok to remove the interaction

dm3<- lmer(log(z_density)~days + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(dm3, dm2)) # significant, don't remove treatment

dm4<- lmer(log(z_density)~treatment + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(dm4, dm2)) # significant, don't remove Time

## don't remove random effects since there's structure, but can still test it w/ a 
#likelihood ratio test, i.e. simple ANOVA

# block effect
dm5<- lmer(z_density~treatment + Time + (1|TimeCat), data=data_zoop)
anova(dm2,dm5)
# block effect significant

# time effect
dm6<- lmer(z_density~treatment + Time + (1|block), data=data_zoop)
anova(dm2,dm6)
# time effect not significant

summary(dm1)

# visualize, export 600*450
bwplot(log(z_density)~days|treatment,data=data_zoop, xlab=list(label="Treatment",cex=1.4),
       ylab=list(label="log zooplankton density (N/mL)",cex=1.4), scales=list(cex=1.1,x=list(rot=45)), 
       par.strip.text=list(cex=1.4))

dodge=position_dodge(width=0.7)
ggplot(data_zoop, aes(x=as.integer(TimeCat), y=z_density, color=treatment))+
  stat_summary(fun.y = mean,geom="point",size=2,position=dodge)+
  stat_summary(fun.y = mean,geom="line",aes(group=treatment),linetype=2,position=dodge)+
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.7,position=dodge)+
  xlab("Date")+
  ylab("Zooplankton density (#/mL)")+
  scale_x_continuous(breaks = breaks,
                     labels = c("May 5", "May 15", "May 27", "June 9"),
                     minor_breaks = NULL) +
  scale_color_discrete(name = "Treatment",
                       labels = c("Control", "Thinning", "Predation"))+
  # scale_y_log10(breaks = c(0.1, 0.25, 0.5, 0.75, 1),
                # minor_breaks = NULL)+
  # coord_trans(y = "log10")+
  theme_bw(base_size=18)+
  expand_limits(y = 0.06) +
  theme(axis.text.x = element_text(hjust=1, angle = 45, vjust=1))

rm(dm1, dm2, dm3, dm4, dm5, dm6, kr.b)

### Richness
rm1<- lmer(richness~treatment * days + (1|block) + (1|TimeCat), data=data_zoop)
rm2<- lmer(richness~treatment + days + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(rm1, rm2)) # not significant, ok to remove the interaction

rm3<- lmer(richness~days + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(rm3, rm2)) # significant, don't remove treatment

rm4<- lmer(richness~treatment + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(rm4, rm2)) # not significant, ok to remove Time

rm5<- lmer(richness~1 + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(rm5, rm4)) # significant, still don't remove treatment

summary(rm4)

Treatment<-factor(c("Control","Thinning","Predation"))
Treatment<- factor(Treatment, levels=c("Control","Thinning","Predation"))
Mean<-c(6.025,6.025+0.5,6.025-0.45) # zoop
SE<-c(0.284,0.2453,0.2453) # zoop
rs<-data.frame(cbind(Treatment,Mean,SE ))
rs$Treatment<-Treatment

## visualize

ggplot(data_zoop, aes(x=richness,fill=treatment))+geom_density(alpha=0.3)+
  ggtitle("Richness by Treatment, all zooplankton")+xlab("Richness")+theme_bw()


ggplot(rs, aes(x=Treatment, y=Mean, ymax=6.5,ymin=5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE),width=0.05)+
  coord_cartesian(ylim = c(4.5, 7)) +
  theme_classic(base_size=18)+
  labs(y = "Richness")+
  scale_x_discrete(labels=Treatment)


# box-&-whiskers plot from lattice
bwplot(richness~treatment, data=data_seeded,xlab=list(label="Treatment",cex=1.4),
       ylab=list(label="Richness",cex=1.4),scales=list(cex=1.1),
       par.strip.text=list(cex=1.3)) # Richness - export 450*450

rm(rm1, rm2, rm3, rm4, rm5, kr.b, Mean, SE, rs)
## keep error in since there's structure, but can still test it - likelihood ratio test, i.e. simple ANOVA

# ## lmertest
# # to test (fixed? random?) effects, turn off REML, use just max likelihood
# # library(lmerTest)
# rm1s<- lmer(richness~treatment * Time + (1|block) + (1|TimeCat), data=data_seeded)
# rm1z<- lmer(richness~treatment * Time + (1|block) + (1|TimeCat), data=data_zoop)
# 
# summary(rm1s)
# summary(rm1z)
# anova(rm1s)
# anova(rm1z)
# 
# if(requireNamespace("pbkrtest", quietly = TRUE))
#   anova(rm1, ddf = "Kenward-Roger")
# anova(rm1, ddf="lme4")
# 
# st.1s <- step(rm1s, ddf="Kenward-Roger")
# plot(st.1s)
# st.1z <- step(rm1z, ddf="Kenward-Roger")
# plot(st.1z)
# st.1z
# 
# rm4s<- lmer(richness~treatment + (1|block) + (1|TimeCat), data=data_seeded)
# rm4z<- lmer(richness~treatment + (1|block) + (1|TimeCat), data=data_zoop)
# summary(rm4)
# anova(rm4)
# if(requireNamespace("pbkrtest", quietly = TRUE))
#   anova(rm4, ddf = "Kenward-Roger")
# anova(rm4, ddf="lme4")
# st4s<- step(rm4s, ddf="Kenward-Roger")
# st4z<- step(rm4z, ddf="Kenward-Roger")
# plot(st4s)
# plot(st4z)

### Rarefied richness
rm1<- lmer(rarefied_S~treatment * days + (1|block) + (1|TimeCat), data=data_zoop)
rm2<- lmer(rarefied_S~treatment + days + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(rm1, rm2)) # not significant, ok to remove the interaction

rm3<- lmer(rarefied_S~days + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(rm3, rm2)) # significant, don't remove treatment

rm4<- lmer(rarefied_S~treatment + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(rm4, rm2)) # not significant, ok to remove Time
summary(rm4)

Mean<-c(4.94741,4.94741+0.07332,4.94741-0.90353) # rarefied by time
SE<-c(0.30518,0.15756,0.15756)
Mean33<-c(4.15394,4.15394+0.07933,4.15394-0.91860) # rarefied to 33
SE33<-c(0.12135,0.12401,0.12401)
rs<-data.frame(cbind(Treatment,Mean,SE, Mean33, SE33))
rs$Treatment<-Treatment

r.s <- ggplot(rs, aes(x=Treatment, y=Mean))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE),width=0.05)+
  coord_cartesian(ylim = c(2.5,5.5))+
  theme_classic(base_size=18)+
  labs(y = "Rarefied richness")+
  scale_x_discrete(labels=Treatment)

r.s.33 <- ggplot(rs, aes(x=Treatment, y=Mean33))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Mean33-SE33, ymax=Mean33+SE33),width=0.05)+
  coord_cartesian(ylim = c(2.5,5.5))+
  theme_classic(base_size=18)+
  labs(y = "Rarefied richness")+
  scale_x_discrete(labels=Treatment)

rsa <- arrangeGrob(r.s, bottom = textGrob("(a)", x = unit(0, "npc"),
                                           y   = unit(2, "npc"), just=c("left","top"),
                                           gp=gpar(col="black", fontsize=14)))
rsb <- arrangeGrob(r.s.33, bottom = textGrob("(b)", x = unit(0, "npc"),
                                          y   = unit(2, "npc"), just=c("left","top"),
                                          gp=gpar(col="black", fontsize=14)))

# arrange panels into single figure
grid.arrange(rsa, rsb, nrow = 1)

## test random effects

# block effect
rm5<- lmer(rarefied_S~treatment + (1|TimeCat), data=data_zoop)
anova(rm4, rm5)
# block effect significant

# time effect
rm6<- lmer(rarefied_S~treatment + (1|block), data=data_zoop)
anova(rm4,rm6)
# time not significant

## visualize
bwplot(rarefied_S~treatment, data=data_seeded,xlab=list(label="Treatment",cex=1.4),
       ylab=list(label="Rarefied richness",cex=1.4),scales=list(cex=1.1),
       par.strip.text=list(cex=1.3)) # Richness - export 450*450

rm(rm1, rm2, rm3, rm4, kr.b, Mean, SE, rs)
rm(r.s, r.s.33, rsa, rsb)

### Inverse Simpson's D
sm1<- lmer(InvSimpson~treatment * days + (1|block) + (1|TimeCat), data=data_zoop)
sm2<- lmer(InvSimpson~treatment + days + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(sm1, sm2)) # significant, retain the interaction of treatment:Time

sm3<-lmer(InvSimpson~days + (1|block) + (1|TimeCat), data=data_seeded)
(kr.b<-KRmodcomp(sm3, sm2)) # significant, retain treatment

sm4<-lmer(InvSimpson~treatment + (1|block) + (1|TimeCat), data=data_seeded)
(kr.b<-KRmodcomp(sm4, sm2)) # not significant, remove time


dm1s<- lmer(InvSimpson~treatment * Time + (1|block) + (1|TimeCat), data=data_seeded)
dm1z<- lmer(InvSimpson~treatment * Time + (1|block) + (1|TimeCat), data=data_zoop)

st.1d.s <- step(dm1s, ddf="Kenward-Roger")
plot(st.1d.s)
st.1d.s

st.1d.z <- step(dm1z, ddf="Kenward-Roger")
plot(st.1d.z)
st.1d.z
summary(dm1z)
anova(dm1z)
anova(dm1z, ddf="Kenward-Roger")
lsmeans(dm1z,list(pairwise~treatment),adjust="tukey")

# # visualize: density plot

# box-&-whiskers plot from lattice
bwplot(InvSimpson~TimeCat|treatment, data=data_zoop, aspect=2,xlab=list(label="Time",cex=1.4),
       ylab=list(label="Inverse Simpson's D",cex=1.4),scales=list(cex=1.1,x=list(rot=45)),
       par.strip.text=list(cex=1.4))

bwplot(InvSimpson~TimeCat|treatment, data=data_seeded, aspect=2,xlab=list(label="Time",cex=1.4),
       ylab=list(label="Inverse Simpson's D",cex=1.4),scales=list(cex=1.1,x=list(rot=45)),
       par.strip.text=list(cex=1.4))

bwplot(InvSimpson~TimeCat|treatment, data=data.noSimo, aspect=2,xlab=list(label="Time",cex=1.4),
       ylab=list(label="Inverse Simpson's D",cex=1.4),scales=list(cex=1.1,x=list(rot=45)),
       par.strip.text=list(cex=1.4))

bwplot(InvSimpson~TimeCat|treatment, data=data.noKera, aspect=2,xlab=list(label="Time",cex=1.4),
       ylab=list(label="Inverse Simpson's D",cex=1.4),scales=list(cex=1.1,x=list(rot=45)),
       par.strip.text=list(cex=1.4))


dodge = position_dodge(width=.7)
breaks = unique(as.integer(data_zoop$TimeCat))

ggplot(data_zoop, aes(x=as.integer(TimeCat), y=InvSimpson, color=treatment))+
  stat_summary(fun.y = mean, geom="point", size=2, position=dodge)+
  stat_summary(fun.y = mean, geom="line", aes(group=treatment), linetype=2, position=dodge)+
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.7, position=dodge)+
  xlab("Date")+
  ylab("Simpson's Diversity")+
  scale_x_continuous(breaks = breaks,
                     labels = c("May 5", "May 15", "May 27", "June 9")) +
  scale_color_discrete(name = "Treatment",
                       labels = c("Control", "Thinning", "Predation"))+
  expand_limits(y = 0) +
  theme_classic(base_size=18)+
  theme(axis.text.x = element_text(vjust=1, hjust = 1, angle = 45))

ggplot(data_seeded, aes(x=treatment,y=InvSimpson))+
  stat_summary(fun.y = mean,geom="point",size=2)+
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.05)+
  scale_x_discrete(name="Treatment",labels=rownames)+
  ylab("Inverse Simpson's D")+
  theme_bw(base_size=18)
  
data.stock.noSimo<-data_seeded
data.stock.noSimo$Simocephalus<-NULL
data.stock.noSimo$nTot<-rowSums(data.stock.noSimo[,2:6]) # columns 2:17 are individual taxon counts after restriction
data.stock.noSimo$InvSimpson<-1/rowSums((data.stock.noSimo[,2:6]/data.stock.noSimo$nTot)^2)

ggplot(data.stock.noSimo, aes(x=treatment,y=InvSimpson))+
  stat_summary(fun.y = mean,geom="point",size=2)+
  stat_summary(fun.data = mean_se,geom="errorbar",width=0)+
  xlab("Treatment")+
  ylab("Inverse Simpson's D")+
  theme_bw(base_size=18)


### Evenness: Simpson's E
data_zoop$evenness<-data_zoop$InvSimpson/data_zoop$richness
em1<- lmer(evenness~treatment * days + (1|block) + (1|TimeCat), data=data_zoop)
em2<- lmer(evenness~treatment + days + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(em1, em2)) # significant, retain the interaction
summary(em1)

em3<- lmer(evenness~days + (1|block) + (1|TimeCat), data=data_seeded)
(kr.b<-KRmodcomp(em3, em2)) # NS, remove treatment

em4<- lmer(evenness~(1|block) + (1|TimeCat), data=data_seeded)
(kr.b<-KRmodcomp(em3, em4)) # NS, remove time

dodge=position_dodge(width=0.7)
ggplot(data_zoop, aes(x=as.integer(TimeCat), y=evenness, color=treatment))+
  stat_summary(fun.y = mean,geom="point",size=2,position=dodge)+
  stat_summary(fun.y = mean,geom="line",aes(group=treatment),linetype=2,position=dodge)+
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.7,position=dodge)+
  scale_x_continuous(breaks = breaks,
                     labels = c("May 5", "May 15", "May 27", "June 9")) +
  xlab("Date")+
  ylab("Simpson's Evenness")+
  scale_color_discrete(name="Treatment",
                       labels=c("Control", "Thinning", "Predation"))+
  theme_classic(base_size=18)+
  expand_limits(y = 0) +
  theme(axis.text.x = element_text(hjust=1, angle = 45, vjust=1))

em1z<- lmer(evenness~treatment * Time + (1|block) + (1|TimeCat), data=data_zoop)
st.1e.z <- step(em1z)
plot(st.1e.z)
st.1e.z

# visualize: density plot
ggplot(data_zoop, aes(InvSimpson/richness))+geom_density(fill="#4271AE",alpha=0.5)+facet_grid(.~Time*treatment)+
  ggtitle("Evenness by Treatment, all plankton")+xlab("Simpson's E")+theme_light()+
  facet_wrap(facets=~Time*treatment,nrow=4)

# box-&-whiskers plot from lattice
bwplot(evenness~TimeCat|treatment, data=data_zoop, aspect=2,xlab=list(label="Time",cex=1.4),
       ylab=list(label="Simpson's E",cex=1.4),scales=list(cex=1.1,x=list(rot=45)),
       par.strip.text=list(cex=1.4)) # export 600*450

em3<- lmer(InvSimpson~Time + (1|block) + (1|TimeCat), data=data_zoop)
(kr.b<-KRmodcomp(em1, em3)) # significant, definitely don't remove both interaction and Time
(kr.b<-KRmodcomp(em2, em3)) # marginally significant


# remove Simocephalus ----
data.noSimo<-data_zoop
data.noSimo$Simocephalus<-NULL
data.noSimo$nTot<-rowSums(data.noSimo[,2:17]) # columns 2:17 are individual taxon counts after restriction
data.noSimo$richness<-rowSums(data.noSimo[,2:17]!=0)
data.noSimo$InvSimpson<-1/rowSums((data.noSimo[,2:17]/data.noSimo$nTot)^2)

# Inverse Simpson's D
sm1<- lmer(InvSimpson~treatment * days + (1|block) + (1|TimeCat), data=data.noSimo)
sm2<- lmer(InvSimpson~treatment + days + (1|block) + (1|TimeCat), data=data.noSimo)
(kr.b<-KRmodcomp(sm1, sm2)) # NS, remove interaction treatment:Time. This is DIFFERENT!

sm3<- lmer(InvSimpson~days + (1|block) +(1|TimeCat), data=data.noSimo)
(kr.b<-KRmodcomp(sm2, sm3)) # significant, don't remove treatment

sm4<- lmer(InvSimpson~treatment + (1|block) +(1|TimeCat), data=data.stock.noSimo)
(kr.b<-KRmodcomp(sm2, sm4)) # NS, remove Time. This is DIFFERENT from full data!

# visualize
ggplot(data.noSimo, aes(InvSimpson))+geom_density(fill="#4271AE",alpha=0.5)+facet_grid(.~treatment)+
  ggtitle("excluding Simocephalus")+coord_flip()+xlab("Diversity")
ggplot(data.noSimo, aes(x=InvSimpson,fill=treatment))+geom_density(alpha=0.3)+
  ggtitle("excluding Simocephalus")+coord_flip()+xlab("Diversity")
bwplot(InvSimpson~TimeCat|treatment, data=data.noSimo, aspect=2,xlab=list(label="Time",cex=1.4),
       ylab=list(label="Inverse Simpson's D",cex=1.4),scales=list(cex=1.1,x=list(rot=45)),
       par.strip.text=list(cex=1.4))
bwplot(InvSimpson~treatment, data=data.noSimo,xlab=list(label="Time",cex=1.4),
       ylab=list(label="Inverse Simpson's D",cex=1.4),scales=list(cex=1.1))

# Evenness: Simpson's E
em1<- lmer(InvSimpson/richness~treatment * days + (1|block) + (1|TimeCat), data=data.noSimo)
em2<- lmer(InvSimpson/richness~treatment + days + (1|block) + (1|TimeCat), data=data.noSimo)
(kr.b<-KRmodcomp(em1, em2)) # NS, remove treatment:Time interaction. This is DIFFERENT!

em3<- lmer(InvSimpson/richness~days + (1|block) + (1|TimeCat), data=data.noSimo)
(kr.b<-KRmodcomp(em3,em2)) # NS, remove treatment. VERY DIFFERENT
em4<- lmer(InvSimpson/richness~treatment + (1|block) + (1|TimeCat), data=data.noSimo)
(kr.b<-KRmodcomp(em4,em2)) # NS, remove Time. DIFFERENT

em5<-lmer(InvSimpson/richness~ (1|block) + (1|TimeCat), data=data.noSimo)
(kr.b<-KRmodcomp(em5,em3))
(kr.b<-KRmodcomp(em5,em4)) # all NS, remove all fixed effects. VERY DIFFERENT
plot(em1)
summary(em5)
# excluding Simocephalus removes treatment*Time interaction and Time main fixed effect from iS
# also removes ALL fixed effects from Simpson's E

# remove Chydorus ----
data.noChyd<-data_zoop
data.noChyd$Chydorus<-NULL
data.noChyd$nTot<-rowSums(data.noChyd[,2:17]) # columns 2:17 are individual taxon counts after restriction
data.noChyd$nTot
data_zoop$nTot
data.noChyd$richness<-rowSums(data.noChyd[,2:17]!=0)
data.noChyd$richness
data_zoop$richness
data.noChyd$InvSimpson<-1/rowSums((data.noChyd[,2:17]/data.noChyd$nTot)^2)
data.noChyd$InvSimpson
data_zoop$InvSimpson

# Inverse Simpson's D
sm1<- lmer(InvSimpson~treatment * Time + (1|block) + (1|TimeCat), data=data.noChyd)
sm2<- lmer(InvSimpson~treatment + Time + (1|block) + (1|TimeCat), data=data.noChyd)
(kr.b<-KRmodcomp(sm1, sm2)) # keep interaction

# Evenness: Simpson's E
em1<- lmer(InvSimpson/richness~treatment * Time + (1|block) + (1|TimeCat), data=data.noChyd)
em2<- lmer(InvSimpson/richness~treatment + Time + (1|block) + (1|TimeCat), data=data.noChyd)
(kr.b<-KRmodcomp(em1, em2)) # keep interaction

# excluding Chydorus does NOT alter significance of differences by treatment or Time

# remove Keratella ----
data.noKera<-data_zoop
data.noKera$Keratella<-NULL
data.noKera$nTot<-rowSums(data.noKera[,2:17]) # columns 2:17 are individual taxon counts after restriction
data.noKera$nTot
data_zoop$nTot
data.noKera$richness<-rowSums(data.noKera[,2:17]!=0)
data.noKera$richness
data_zoop$richness
data.noKera$InvSimpson<-1/rowSums((data.noKera[,2:17]/data.noKera$nTot)^2)
data.noKera$InvSimpson
data_zoop$InvSimpson


# Inverse Simpson's D
sm1<- lmer(InvSimpson~treatment * Time + (1|block) + (1|TimeCat), data=data.noKera)
sm2<- lmer(InvSimpson~treatment + Time + (1|block) + (1|TimeCat), data=data.noKera)
(kr.b<-KRmodcomp(sm1, sm2)) # no longer significant - DIFFERENT

sm3<- lmer(InvSimpson~Time + (1|block) + (1|TimeCat), data=data.noKera)
(kr.b<-KRmodcomp(sm3, sm2)) # keep treatment

sm4<- lmer(InvSimpson~treatment + (1|block) + (1|TimeCat), data=data.noKera)
(kr.b<-KRmodcomp(sm4, sm2)) # NS, remove time - DIFFERENT

# visualize
ggplot(data.noKera, aes(InvSimpson))+geom_density(fill="#4271AE",alpha=0.5)+facet_grid(.~treatment)
ggplot(data.noKera, aes(x=InvSimpson,fill=treatment))+geom_density(alpha=0.3)+
  ggtitle("excluding Keratella")+coord_flip()+xlab("Diversity")

# Evenness: Simpson's E
em1<- lmer(InvSimpson/richness~treatment * Time + (1|block) + (1|TimeCat), data=data.noKera)
em2<- lmer(InvSimpson/richness~treatment + Time + (1|block) + (1|TimeCat), data=data.noKera)
(kr.b<-KRmodcomp(em1, em2)) # NS, remove interaction - DIFFERENT

em3<- lmer(InvSimpson/richness~Time + (1|block) + (1|TimeCat), data=data.noKera)
(kr.b<-KRmodcomp(sm3, sm2)) # keep treatment

em4<- lmer(InvSimpson/richness~treatment + (1|block) + (1|TimeCat), data=data.noKera)
(kr.b<-KRmodcomp(sm4, sm2)) # NS, remove time - DIFFERENT
summary(em4)

# visualize
ggplot(data.noKera, aes(InvSimpson/richness))+geom_density(fill="#4271AE",alpha=0.5)+facet_grid(.~treatment)
ggplot(data.noKera, aes(x=InvSimpson/richness,fill=treatment))+geom_density(alpha=0.3)+
  ggtitle("excluding Keratella")+coord_flip()+xlab("Simpson's Evenness")

# removing Keratella also eliminates significant effect of Time on iS and Simpson's E 


# remove Cyclopoid
data.noCycl<-data_zoop
data.noCycl$Cyclopoid<-NULL
data.noCycl$nTot<-rowSums(data.noCycl[,2:17]) # columns 2:17 are individual taxon counts after restriction
data.noCycl$nTot
data_zoop$nTot
data.noCycl$richness<-rowSums(data.noCycl[,2:17]!=0)
data.noCycl$richness
data_zoop$richness
data.noCycl$InvSimpson<-1/rowSums((data.noCycl[,2:17]/data.noCycl$nTot)^2)
data.noCycl$InvSimpson
data_zoop$InvSimpson

# Inverse Simpson's D
sm1<- lmer(InvSimpson~treatment * Time + (1|block) + (1|TimeCat), data=data.noCycl)
sm2<- lmer(InvSimpson~treatment + Time + (1|block) + (1|TimeCat), data=data.noCycl)
(kr.b<-KRmodcomp(sm1, sm2)) # significant, don't remove interaction


ggplot(data.noCycl, aes(x=InvSimpson,fill=treatment))+geom_density(alpha=0.3)+
  ggtitle("excluding Cyclopoid copepods")+coord_flip()
ggplot(data_zoop, aes(x=InvSimpson,fill=treatment))+geom_density(alpha=0.3)+
  ggtitle("with Cyclopoid copepods")+coord_flip()


# Evenness: Simpson's E
em1<- lmer(InvSimpson/richness~treatment * days + (1|block) + (1|TimeCat), data=data.noCycl)
em2<- lmer(InvSimpson/richness~treatment + days + (1|block) + (1|TimeCat), data=data.noCycl)
(kr.b<-KRmodcomp(em1, em2)) # NS, remove treatment:Time interaction - DIFFERENT

em3<- lmer(InvSimpson/Time~days + (1|block) + (1|TimeCat), data=data.noCycl)
(kr.b<-KRmodcomp(em2, em3)) # significant, don't remove treatment

em4<- lmer(InvSimpson~treatment + (1|block) + (1|TimeCat), data=data.noCycl)
(kr.b<-KRmodcomp(em4, em2)) # NS, ok to remove Time - DIFFERENT

# excluding Cyclopoid copepods eliminates significant effect of Time on iS and Simpson's E



rownames(community_zoop)<-data_zoop$V1
zoop.mds2<-metaMDS(community_zoop,k=2,trymax=80)
### check stress plot
stressplot(zoop.mds2)

zoop.mds3<-metaMDS(community_zoop,k=3)
stressplot(zoop.mds3)

### visualize plot
plot(zoop.mds3)


treat<-c(rep(c("1. Control","1. Thinning","1. Predation"),40))#,
        # rep(c("2. Control","2. Thinning","2. Predation"),10),
        # rep(c("3. Control","3. Thinning","3. Predation"),10),
        # rep(c("4. Control","4. Thinning","4. Predation"),10))
colors<-c(rep(c("gray90","light blue","tomato"),10),
          rep(c("gray60","cornflower blue","orangered"),10),
          rep(c("gray40","blue","red3"),10),
          rep(c("black","dark blue","red4"),10))
ordiplot(zoop.mds2,type="n")
ordihull(zoop.mds2,groups=treat,draw="polygon",col="grey90",label=F)
orditorp(zoop.mds2,display="sites",col=colors,air=0.01,cex=1)
orditorp(zoop.mds3,display="species",col="red",air=0.01)

library(vegan3d)
ordiplot3d(zoop.mds3,choices=c(1:3),col=colors)
orglpoints(zoop.mds3,display="sites",col=colors,air=0.01,cex=1.25)


## investigate chydorus, simocephalus, and keratella

ggplot(data = data_zoop, aes(x = TimeCat, y = Chydorus/sample_vol_flush, col=treatment)) +
  stat_summary(fun.y = mean,geom="point",size=2,position=dodge) +
  stat_summary(fun.y = mean, geom = "line", aes(group = treatment), linetype = 2, position = dodge) +
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.1,position=dodge) +
  scale_color_discrete(name="Treatment",labels=rownames) +
  theme_bw(base_size=18) +
  theme(axis.text.x = element_text(hjust=1, angle = 45, vjust=1)) +
  ylab("Chyd density")+
  xlab("Date")
# facet_wrap(~Species)+

ggplot(data = data_zoop, aes(x = TimeCat, y = Chydorus/nTot, col=treatment)) +
  stat_summary(fun.y = mean,geom="point",size=2,position=dodge) +
  stat_summary(fun.y = mean, geom = "line", aes(group = treatment), linetype = 2, position = dodge) +
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.1,position=dodge) +
  scale_color_discrete(name="Treatment",labels=rownames) +
  theme_bw(base_size=18) +
  theme(axis.text.x = element_text(hjust=1, angle = 45, vjust=1)) +
  ylab("Chyd fraction")+
  xlab("Date")
# facet_wrap(~Species)+

ggplot(data = data_zoop, aes(x = TimeCat, y = Simocephalus/sample_vol_flush, col=treatment)) +
  stat_summary(fun.y = mean,geom="point",size=2,position=dodge) +
  stat_summary(fun.y = mean, geom = "line", aes(group = treatment), linetype = 2, position = dodge) +
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.1,position=dodge) +
  scale_color_discrete(name="Treatment",labels=rownames) +
  theme_bw(base_size=18) +
  theme(axis.text.x = element_text(hjust=1, angle = 45, vjust=1)) +
  ylab("Simo density")+
  xlab("Date")
  # facet_wrap(~Species)+

ggplot(data = data_zoop, aes(x = TimeCat, y = Simocephalus/nTot, col=treatment)) +
  stat_summary(fun.y = mean,geom="point",size=2,position=dodge) +
  stat_summary(fun.y = mean, geom = "line", aes(group = treatment), linetype = 2, position = dodge) +
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.1,position=dodge) +
  scale_color_discrete(name="Treatment",labels=rownames) +
  theme_bw(base_size=18) +
  theme(axis.text.x = element_text(hjust=1, angle = 45, vjust=1)) +
  ylab("Simo fraction")+
  xlab("Date")
# facet_wrap(~Species)+


ggplot(data = data_zoop, aes(x = TimeCat, y = Keratella/sample_vol_flush, col=treatment)) +
  stat_summary(fun.y = mean,geom="point",size=2,position=dodge) +
  stat_summary(fun.y = mean, geom = "line", aes(group = treatment), linetype = 2, position = dodge) +
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.1,position=dodge) +
  scale_color_discrete(name="Treatment",labels=rownames) +
  theme_bw(base_size=18) +
  theme(axis.text.x = element_text(hjust=1, angle = 45, vjust=1)) +
  ylab("Kera density")+
  xlab("Date")


ggplot(data = data_zoop, aes(x = TimeCat, y = Keratella/nTot, col=treatment)) +
  stat_summary(fun.y = mean,geom="point",size=2,position=dodge) +
  stat_summary(fun.y = mean, geom = "line", aes(group = treatment), linetype = 2, position = dodge) +
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.1,position=dodge) +
  scale_color_discrete(name="Treatment",labels=rownames) +
  theme_bw(base_size=18) +
  theme(axis.text.x = element_text(hjust=1, angle = 45, vjust=1)) +
  ylab("Kera fraction")+
  xlab("Date")

