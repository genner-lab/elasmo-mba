#!/usr/bin/env Rscript

#Activate packages
library(here)
library(vegan)
library(ape)
library(sjPlot)
library(ggpubr)
library(mgcv)
library(mgcViz)
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(dplyr)
library(reshape2)
library(indicspecies)

#Read in OTU data, file contains all 209 samples, even if zero elasmobrachs present
Elasmo_Raw <- read.table(here("assets/otu-table.txt"),header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

#Ensure Depth is a factor
Elasmo_Raw$Depth <- as.factor(Elasmo_Raw$Depth)

#Retain a file with all samples, even zeros counts
Elasmo_zeros <- Elasmo_Raw

#Remove rows with all zero values, leaves samples with at least one elasmobranch presnt
Elasmo_counts <- Elasmo_Raw[rowSums(Elasmo_Raw[,6:18])>0,]

#Dataframe of the count data that has been Wisconsin transformed
Elasmo_Wisc <- Elasmo_counts[,6:18]
Elasmo_Wisc <- wisconsin(Elasmo_Wisc)
Elasmo_Wisc <- cbind(Elasmo_counts[,1:5],Elasmo_Wisc)

#Dataframe of the count data that has been Hellinger transformed
Elasmo_Hell <- decostand (Elasmo_counts[,6:18], method="hellinger")
Elasmo_Hell <- cbind(Elasmo_counts[,1:5],Elasmo_Hell)

#Principal Coordinate Analyses (PCOA) on the data matrices
Elasmo_Wisc.D <- vegdist(Elasmo_Wisc[,6:18], "jaccard")
Elasmo_Hell.D <- vegdist(Elasmo_Hell[,6:18], "jaccard")

PCOA_Elasmo_Wisc <- pcoa(Elasmo_Wisc.D, correction="none", rn=NULL)
PCOA_Elasmo_Hell <- pcoa(Elasmo_Hell.D, correction="none", rn=NULL)

#To determine the variance captured by each PCOA axis, see the "Rel_corr_eig" column
PCOA_Elasmo_Wisc$values
PCOA_Elasmo_Hell$values

#Output the PCOA scores for each sample
PCOA_Elasmo_Wisc_Data <- cbind(Elasmo_counts[,1:5],PCOA_Elasmo_Wisc$vectors)
PCOA_Elasmo_Hell_Data <- cbind(Elasmo_counts[,1:5],PCOA_Elasmo_Hell$vectors)

#Plot the PCOA scores (manually enter the variance captured on axis labels)

SpaceWisc <- ggscatter(PCOA_Elasmo_Wisc_Data, x = "Axis.1", y = "Axis.2", color = "localityID",
                      palette = c("#999999", "#E69F00", "#56B4E9"),
                      ellipse = FALSE, ellipse.type = "convex", mean.point = TRUE,
                      star.plot = TRUE) +
  labs(x ="PCOA1 (8.3% of variation)", y = "PCOA2 (7.5% of variation)") +
  theme(legend.position = "none") 
SpaceWisc

SpaceHell <- ggscatter(PCOA_Elasmo_Hell_Data, x = "Axis.1", y = "Axis.2", color = "localityID",
                     palette = c("#999999", "#E69F00", "#56B4E9"),
                     ellipse = FALSE, ellipse.type = "convex", mean.point = TRUE,
                     star.plot = TRUE) +
  labs(x ="PCOA1 (7.1% of variation)", y = "PCOA2 (5.6% of variation)") +
  theme(legend.position = "none") 
SpaceHell

#Combine the two PCOA plots

ggarrange(SpaceWisc, SpaceHell,
          labels = c("a", "b"),
          ncol = 1, nrow = 2, legend = "right",   common.legend = TRUE)

#Statistically test for differences among localities, depths and months, using Permanova

PCOA_Elasmo_Wisc_Data$Month <- as.factor(PCOA_Elasmo_Wisc_Data$Month)
PCOA_Elasmo_Hell_Data$Month <-as.factor(PCOA_Elasmo_Hell_Data$Month)

SpaceTest_Wisc <- adonis2(Elasmo_Wisc.D ~ localityID*Depth*Month, data = Elasmo_Wisc)
SpaceTest_Wisc

SpaceTest_Hell <- adonis2(Elasmo_Hell.D ~ localityID*Depth*Month, data = Elasmo_Wisc)
SpaceTest_Hell

#Pairwise test of differences among site

pairwise.adonis(Elasmo_Wisc[,6:18],Elasmo_Wisc$localityID,sim.method = "jaccard")
pairwise.adonis(Elasmo_Hell[,6:18],Elasmo_Hell$localityID,sim.method = "jaccard")

#Exploring indicator species of each habitat and depth

Wisc_Abund = Elasmo_Wisc[,6:ncol(Elasmo_Wisc)]

Wisc_Depth = Elasmo_Wisc$Depth
Wisc_Depth_sp = multipatt(Wisc_Abund, Wisc_Depth, func = "IndVal", duleg=TRUE, control = how(nperm=9999))
Wisc_Depth_sp

Wisc_Locality = Elasmo_Wisc$localityID
Wisc_Locality_sp = multipatt(Wisc_Abund, Wisc_Locality, func = "IndVal", duleg=TRUE,  control = how(nperm=9999))
Wisc_Locality_sp

Hell_Abund = Elasmo_Hell[,6:ncol(Elasmo_Hell)]

Hell_Depth = Elasmo_Hell$Depth
Hell_Depth_sp = multipatt(Hell_Abund, Hell_Depth, func = "IndVal", duleg=TRUE, control = how(nperm=9999))
Hell_Depth_sp

Hell_Locality = Elasmo_Hell$localityID
Hell_Locality_sp = multipatt(Hell_Abund, Hell_Locality, func = "IndVal", duleg=TRUE,  control = how(nperm=9999))
Hell_Locality_sp

#Plot temporal variation in community structure along the primary axis of PCOA variation

PCOA_Elasmo_Wisc_Data$Month <- as.numeric(PCOA_Elasmo_Wisc_Data$Month)
PCOA_Elasmo_Hell_Data$Month <- as.numeric(PCOA_Elasmo_Hell_Data$Month)

TimeWisc <- ggscatter(PCOA_Elasmo_Wisc_Data, x = "Month", y = "Axis.1", color = "Year") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) + 
  ylim (-0.8,0.8) +
  scale_x_continuous(breaks = get_breaks(n = 10)) +
  theme(legend.position = "none")  + labs(x ="Month", y = "PCOA1")
TimeWisc

TimeHell <- ggscatter(PCOA_Elasmo_Hell_Data, x = "Month", y = "Axis.1", color = "Year") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) + 
  ylim (-0.8,0.8) +
  scale_x_continuous(breaks = get_breaks(n = 10)) +
  theme(legend.position = "none")  + labs(x ="Month", y = "PCOA1")
TimeHell

TimeWisc2 <- ggscatter(PCOA_Elasmo_Wisc_Data, x = "Month", y = "Axis.2", color = "Year") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) + 
  ylim (-0.8,0.8) +
  scale_x_continuous(breaks = get_breaks(n = 10)) +
  theme(legend.position = "none")  + labs(x ="Month", y = "PCOA2")
TimeWisc2

TimeHell2 <- ggscatter(PCOA_Elasmo_Hell_Data, x = "Month", y = "Axis.2", color = "Year") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE) + 
  ylim (-0.8,0.8) +
  scale_x_continuous(breaks = get_breaks(n = 10)) +
  theme(legend.position = "none")  + labs(x ="Month", y = "PCOA2")
TimeHell2

#Final plot 8x8

ggarrange(TimeWisc, TimeWisc2, TimeHell, TimeHell2,
          labels = c("a", "c", "b","d"), 
          ncol = 2, nrow = 2)

#Quantify seasonal pattern in Hellinger data using GAMs (negative binomial), acounting for locality and depth.
#We use the full data for this, including zeros elasmobranch counts.

#Make a dataframe that has been Hellinger transformed
Elasmo_Hell_zeros <- decostand (Elasmo_zeros[,6:18], method="hellinger")
Elasmo_Hell_zeros <- cbind(Elasmo_zeros[,1:5],Elasmo_Hell_zeros)

S_can_eDNA_gam <- gam(Scyliorhinus_canicula ~ s(Month, k=5) + Year + Depth + localityID, family=nb, select=TRUE ,method="REML", data=Elasmo_Hell_zeros)
summary(S_can_eDNA_gam)
S_can_viz <- getViz(S_can_eDNA_gam)
print(plot(S_can_viz))
print(plot(S_can_viz, allTerms = T), pages = 1) 

M_ast_eDNA_gam <- gam(Mustelus_asterias~  s(Month, k=5) + Year + Depth + localityID, family=nb, select=TRUE,method="REML", data=Elasmo_Hell_zeros)
summary(M_ast_eDNA_gam)
M_ast_viz <- getViz(M_ast_eDNA_gam)
print(plot(M_ast_viz), pages = 1)
print(plot(M_ast_viz, allTerms = T), pages = 1) 

R_bra_eDNA_gam <- gam(Raja_brachyura~ s(Month, k=5) + Year + Depth + localityID, family=nb, select=TRUE,method="REML", data=Elasmo_Hell_zeros)
summary(R_bra_eDNA_gam)
R_bra_viz <- getViz(R_bra_eDNA_gam)
print(plot(R_bra_viz), pages = 1)
print(plot(R_bra_viz, allTerms = T), pages = 1) 

R_cla_eDNA_gam <- gam(Raja_clavata~ ss(Month, k=5) + Year + Depth + localityID, family=nb, select=TRUE,method="REML", data=Elasmo_Hell_zeros)
summary(R_cla_eDNA_gam)
R_cla_viz <- getViz(R_cla_eDNA_gam)
print(plot(R_cla_viz), pages = 1)
print(plot(R_cla_viz, allTerms = T), pages = 1) 

R_mic_eDNA_gam <- gam(Raja_microocellata~ s(Month, k=5) + Year + Depth + localityID, family=nb, select=TRUE,method="REML", data=Elasmo_Hell_zeros)
summary(R_mic_eDNA_gam)
R_mic_viz <- getViz(R_mic_eDNA_gam)
print(plot(R_mic_viz), pages = 1)
print(plot(R_mic_viz, allTerms = T), pages = 1) 

R_mon_eDNA_gam <- gam(Raja_montagui~  s(Month, k=5) + Year + Depth + localityID, family=nb, select=TRUE,method="REML", data=Elasmo_Hell_zeros)
summary(R_mon_eDNA_gam)
R_mon_viz <- getViz(R_mon_eDNA_gam)
print(plot(R_mon_viz), pages = 1)
print(plot(R_mon_viz, allTerms = T), pages = 1) 

#Correlating trawl catches and eDNA results.
#Based on pre-prepared file with total eDNA reads (4th rooted), mean CPUE in hauls (4th rooted) and proportion of hauls present, for key periods of time

eDNAvsTrawls <- read.table(here("assets/eDNAvsTrawls.txt"), header=TRUE, fill=TRUE, sep="\t", check.names=FALSE)

#Build the linear models, including the statistical significance of the associations

Reg1 <- lm (TotalReads_4thRoot ~ All_CPUE4th, data = eDNAvsTrawls)
summary (Reg1)

Reg2 <- lm (TotalReads_4thRoot ~ All_PropHaulsPresent, data = eDNAvsTrawls)
summary (Reg2)

#Now we build the plots

set_theme(
  geom.outline.color = "antiquewhite4", 
  geom.outline.size = 1, 
  geom.label.size = 2,
  geom.label.color = "white",
  title.color = "black", 
  axis.title.size = 1.5,
  axis.textcolor = "black", 
  base = theme_classic(base_size = 8))

Reg1ModelPlot <- ggplot(eDNAvsTrawls, aes(y=TotalReads_4thRoot, x=All_CPUE4th)) +
  geom_point(color="black")+
  geom_smooth(method=lm, color="black", fill="steelblue") + ylim (0,100)
Reg1_plot <- Reg1ModelPlot + labs (x="CPUE 1911-2017 (individuals per hour, 4th root)", y="eDNA reads (4th root)") + ggtitle(" ")
Reg1_plot

Reg2ModelPlot <- ggplot(eDNAvsTrawls, aes(y=TotalReads_4thRoot, x=All_PropHaulsPresent)) +
  geom_point(color="black")+
  geom_smooth(method=lm, color="black", fill="steelblue") + ylim (0,100) + xlim (0,80)
Reg2_plot <- Reg2ModelPlot + labs (x="Frequency of occurrence 1911-2017 (% hauls present)", y="eDNA reads (4th root)") + ggtitle(" ")
Reg2_plot

#10x5

figure <- ggarrange(Reg1_plot, Reg2_plot, 
                    labels = c(" ", " "),
                    ncol = 2, nrow = 1)
figure

#Now just a barchart of the species names (save at 9 x 5)
#Based on pre-prepared file with the total reads, and full species names.

TotalReads <- read.table(here("assets/TotalReads.txt"), header=TRUE, fill=TRUE, sep="\t", check.names=FALSE)

set_theme(
  geom.outline.color = "antiquewhite4", 
  geom.outline.size = 1, 
  geom.label.size = 6,
  geom.label.color = "white",
  title.color = "black", 
  axis.title.size = 1.5,
  axis.textcolor = "black", 
  base = theme_classic(base_size = 8))

p <- ggplot(TotalReads, aes(x = reorder(Species, TotalReads), y = TotalReads))
p <- p + geom_bar(stat="identity", color='steelblue',fill='steelblue', width = 0.8,)
p <- p + theme(axis.text.x=element_text(angle=0, hjust=1)) + coord_flip() + scale_y_continuous(labels = scales::comma)
p <- p +labs (x="Species", y="eDNA reads (Total in all samples)")
p

#Now to made the species accumulation per sample images
#First the eDNA accumulation (axes set so these images can combined outside R later)

eDNA_raw <- specaccum(Elasmo_zeros[,6:18], "random")
summary(eDNA_raw)
plot(eDNA_raw, ci.type="poly", col="black", ci.lty=0, ci.col="darkorange", las=1, xlim=c(1, 210), ylim=c(0,30),bg = "transparent")

#Second the trawl accumulation data for 2017, to align the timescale with eDNA sampling (axes set to above)
SpeciesRaw <- read.table(here("assets/Trawls_2017.txt"), header=TRUE, fill=TRUE, sep="\t", check.names=FALSE)

trawls_raw <- specaccum(SpeciesRaw, "random")
summary(trawls_raw)
plot(trawls_raw, ci.type="poly", col="black", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 210), ylim=c(0,30), bg = "transparent")

#Finally, a heatmap of read abundance, based on hellinger transformed data

Elasmo_Summary <- Elasmo_Hell_zeros  %>%
  group_by(Year) %>%
  group_by(Month, .add = TRUE) %>%
  summarise(across(Alopias_vulpinus:Squalus_acanthias, mean, na.rm= TRUE))

Elasmo_Summary$Year <- as.factor(Elasmo_Summary$Year)
Elasmo_Summary$Month <- as.factor(Elasmo_Summary$Month)
Elasmo_Summary_long <- melt(Elasmo_Summary, id.vars = c("Month", "Year"))

Elasmo_Summary_long$Year_Month = paste(Elasmo_Summary_long$Year,Elasmo_Summary_long$Month)
Elasmo_Summary_long$Year_Month <- as.factor(Elasmo_Summary_long$Year_Month)

library(plyr)

revalue(Elasmo_Summary_long$Year_Month, c("2017 2" = "2017_02")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2017 3" = "2017_03")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2017 4" = "2017_04")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2017 5" = "2017_05")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2017 6" = "2017_06")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2017 7" = "2017_07")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2017 8" = "2017_08")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2017 9" = "2017_09")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2017 10" = "2017_10")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2017 11" = "2017_11")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2018 2" = "2018_02")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2018 3" = "2018_03")) -> Elasmo_Summary_long$Year_Month
revalue(Elasmo_Summary_long$Year_Month, c("2018 4" = "2018_04")) -> Elasmo_Summary_long$Year_Month

#Need to load in file of names
species <- read.table(here::here("assets/species.txt"), header=TRUE, fill=TRUE, sep="\t", check.names=FALSE)

Elasmo_Summary_long <- cbind(Elasmo_Summary_long, species)

Elasmo_Summary_long$species <- forcats::fct_rev(factor(Elasmo_Summary_long$species))

#Generate the plot
Heatmap <-  ggplot(Elasmo_Summary_long,aes(y=Year_Month,x=species,fill=value))+
  scale_y_discrete(limits=c("2017_02","2017_03","2017_04","2017_05","2017_06","2017_07","2017_08","2017_09","2017_10","2017_11","2018_02","2018_03","2018_04")) +
  geom_tile(colour="grey",size=0.25)+
  labs(x="",y="")+ 
  scale_fill_gradient(low = "white", high = "red", na.value = NA)+
  theme_grey(base_size=8)+ coord_flip() +
  theme(
    legend.text=element_text(face="bold"),
    axis.ticks=element_line(size=0.3),
    axis.text.x=element_text(angle = -45, hjust = 0, size =8),
    axis.text.y=element_text(size =8),
    plot.background=element_blank(),
    panel.border=element_blank())
Heatmap

#Analyses done
