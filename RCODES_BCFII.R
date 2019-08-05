##### BALD CYPRESS SEEDLING FLOODING EXPERIMENT DATA ANALYSES#######

#* Setting your working directory
setwd("~/Desktop/GrowthAnalysis_BCFII/GrowthAnalyses") # This code sets your working directory

############################################################################################# 
## PACKAGES AND LIBRARIES NEEDED 
#############################################################################################

install.packages("Rmisc")
install.packages("ggplot2", dependencies=TRUE)
install.packages ("grid")
install.packages ("gridExtra")
install.packages ("grid")
install.packages("MASS")
install.packages("lme4")
install.packages("agridat")
install.packages("car")
install.packages("lsmeans")
install.packages("lmerTest")
install.packages("lme4")
install.packages("agricolae")
install.packages( "multcomp")
install.packages ("emmeans")
install.packages("pls")

library(ggplot2)
library(gridExtra)
library(grid)
install.packages("ggpubr")
library(ggpubr)
install.packages("Rmisc")
library(Rmisc)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library (car) ## Needed for testing homogeneity of variances
library(agricolae) # Needed to obtain Tukey groupings
library(MASS)
library(lme4)
library(lsmeans)
library(lmerTest)
library(foreign)
library(multcomp)
library(multcompView)
library(pls) ## Needed for partial least-square analyses

############################################################################################# 
## DATA SETS TO BE USED:
#############################################################################################

#****** Growth after flooding - Total differences between initial and final growth measurements

GAF <- read.table("HeightRelChange_Final.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 

#****** Plant performance i.e. dry biomass #########

P <- read.table("PlantPerformance_Allvariables.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 
P$DryBiomassRatio <- P$ADW/P$BDW  ### Calculating Shoot (Above) to root (belowground) ratios
P$TDB <- P$ADW + P$BDW  #### Calculating Total Dry Biomass (TDB)
P$RADW <- P$ADW/P$IH    #### Calcutlating Relative Aboveground dry biomass with respect to the initial height of plant before experiment (H)   
P$RBDW <- P$BDW/P$IH    #### Calculating Relative Belowground dry biomass
P$RatioRAB <- P$RADW/P$RBDW  ### Calculating Relative Above-Belowground biomass ratio

#****** FUNGI (AMF & DSE) #########

E <- read.table("AMF_datamatrix_noreps.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 

E2 <- read.table("AMF_datamatrix_noreps_nooutlier.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 


############################################################################
## TESTING ASSUMPTIONS OF NORMALITY & HETEROSCEDASCITY
############################################################################
install.packages("car")  # For testing homogeneity of variances
library (car)

#****************************************************************************
#             GROWTH VARIABLES: CHANGE IN HEIGHT, SLA
#****************************************************************************
# Relative Total Change in Height (RTHO)
hist(GAF$RTHO, probability = T, main= "Histogram of normal data - Relative Total Change in Height", xlab="Approximately normally distributed data")
lines(density(GAF$RTHO), col=2)
qqnorm(GAF$RTHO,main="QQ plot of normal data- Relative Total Change in Height",pch=19) 
qqline(GAF$RTHO) 

shapiro.test(GAF$RTHO) # To get shapiro-wilk tests

leveneTest(RTHO ~ Soil *Hydrology, data=GAF) # to get levene's tests of homogeneity of variances
            

            #* Transforming data & testing assumptions using Linear Mixed Model:
            GAF$LogRTHO <- log(GAF$RTHO + 1)
            lmmRTHO <- lmer(LogRTHO ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                           REML = FALSE)
            hist(resid(lmmRTHO,type="deviance"))
            qqnorm(resid(lmmRTHO,type="deviance"))
            qqline(resid(lmmRTHO,type="deviance"))
            
            Anova(lmmRTHO)
            #* Conclusion: residuals are not perfect but very close to normal. Proceed with this model!

#SLAFL - specific leave area of firts leave
hist(P$SLAFL, probability = T, main= "Histogram of normal data - SLAFL", xlab="Approximately normally distributed data")
lines(density(P$SLAFL), col=2)
qqnorm(P$SLAFL,main="QQ plot of normal data- SLAFL",pch=19) 
qqline(P$SLAFL) 
shapiro.test(P$SLAFL) # Test of Normality - Not normality in this data
leveneTest(SLAFL ~ Soil*Inundation, data=P) ## There is homogeneity of variance
plot(SLAFL ~ Soil*Inundation, data=P)
            
            #* Transforming data & testing assumptions using Linear Mixed Model:
            P$LogSLA <- log(P$SLAFL)
            lmmSLA <- lmer(LogSLA ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = FALSE)
            hist(resid(lmmSLA,type="deviance"))
            qqnorm(resid(lmmSLA,type="deviance"))
            qqline(resid(lmmSLA,type="deviance"))
            #* Conclusion: residuals are normal. Proceeded with LM!
            
            
# RHOAF = Change in Height after flooding
hist(GAF$RHOAF, probability = T, main= "Histogram of normal data - Relative Total Change in Height", xlab="Approximately normally distributed data")
lines(density(GAF$RHOAF), col=2)
qqnorm(GAF$RHOAF,main="QQ plot of normal data- Relative Total Change in Height",pch=19) 
qqline(GAF$RHOAF) 
shapiro.test(GAF$RHOAF)
leveneTest(RHOAF ~ Soil*Inundation, data=GAF) # Testing heteroscesdacity

          #* Transforming data & testing assumptions using Linear Mixed Model:
          GAF$LogRHOAF <- log(GAF$RHOAF + 1)
          lmmRHOAF <- lmer(LogRHOAF ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                REML = FALSE)
          hist(resid(lmmRHOAF,type="deviance"))
          qqnorm(resid(lmmRHOAF,type="deviance"))
          qqline(resid(lmmRHOAF,type="deviance"))
          
          #* Conclusion: residuals are not normal. Proceeded with GLMMS!


# RHIN = Change in height after inoculation
hist(GAF$RHIN, probability = T, main= "Histogram of normal data - Relative Total Change in Height", xlab="Approximately normally distributed data")
lines(density(GAF$RHIN), col=2)
qqnorm(GAF$RHIN,main="QQ plot of normal data- Relative Total Change in Height",pch=19) 
qqline(GAF$RHIN) 
shapiro.test(GAF$RHIN)
leveneTest(RHIN ~ Soil*Inundation, data=GAF)

          #* Transforming data & testing assumptions using Linear Mixed Model:
          GAF$LogRHIN <- log(GAF$RHIN + 1)
          lmmRHIN <- lmer(LogRHIN ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                 REML = FALSE)
          hist(resid(lmmRHIN,type="deviance"))
          qqnorm(resid(lmmRHIN,type="deviance"))
          qqline(resid(lmmRHIN,type="deviance"))
          #* Conclusion: residuals are normal. Proceeded with LM!


#IH = Initial Height of plants
hist(P$IH, probability = T, main= "Histogram of normal data - Leave Dry Biomass", xlab="Approximately normally distributed data")
lines(density(P$IH), col=2)
qqnorm(P$IH,main="QQ plot of normal data- Leave Dry Biomass",pch=19) 
qqline(P$IH) 
shapiro.test(P$IH) # Test of Normality: met assumptions
leveneTest(IH ~ Microbe*Inundation, data=P)

#FH = Final Height of plants
hist(P$IH, probability = T, main= "Histogram of normal data - Leave Dry Biomass", xlab="Approximately normally distributed data")
lines(density(P$FH), col=2)
qqnorm(P$FH,main="QQ plot of normal data- Leave Dry Biomass",pch=19) 
qqline(P$FH) 
shapiro.test(P$FH) # Test of Normality : met assumptions
leveneTest(FH ~ Microbe*Inundation, data=P) # 

#****************************************************************************
#            PERFORMANCE VARIABLES: DRY BIOMASS MEASUREMENTS
#****************************************************************************

# Aboveground biomass (AB)
hist(P$ADW, probability = T, main= "Histogram of normal data - Aboveground Dry Biomass", xlab="Approximately normally distributed data")
lines(density(P$SDW), col=2)
qqnorm(P$ADW,main="QQ plot of normal data- Aboveground Dry Biomass",pch=19) 
qqline(P$ADW) 
shapiro.test(P$ADW) # Test of Normality
leveneTest(ADW ~ Microbe*Inundation, data=P)

          #* Transforming data & testing assumptions using Linear Mixed Model:
          P$LogADW <- log(P$ADW)
          lmmADW <- lmer(LogADW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                REML = FALSE)
          hist(resid(lmmADW,type="deviance"))
          qqnorm(resid(lmmADW,type="deviance"))
          qqline(resid(lmmADW,type="deviance"))
          #* Conclusion: residuals are normal. Proceeded with LM!


# Belowground biomass (BB)
hist(P$BDW, probability = T, main= "Histogram of normal data - Belowground Dry Biomass", xlab="Approximately normally distributed data")
lines(density(P$BDW), col=2)
qqnorm(P$BDW,main="QQ plot of normal data- Belowground Dry Biomass",pch=19) 
qqline(P$BDW) 
shapiro.test(P$BDW) # Test of Normality
leveneTest(BDW ~ Microbe*Inundation, data=P)

          #* Transforming data & testing assumptions using Linear Mixed Model:
          P$LogBDW <- log(P$BDW)
          lmmBDW <- lmer(LogBDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
               REML = TRUE)
          hist(resid(lmmBDW,type="deviance"))
          qqnorm(resid(lmmBDW,type="deviance"))
          qqline(resid(lmmBDW,type="deviance"))
          #* Conclusion: residuals are normal. Proceeded with LM!
          
          Anova(lmmBDW)



# Total Dry Biomass (TDB)
hist(P$TDB, probability = T, main= "Histogram of normal data - Shoot to root ratio", xlab="Approximately normally distributed data")
lines(density(P$TDB), col=2)
qqnorm(P$TDB,main="QQ plot of normal data- Shoot to root ratio",pch=19) 
qqline(P$TDB) 
shapiro.test(P$TDB) # Test of Normality
leveneTest(TDB ~ Microbe*Inundation, data=P)

        #* Transforming data & testing assumptions using Linear Mixed Model:
        P$LogTDB <- log(P$TDB)
        lmmTDB <- lmer(LogTDB ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                   REML = FALSE)
        hist(resid(lmmTDB,type="deviance"))
        qqnorm(resid(lmmTDB,type="deviance"))
        qqline(resid(lmmTDB,type="deviance"))
        #* Conclusion: residuals are normal. Proceeded with LM!

# DryBiomassRatio (Shoot:Root ratio)
hist(P$DryBiomassRatio, probability = T, main= "Histogram of normal data - Shoot to root ratio", xlab="Approximately normally distributed data")
lines(density(P$BDW), col=2)
qqnorm(P$DryBiomassRatio,main="QQ plot of normal data- Shoot to root ratio",pch=19) 
qqline(P$DryBiomassRatio) 
shapiro.test(P$DryBiomassRatio) # Test of Normality
leveneTest(DryBiomassRatio ~ Microbe*Inundation, data=P)
        
        #* Transforming data & testing assumptions using Linear Mixed Model:
        P$LogABratio <- log(P$DryBiomassRatio)
        lmmABratio <- lmer(LogABratio ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = FALSE)
        hist(resid(lmmABratio,type="deviance"))
        qqnorm(resid(lmmABratio,type="deviance"))
        qqline(resid(lmmABratio,type="deviance"))
        #* Conclusion: residuals are normal. Proceeded with LM!
        
        Anova(lmmABratio)
        
# Leaf biomass
hist(P$LDW, probability = T, main= "Histogram of normal data - Shoot to root ratio", xlab="Approximately normally distributed data")
lines(density(P$LDW), col=2)
qqnorm(P$LDW,main="QQ plot of normal data- Shoot to root ratio",pch=19) 
qqline(P$LDW) 
shapiro.test(P$LDW) # Test of Normality
leveneTest(LDW ~ Soil*Inundation, data=P)
        
        #* Transforming data & testing assumptions using Linear Mixed Model: 
        P$LogLDW <- log(P$LDW)
        lmmLDW <- lmer(LogLDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = TRUE)
        hist(resid(lmmLDW,type="deviance"))
        qqnorm(resid(lmmLDW,type="deviance"))
        qqline(resid(lmmLDW,type="deviance"))
        #* Conclusion: residuals are normal. Proceeded with LM!
        
        
# Stem biomass
hist(P$SDW, probability = T, main= "Histogram of normal data - Shoot to root ratio", xlab="Approximately normally distributed data")
lines(density(P$SDW), col=2)
qqnorm(P$SDW,main="QQ plot of normal data- Shoot to root ratio",pch=19) 
qqline(P$SDW) 
shapiro.test(P$SDW) # Test of Normality
leveneTest(SDW ~ Soil*Inundation, data=P)
        
        #* Transforming data & testing assumptions using Linear Mixed Model: 
        P$LogSDW <- log(P$SDW)
        lmmSDW <- lmer(LogSDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                       REML = TRUE)
        hist(resid(lmmSDW,type="deviance"))
        qqnorm(resid(lmmSDW,type="deviance"))
        qqline(resid(lmmSDW,type="deviance"))
        #* Conclusion: residuals are normal. Proceeded with LM!
        

        
############################################################################################################
# Running Linear Mixed Models & extracting LSmeans
#############################################################################################################

options("contrasts")  # To veryfy which type of contrast is R using for the ANOVAs tables
options(contrasts=c("contr.helmert","contr.poly")) ### To set up correct contrasts for Type III ANOVA's
#options(contrasts=c("contr.treatment","contr.poly"))
#To comeback to default options: backup_options <- options() and then type options(backup_options). Do this before changing any options when you start R

#RTHO= Relative Total Change in Height
lmmRTHO <- lmer(LogRTHO ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                REML = TRUE) #Linear Mixed Model
Anova(lmmRTHO, type = "III")  # These are the results to report

          #*Posthoc tests & LSmeans
          library(multcomp)
          library(multcompView)
          lsmeans(lmmRTHO, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
          lsmeans(lmmRTHO, pairwise ~ Inundation|Soil, adjust ="bonferroni") #within soil treatments among hydrological regimes

            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            RTHO.lsm <- lsmeans(lmmRTHO, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(RTHO.lsm, alpha = .50)
            
            #*Testing random factor effect
            lmmRTHO1 <- lmer(LogRTHO ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                            REML = FALSE)
            lmmRTHO0 <- lm(LogRTHO ~ Inundation + Soil + Soil*Inundation, data = GAF)
            anova(lmmRTHO1, lmmRTHO0)
            
            
#SLA = Specific Leaf Area
P$LogSLA <- log(P$SLAFL)
lmmSLA <- lmer(LogSLA ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = FALSE)
Anova(lmmSLA, type = "III")

            #*Posthoc tests & LSmeans
            library(multcomp)
            library(multcompView)
            lsmeans(lmmSLA, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            lsmeans(lmmSLA, pairwise ~ Inundation|Soil, adjust ="bonferroni") #within soil treatments among hydrological regimes
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            SLA.lsm <- lsmeans(lmmSLA, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(SLA.lsm, alpha = .50)
            
            #*Testing random factor effect
            lmmSLA1 <- lmer(LogSLA ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                             REML = FALSE)
            lmmSLA0 <- lm(LogSLA ~ Inundation + Soil + Soil*Inundation, data = P)
            anova(lmmSLA1, lmmSLA0)
            
            
            
            

#TDB = Total Dry Biomass
P$LogTDB <- log(P$TDB)

lmmTDB <- lmer(LogTDB ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P, REML = FALSE)   

            #*Posthoc tests & LSmeans
            library(multcomp)
            library(multcompView)
            lsmeans(lmmTDB, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            lsmeans(lmmTDB, pairwise ~ Inundation|Soil, adjust ="bonferroni") #within soil treatments among hydrological regimes
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            TDB.lsm <- lsmeans(lmmTDB, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(TDB.lsm, alpha = .50)
            
            #*Testing random factor effect
            lmmTDB1 <- lmer(LogTDB ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                            REML = FALSE)
            lmmTDB0 <- lm(LogTDB ~ Inundation + Soil + Soil*Inundation, data = P)
            anova(lmmTDB1, lmmTDB0)
            

#Shoot:Root ratio
            
P$LogABratio <- log(P$DryBiomassRatio)
lmmABratio <- lmer(LogABratio ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = FALSE)   
            
            #*Posthoc tests & LSmeans
            library(multcomp)
            library(multcompView)
            lsmeans(lmmABratio, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            lsmeans(lmmABratio, pairwise ~ Inundation|Soil, adjust ="bonferroni") #within soil treatments among hydrological regimes
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            ABratio.lsm <- lsmeans(lmmABratio, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(ABratio.lsm, alpha = .50)
            
            #*Testing random factor effect
            lmmABratio1 <- lmer(LogABratio ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                            REML = FALSE)
            lmmABratio0 <- lm(LogABratio ~ Inundation + Soil + Soil*Inundation, data = P)
            anova(lmmABratio1, lmmABratio0)
            

#########################################################################################
            ##  ANALYSES OF AMF & DSE USING GLMS
##########################################################################################
            
E <- read.table("AMF_datamatrix_noreps.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 
            
            
mean(E$TAMF2)
var(E$TAMF2)
            
mean(E$TDSE)
var(E$TDSE)
            
table(E$TAMF)
table(E$TDSE)

            
#************* For total AMF counts - Data was categorized based on counts (e.g. TAMF2)
            
#Generalized Mixed model w Seed as random effect
            AMFfit_full <-  glmer(TAMF2 ~ Soil + Inundation + Soil*Inundation + (1|Seed), family=poisson(link="log"), data = E) #* Full possible model
            summary(AMFfit_full)
            Anova(AMFfit_full, type = "III")
            
            #*Testing significance of each factor
            AMFfit0 <-  glmer(TAMF2 ~ Soil + (1|Seed), family=poisson(link="log"), data = E) 
            AMFfit1 <-  glmer(TAMF2 ~ Soil + Inundation + (1|Seed), family=poisson(link="log"), data = E) 
            AMFfit2 <-  glmer(TAMF2 ~ Soil + Inundation + Soil*Inundation + (1|Seed), family=poisson(link="log"), data = E) 
            anova(AMFfit0,AMFfit1,AMFfit2)
            
            #*Testing significance of random factor
            AMFfit_nram <-  glm(TAMF2 ~ Soil + Inundation + Soil*Inundation, family=poisson(link="log"), data = E)
            anova(AMFfit2,AMFfit_nram,test="Chisq")
            
            #*Posthoc tests
            lsmeans(AMFfit2, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            lsmeans(AMFfit2, pairwise ~ Inundation|Soil, adjust ="bonferroni") #within soil treatments among hydrological regimes
            
            #*Tukey-groupings
            AMFfit_full <-  glmer(TAMF2 ~ Soil + Inundation + Soil*Inundation + (1|Seed), family=poisson(link="log"), data = E) #* Full possible model
            summary(AMFfit_full)
            AMF.lsm <- lsmeans(AMFfit_full, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(AMF.lsm, alpha = .50)
            
#******************* For total DSE counts - Data categorized into TDSE2
            
            #Generalized Mixed model w Seed as random effect
            DSEfit_full <-  glmer(TDSE2 ~ Soil + Inundation + Soil*Inundation + (1|Seed), family=poisson(link="log"), data = E) #* Full possible model
            summary(DSEfit_full)
            Anova(DSEfit_full)
            
            #*Testing significance of each factor
            DSEfit0 <-  glmer(TDSE2 ~ Soil + (1|Seed), family=poisson(link="log"), data = E) 
            DSEfit1 <-  glmer(TDSE2 ~ Soil + Inundation + (1|Seed), family=poisson(link="log"), data = E) 
            DSEfit2 <-  glmer(TDSE2 ~ Soil + Inundation + Soil*Inundation + (1|Seed), family=poisson(link="log"), data = E) 
            anova(DSEfit0,DSEfit1,DSEfit2)
            
            #*Testing significance of random factor
            DSEfit_nram <-  glm(TDSE2 ~ Soil + Inundation + Soil*Inundation, family=poisson(link="log"), data = E)
            anova(DSEfit2,DSEfit_nram,test="Chisq")
            
            #*Posthoc tests
            lsmeans(DSEfit2, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            lsmeans(DSEfit2, pairwise ~ Inundation|Soil, adjust ="bonferroni") #within soil treatments among hydrological regimes
            
            
#******************************************************************************************************************************************************* 
            #*        Analyses in total % fungal colonization for non-sterile treatment only
#*******************************************************************************************************************************************************
            
            #Dataset:
            RootFungi <- read.csv("RootFungi.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 
            
            
            #*Testing for Normality
            hist(RootFungi$Colonization, probability = T, main= "Histogram of normal data - SLALL", xlab="Approximately normally distributed data")
            lines(density(RootFungi$Colonization), col=2)
            qqnorm(RootFungi$Colonization,main="QQ plot of normal data- SLAFL",pch=19) 
            qqline(RootFungi$Colonization) 
            shapiro.test(RootFungi$Colonization) # Test of Normality Not met (p= 0.04 ; so not too bad)
            leveneTest(Colonization ~ Inundation*Fungi, data=RootFungi) ## homogeneity of variance was met
            
            #transforming data using logit
            RootFungi$logitColo <- logit(RootFungi$Colonization,adjust = 0.025)
            
            #testing normality of logit-transformed data
            hist(RootFungi$logitColo, probability = T, main= "Histogram of normal data - SLALL", xlab="Approximately normally distributed data")
            lines(density(RootFungi$logitColo), col=2)
            qqnorm(RootFungi$logitColo,main="QQ plot of normal data- SLAFL",pch=19) 
            qqline(RootFungi$logitColo) 
            shapiro.test(RootFungi$logitColo) # Test of Normality - Transformed data did not met assumptions of nornality; it made it worse
            leveneTest(logitColo ~ Inundation, data=RootFungi) ## There is not homogeneity of variance
            
            #transforming data using arcsine
            RootFungi$arcColo <- asin(sqrt(RootFungi$Colonization))
            
            #testing normality of arcsine-transformed data
            hist(RootFungi$arcColo, probability = T, main= "Histogram of normal data - SLALL", xlab="Approximately normally distributed data")
            lines(density(RootFungi$arcColo), col=2)
            qqnorm(RootFungi$arcColo,main="QQ plot of normal data- SLAFL",pch=19) 
            qqline(RootFungi$arcColo) 
            shapiro.test(RootFungi$arcColo) # Test of Normality - Transformed data did not met assumptions of nornality
            leveneTest(arcColo ~ Inundation, data=RootFungi) ## There is not homogeneity of variance
            #* Normality and homogeneity were made worse by the transformations - Using GLMMS next
            
            
#******************************* GLIM Poisson model with categorical data
            
            sumColonization = summarySE(data=P, measurevar = "TDB", groupvars = c("Soil", "Hydrology"))
            #Colcat is percent colonization in categories
            
            #Generalized Mixed model w Seed as random effect
            Fungifit_full <-  glmer(ColCat2 ~ Fungi + Inundation + Fungi*Inundation + (1|Seed), family=poisson(link="log"), data = RootFungi) #* Full possible model
            summary(Fungifit_full)
            Anova(Fungifit_full)
            
            #*Testing significance of each factor
            DSEfit0 <-  glmer(ColCat2 ~ Fungi + (1|Seed), family=poisson(link="log"), data = E) 
            DSEfit1 <-  glmer(ColCat2 ~ Fungi + Inundation + (1|Seed), family=poisson(link="log"), data = E) 
            DSEfit2 <-  glmer(TDSE2 ~ Soil + Inundation + Soil*Inundation + (1|Seed), family=poisson(link="log"), data = E) 
            anova(DSEfit0,DSEfit1,DSEfit2)
            
            #*Testing significance of random factor
            DSEfit_nram <-  glm(TDSE2 ~ Soil + Inundation + Soil*Inundation, family=poisson(link="log"), data = E)
            anova(DSEfit2,DSEfit_nram,test="Chisq")
            
            #*Posthoc tests
            lsmeans(DSEfit2, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            lsmeans(DSEfit2, pairwise ~ Inundation|Soil, adjust ="bonferroni") #within soil treatments among hydrological regimes
            
##################################################################################################################################################################### 
            #######     PARTIAL-LEAST SQUARES ANALYSES ####### 
#####################################################################################################################################################################
            #Data sets:
            
            ##this will download and install the pls package
            install.packages("pls")
            
            ##this will load the pls package
            library(pls)
            
            ##this will read the csv file with your data
            df <- read.csv("AMF_datamatrix_noreps-1.csv")
            
            #PLS CODE
            ###scale = TRUE; means that each variable is standardised by dividing it by its standard deviation, and if scale is a numeric vector, each variable is divided by the corresponding number.
            ##ncomp = #; means number of components to use. Leaving this out of the code means that PLS will run the max number of components.
            ##validation = "CV"; means that the program will perform Cross-Validation while "LOO" will perform Leave One Out cross validation
            ##method = this just selects the algorithm that pls will use. There are two more algorithms in the pls package but this one worked best for my data
            
            AMF_IN <- subset(df, Soil == "IN", select=ID:FH) ## Extracting only the data for the IN treatment
            
            df.metric <- plsr(TDB  ~ FUNGI + SLA + ADW + BDW + RTHO, scale = TRUE, validation = "LOO", method = "oscorespls", data = df)
            
            # If wanting to test plant growth and performance
            #df.metric <- plsr(TDB + RTHO ~ TAMF + TDSE + SLAFL + ADW + BDW, scale = TRUE, validation = "CV", method = "oscorespls", data = df)
            
            # For testing associations in seedlings reared in non-sterile soils
            df.metric <- plsr(TDB  ~ TAMF + TDSE + SLAFL + ADW + BDW + RTHO, scale = TRUE, validation = "LOO", method = "oscorespls", data = AMF_IN)
            
            
            ##summary results
            summary(df.metric)
            
            
            #prediction plot:
            #This shows the cross-validated predictions with two components versus measured values
            #We have chosen an aspect ratio of 1, and to draw a target line.  The points follow
            #the target line quite nicely, and there is no indication of a curvature or other anomalies
            
            plot(df.metric, ncomp=2, asp=1, line=TRUE)
            
            ##explained variance:  relative amount of X variance explained by each componen
            exp.var <- explvar(df.metric)
            exp.var
            
        
            ##plot Root Mean Squared Error Prediction (RMSEP) - to see how many components are needed to fit the model
            plot(RMSEP(df.metric), legendpos = "topright")
            
            # To test how many components are appropiate for fitting the model:
            ncomp.onesigma <- selectNcomp(df.metric, method="onesigma", plot=TRUE) # This develops the onesigma procedure
            ncomp.permut <- selectNcomp(df.metric, method="randomization", plot=TRUE) # This develops the a permutation approach
                                  # According to these results, 6 components are needed for the best fit of the model
            
            ###* Plotting
            par(c(2,2))
            plot(RMSEP(df.metric), legendpos = "topright")
            ncomp.onesigma 
            
            ## correlation loadings plot, also known as loading weights 
            plsIN <- plot(df.metric, "correlation", comps = 1:2,labels ="names", cex = 0.8)
            
            ##biplot <- the biplot is a summary map including the X- and Y-variables as well as the observations
            biplot(df.metric, comps = 1:2, which = "loadings", cex = 0.8, expand = .7)
            
          
            ##prediction plot
            plot(df.metric, ncomp = 5, asp = 1, line = TRUE)
            
            ## scores plot
            plot(df.metric, plottype = "scores", comps = 1:2)
            
            ## correlation loadings plot, also known as loading weights
            plot(df.metric, "correlation", comps = 1:2,labels ="names")
            abline(h=0)
            
            plot(df.metric, "correlation", comps = 1:2,type = "points")
            abline(h=0)
            
            ##loadings 
            loadings(df.metric)
            
            ##correlation loadings plot
            corrplot(df.metric, comps = 1:2,labels = "names")
            
            
            
            ##Regression coefficients
            coef(df.metric)
            
            plot(df.metric, "biplot")
            
            #correlation plot
            corrplot(df.metric, loadings, comps = 1:2, radii = c(sqrt(1/2), 1), identify = TRUE, labels = "names")
            
            #scoreplot
            plot(df.metric, plottype = "loadings", comps = 1:2, labels = "names" )
            
            #the explained variances
            explvar(df.metric)
            
            lw <- loading.weights(df.metric)
            lw
            
            ########VIP
            ## VIP returns all Variable Iimportance to the Projections values for all variables and all number of components,
            ## as a ncomp x nvars matrix.
            
            VIP <- function(df.metric) {
              if (df.metric$method != "oscorespls")
                stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
              if (nrow(df.metric$Yloadings) > 1)
                stop("Only implemented for single-response models")
              
              SS <- c(df.metric$Yloadings)^2 * colSums(df.metric$scores^2)
              Wnorm2 <- colSums(df.metric$loading.weights^2)
              SSW <- sweep(df.metric$loading.weights^2, 2, SS / Wnorm2, "*")
              sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
            }
            
            VIP<-VIP(df.metric)
            VIP
            
###############################################################################
            ###### MAIN MANUSCRIPT FIGURES
################################################################################
            
            # To change order of Inundation levels
            GAF$Hydrology <-factor(GAF$Inundation,levels=levels(GAF$Inundation)[c(2,1,3)])
            P$Hydrology <-factor(P$Inundation,levels=levels(P$Inundation)[c(2,1,3)])
            
            
            ## RTHO PLOT
            
            #Need to recalculate lsmeans
            glmm_RTHO <- glmer(RTHO.t ~ Hydrology + Soil + Hydrology*Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            library(lsmeans)
            lmmRTHO.rg1 <- ref.grid(glmm_RTHO, type="response") # here we type = "response" to obtain back-transformed lsmeans
            lsmeansRTHO <- summary(lmmRTHO.rg1)
            
            # Making plot in ggplot2
            sumRTHOplot <- ggplot(lsmeansRTHO, aes(x= Soil, y= response, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "RTHO (cm)") 
            p5 <- sumRTHOplot + theme_bw()
            #p5.1 <- p5 + theme(legend.position="none")
            p5.3 <- p5 + ylim(0,5)
            p5.2 <- p5.3 + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                 axis.text.x = element_blank(),
                                 axis.text.y = element_text(size=14, colour ="black"))
            pRTHO <- p5.2 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                  panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            ggsave("pRTHO.pdf", plot=pRTHO, device="pdf", scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
            
            
            
            ########## SLA plot
            
            #Need to recalculate lsmeans
            lmmSLAFL <- lmer(SLAFL ~ Hydrology + Soil + Soil*Hydrology +  (1 | Seed), data = P,
                             REML = FALSE) # Full model for SLAFL
            library(lsmeans)
            lmmSLA.rg1 <- ref.grid(lmmSLAFL) 
            lsmeansSLA <- summary(lmmSLA.rg1)
            
            # plot
            SLAFLplot <- ggplot(lsmeansSLA, aes(x= Soil, y= prediction, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=prediction-SE, ymax=prediction+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y="SLA (m2 Kg-1)") 
            p11 <- SLAFLplot + theme_bw() 
            p12 <- p11 + ylim(0,2.5)
            p13 <- p12 + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                               axis.title.x = element_blank(),
                               axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                               axis.text.x = element_blank(),
                               axis.text.y = element_text(size=14, colour ="black"))
            pSLAFLplot <- p13 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                      panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            ggsave("pSLAFLplot.pdf", plot=pSLAFLplot, device="pdf", scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
            
            
            #### TOTAL DRY BIOMASS PLOT (TDB)
            
            #Need to recalculate lsmeans
            glmm_TDB2 <- glmer(TDB ~  Soil + Hydrology + Soil*Hydrology + (1 | Seed), data = P, family = gaussian(link = "log")) 
            lmmTDB.rg1 <- ref.grid(glmm_TDB2, type="response") # here we type = "response" to obtain back-transformed lsmeans
            lsmeansTDB <- summary(lmmTDB.rg1)
            
            # Plot TDB
            TDBplot <- ggplot(lsmeansTDB, aes(x= Soil, y= response, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "Total Dry Biomass (gr)") 
            TDBplot1 <- TDBplot + ylim(0,3)
            pTDB <- TDBplot1 + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                     axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                     axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                     axis.text.x = element_text(size=14, colour ="black"),
                                     axis.text.y = element_text(size=14, colour ="black"))
            p6 <- pTDB + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pTDW <- p6 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                               panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            ggsave("pTDW.pdf", plot=pTDW, device="pdf", scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
            
            ################ PLOT A:B ratio 
            
            #Need to recalculate lsmeans
            library(lsmeans)
            glmm_ABratio <- glmer(DryBiomassRatio ~ Hydrology + Soil + Hydrology*Soil + (1 | Seed), data = P, family = gaussian(link = "log"))  
            lmmABratio.rg1 <- ref.grid(glmm_ABratio, type="response") # here we type = "response" to obtain back-transformed lsmeans
            lsmeansABratio <- summary(lmmABratio.rg1)
            
            #Plotting A:B ratio
            
            DryBiomassRatioplot <- ggplot(lsmeansABratio, aes(x= Soil, y= response, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "A:B ratio (gr)") 
            p8 <- DryBiomassRatioplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                              axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                              axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                              axis.text.x = element_text(size=14, colour ="black"),
                                              axis.text.y = element_text(size=14, colour ="black"))
            p9 <- p8 + ylim(0,3)
            p10 <- p9 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pRatio <- p10 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                  panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            # **Plotting all graphs together with a shared legend  
            Growth <- ggarrange(pIH, pFH, pRTHO, ncol=3, common.legend = TRUE, legend="top")
            #grid.arrange(pIH, pFH, pRTHO, pSLAFLplot, ncol=4)
            ggsave("Growth.pdf", plot=Growth, device="pdf", scale = 1, width = 20, height = 10, units = "cm", dpi = 300)
            
            Figure4 <- ggarrange(pRTHO,pSLAFLplot,pTDW, pRatio, ncol=2,nrow=2, common.legend = TRUE, legend="top")
            ggsave("Figure3.pdf", plot=Figure3, device="pdf", scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
            
            
            
##### TDSE & TAMF for Inoculant only
            
            RootFungi <- read.csv("RootFungi.csv", header = TRUE, sep=",", strip.white = TRUE, na.strings = "EMPTY") 
            RootFungi$Hydrology <- factor(RootFungi$Inundation,levels=levels(RootFungi$Inundation)[c(2,1,3)]) # To arrange inundation at specific order
            sumRootFungi = summarySE(data=RootFungi, measurevar= "Colonization", groupvars = c("Fungi", "Hydrology"))
            sumRootFungiType = summarySE(data=RootFungi, measurevar= "Colonization", groupvars = c("Fungi"))
            
            
            
            Fungiplot <- ggplot(sumRootFungi, aes(x= Fungi, y= Colonization, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=Colonization-se, ymax=Colonization+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Type of Fungi", y= "Root Colonization (%)") 
            p8 <- Fungiplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                    axis.text.x = element_text(size=14, colour ="black"),
                                    axis.text.y = element_text(size=14, colour ="black"),
                                    legend.position = "top")
            p9 <- p8 + ylim(0,0.45)
            p10 <- p9 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pFungi <- p10 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                  panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            Figure3 <- ggarrange(pFungi, ncol=1, common.legend = TRUE, legend="top")
            ggsave("Figure3_RootFungi.pdf", plot=pFungi, device="pdf", scale = 1, width = 20, height = 20, units = "cm", dpi = 300)
            
            
            
            
            
            
                     
###############################################################################
            ###### SUPPLEMENTAL MATERIALS ANALYSES
################################################################################          
           
            
            head(GAF)
#IH - models
            lmmIH <- lmer(IH ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                          REML = FALSE) # Simple model for IH
            summary(lmmIH)
            Anova(lmmIH)
            
            #* testing the effect adding inundation to the model
            lmmIH.null <- lmer(IH ~ Soil +  (1 | Seed), data = GAF,
                               REML = FALSE) # Simple model for IH
            
            lmmIH.inu <- lmer(IH ~ Soil + Inundation +  (1 | Seed), data = GAF,
                              REML = FALSE) # Adding inundation
            
            anova (lmmIH.null,lmmIH.inu) # Likelihood ratio test to demonstrate no effect of inundation
            
            #* testing the effect adding soil to the model
            lmmIH.0 <- lmer(IH ~ 1 + (1 | Seed), data = GAF,
                            REML = FALSE) # Simple model for IH
            
            lmmIH.soil <- lmer(IH ~ Soil + (1 | Seed), data = GAF,
                               REML = FALSE) # Adding inundation
            
            anova (lmmIH.0,lmmIH.soil)
            
            #* testing the effect of random factors in the model
            
            lmmIH.0 <- lmer(IH ~ Soil + Hydrology + (1 | Seed), data = GAF,
                            REML = FALSE) # Simple model for IH
            
            lmmIH.1 <- lmer(IH ~ Soil + Hydrology + , data = GAF,
                            REML = FALSE)
            
            anova(lmmIH.1,lmmIH.0)
            
            
            # lsmeans IH
            library(lsmeans)
            lmmIH.rg1 <- ref.grid(lmmIH)
            summary(lmmIH.rg1)
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            IH.lsm <- lsmeans(lmmIH, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(IH.lsm, alpha = .50) # To obtain groupings! 
            
            
#FH -model
            lmmFH <- lmer(FH ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = GAF,
                          REML = FALSE) # Simple model FH
            summary(lmmFH)
            Anova(lmmFH)
            
            #* testing the effect adding inundation to the model
            lmmFH.null <- lmer(FH ~ Soil +  (1 | Seed), data = GAF,
                               REML = FALSE) # Simple model for IH
            
            lmmFH.inu <- lmer(FH ~ Soil + Inundation +  (1 | Seed), data = GAF,
                              REML = FALSE) # Adding inundation
            
            anova (lmmFH.null,lmmFH.inu) # Likelihood ratio test to demonstrate no effect of inundation
            
            #* testing the effect adding soil to the model
            lmmFH.0 <- lmer(FH ~ 1 + (1 | Seed), data = GAF,
                            REML = FALSE) # Simple model for IH
            
            lmmFH.soil <- lmer(FH ~ Soil + (1 | Seed), data = GAF,
                               REML = FALSE) # Adding inundation
            
            anova (lmmFH.0,lmmFH.soil)
            summary(lmmFH.soil)
            
            
            # easy backward selection
            
            stimm.step1 <- step(lmmFH) #backward step model selection -random factor dropped
            stimm.step1
            
            #lsmeans FH
            library(lsmeans)
            lmmFH.rg1 <- ref.grid(lmmFH)
            summary(lmmFH.rg1)
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            FH.lsm <- lsmeans(lmmFH, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(FH.lsm, alpha = .50) # To obtain groupings! 
            
            
            
#######    RTIHN --- GLMM ###################

head(GAF)
            
GAF$RHIN.t <- GAF$RHIN + 1
            
            
            #full model including random effect but using Laplace approximation instead 
            glmm_RHIN.t <- glmer(RHIN.t ~ Inundation + Soil + Inundation*Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            summary(glmm_RHIN.t)
            Anova(glmm_RHIN.t)
            
            #*testing relevance of inundation and interaction term
            glmm_RHIN.t0 <- glmer(RHIN.t ~  Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            glmm_RHIN.t1 <- glmer(RHIN.t ~  Soil + Inundation + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            glmm_RHIN.t2 <- glmer(RHIN.t ~  Soil + Inundation + Soil*Inundation + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            anova(glmm_RHIN.t0,glmm_RHIN.t1)
            anova(glmm_RHIN.t0,glmm_RHIN.t1,glmm_RHIN.t2, test="Chisq") 
            
            #*testing importance of random effect
            RHIN.nrandom <- glm(RHIN.t ~ Inundation + Soil + Inundation*Soil, family = gaussian(link = "log"),data = GAF) 
            RHIN.random <- glmer(RHIN.t ~ Inundation + Soil + Inundation*Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            anova(RHIN.random, RHIN.nrandom,test="Chisq") #* X2=0, P=1
            #*Posthoc tests
            lsmeans(glmm_RHIN.t2, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            lsmeans(glmm_RHIN.t2, pairwise ~ Inundation|Soil, adjust ="bonferroni") #within soil treatments among hydrological regimes
            
            
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            RIHN.lsm <- lsmeans(glmm_RHIN.t, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(RIHN.lsm, alpha = .50) # To obtain groupings!

            
#######    RHOAF --- GLMM ###################
            
            GAF$RHOAF.t
            
            GAF$RHOAF.t <- GAF$RHOAF + 1
            
            #full model including random effect but using Laplace approximation instead 
            glmm_RHOAF.t <- glmer(RHOAF.t ~ Inundation + Soil + Inundation*Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            summary(glmm_RHOAF.t)
            Anova(glmm_RHOAF.t)
            
            #*testing relevance of inundation and interaction term
            glmm_RHOAF.t0 <- glmer(RHOAF.t ~  Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            glmm_RHOAF.t1 <- glmer(RHOAF.t ~  Soil + Inundation + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            glmm_RHOAF.t2 <- glmer(RHOAF.t ~  Soil + Inundation + Soil*Inundation + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            anova(glmm_RHOAF.t0,glmm_RHOAF.t1)
            anova(glmm_RHOAF.t0,glmm_RHOAF.t1,glmm_RHOAF.t2, test="Chisq") 
            
            #*testing importance of random effect
            RHOAF.nrandom <- glm(RHOAF.t ~ Inundation + Soil + Inundation*Soil, family = gaussian(link = "log"),data = GAF) 
            RHOAF.random <- glmer(RHOAF.t ~ Inundation + Soil + Inundation*Soil + (1 | Seed), data = GAF, family = gaussian(link = "log")) 
            anova(RHIN.random, RHIN.nrandom,test="Chisq") #* X2=0, P=1
            #*Posthoc tests
            lsmeans(glmm_RHOAF.t2, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            lsmeans(glmm_RHOAF.t2, pairwise ~ Inundation|Soil, adjust ="bonferroni") #within soil treatments among hydrological regimes
           
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            RHOAF.lsm <- lsmeans(glmm_RHOAF.t, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(RHOAF.lsm, alpha = .50) # To obtain groupings! 
            
           
####### Abroveground biomass (AB) 
            
            P$LogADW <- log(P$ADW)
            lmmADW <- lmer(LogADW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = TRUE)
            Anova(lmmADW)
            
            #*LSmeans
            lsmeans(lmmADW, pairwise ~ Inundation*Soil, adjust ="bonferroni")
            lsmeans(lmmADW, pairwise ~ Inundation|Soil, adjust ="bonferroni")
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            ADW.lsm <- lsmeans(lmmADW, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(ADW.lsm, alpha = .50) # To obtain groupings! 

#### Belowground Biomass (BB) 
            
            P$LogBDW <- log(P$BDW)
            lmmBDW <- lmer(LogBDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = TRUE)
            Anova(lmmBDW, type ="III")
            
            #*LSmeans
            lsmeans(lmmBDW, pairwise ~ Inundation*Soil, adjust ="bonferroni")
            lsmeans(lmmBDW, pairwise ~ Inundation|Soil, adjust ="bonferroni")
            
            #*Obtaining Tukey-groupings:
            #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
            BDW.lsm <- lsmeans(lmmBDW, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
            cld(BDW.lsm, alpha = .50) # To obtain groupings! 
            
###### Leaf biomass
            
P$LogLDW <- log(P$LDW)
lmmLDW <- lmer(LogLDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
                           REML = TRUE)
Anova(lmmLDW, type ="III")

#*LSmeans
lsmeans(lmmLDW, pairwise ~ Inundation*Soil, adjust ="bonferroni")
lsmeans(lmmLDW, pairwise ~ Inundation|Soil, adjust ="bonferroni")

#*Obtaining Tukey-groupings:
#-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
LDW.lsm <- lsmeans(lmmLDW, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
cld(LDW.lsm, alpha = .50) # To obtain groupings! 


##### Stem biomass
P$LogSDW <- log(P$SDW)
lmmSDW <- lmer(LogSDW ~ Inundation + Soil + Soil*Inundation +  (1 | Seed), data = P,
               REML = TRUE)
Anova(lmmSDW, type = "III")

          #*LSmeans
          lsmeans(lmmSDW, pairwise ~ Inundation*Soil, adjust ="bonferroni")
          lsmeans(lmmSDW, pairwise ~ Inundation|Soil, adjust ="bonferroni")

          #*Obtaining Tukey-groupings:
          #-------Two LS means that share one or more of the same grouping symbols are not significantly differentat the stated value of alpha
          SDW.lsm <- lsmeans(lmmSDW, pairwise ~ Inundation*Soil, adjust ="bonferroni") #all pairwise comparisons
          cld(SDW.lsm, alpha = .50) # To obtain groupings! 

            
###############################################################################
            ###### SUPPLEMENTAL MATERIALS FIGURES
################################################################################

# To change order of Inundation levels
GAF$Hydrology <-factor(GAF$Inundation,levels=levels(GAF$Inundation)[c(2,1,3)])
P$Hydrology <-factor(P$Inundation,levels=levels(P$Inundation)[c(2,1,3)])

            
            
            ## INITIAL HEIGHT PLOT
            
            # Need to re-calculate lsmeans from new re-organized data set
            lmmIH <- lmer(IH ~ Hydrology + Soil + Soil*Hydrology +  (1 | Seed), data = GAF,REML = FALSE) # Simple model for IH
            lmmIH.rg1 <- ref.grid(lmmIH)
            lsmeansHI <- summary(lmmIH.rg1) #This way ggplot2 can recognize the output as a dataframe
            
            # Making plot in ggplot2
            sumIHplot <- ggplot(lsmeansHI, aes(x= Soil, y= prediction, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=prediction-SE, ymax=prediction+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "Initial Height (cm)") 
            IHplot <- sumIHplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                        axis.text.x = element_text(size=14, colour ="black"),
                                        axis.text.y = element_text(size=14, colour ="black"),
                                        legend.position = "top")
            p7 <-  IHplot + scale_x_discrete(name = "Soil type", labels = c("IN" = "Non sterile", "NO"= "Sterile"))
            #IHplot1<- IHplot + theme(legend.position="none")
            IHplot2 <- p7 + ylim(0,30) 
            pX <- IHplot2 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pIH <- pX + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                              panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            ## FINAL HEIGHT PLOT
            
            # Need to re-calculate lsmeans from new re-organized data set
            lmmFH <- lmer(FH ~ Hydrology + Soil + Soil*Hydrology +  (1 | Seed), data = GAF,REML = FALSE) # Simple model for IH
            lmmFH.rg1 <- ref.grid(lmmFH)
            lsmeansFI <- summary(lmmFH.rg1) #This way ggplot2 can recognize the output as a dataframe
            
            ### FINAL HEIGHT PLOT
            sumFHplot <- ggplot(lsmeansFI, aes(x= Soil, y= prediction, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=prediction-SE, ymax=prediction+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "Final Height (cm)") 
            FHplot <- sumFHplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                        axis.text.x = element_text(size=14, colour ="black"),
                                        axis.text.y = element_text(size=14, colour ="black"),
                                        legend.position = "top")
            p7 <-  FHplot + scale_x_discrete(name = "Soil type", labels = c("IN" = "Non sterile", "NO"= "Sterile"))
            #FHplot1<- FHplot + theme(legend.position="none")
            FHplot2 <- p7 + ylim(0,30) 
            pX <- FHplot2 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pFH <- pX + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                              panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            # RHOAF PLOT
            
            # Need to re-calculate lsmeans from new re-organized data set
            library(lsmeans)
            PQL_RTHO <- glmmPQL(RTHO.t ~ Hydrology + Soil + Hydrology*Soil, ~1 | Seed, family = gaussian(link = "log"),
                                data = GAF, verbose = FALSE) #Note that I am using RTHO.t due to zeros
            
            lmmRTHO.rg1 <- ref.grid(glmm_RTHO2, type="response") # here we type = "response" to obtain back-transformed lsmeans
            summary(lmmRTHO.rg2)
            
            sumRHOAFplot <- ggplot(lsmeansRTHO, aes(x= Soil , y= response, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Inoculation Treatments", y= "RHOAF (mm)") 
            p6 <- sumRHOAFplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                       axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                       axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                       axis.text.x = element_text(size=14, colour ="black"),
                                       axis.text.y = element_text(size=14, colour ="black"),
                                       legend.position = "top")
            p7 <- p6 + scale_x_discrete(name = "Soil type", labels = c("IN" = "Non sterile", "NO"= "Sterile"))
            p12 <- p7 + ylim(0,5.0)
            p13 <- p12 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pRHOAF <- p13 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                  panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            # RHIN PLOT
            
            # Need to re-calculate lsmeans from new re-organized data set
            ????
              
              
              sumRHINplot <- ggplot(sumRHIN, aes(x= Soil , y= RHIN, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=RHIN-se, ymax=RHIN+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("white","light gray","dark gray")) +
              labs(x="Inoculation Treatments", y= "RHIN (mm)") 
            p5 <- sumRHINplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                      axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                      axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                      axis.text.x = element_text(size=14, colour ="black"),
                                      axis.text.y = element_text(size=14, colour ="black"),
                                      legend.position = "top")
            p7 <- p5 + scale_x_discrete(name = "Soil type", labels = c("IN" = "Non sterile", "NO"= "Sterile"))
            p12 <- p7 + ylim(0,3.0)
            p13 <- p12 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pRHIN <- p13 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                 panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            
#################################### FUNGI ###########################################################################

E$Hydrology <-factor(E$Inundation,levels=levels(E$Inundation)[c(2,1,3)])           
#*Re-calculating raw means not categories- basic percentages
            sumTAMF = summarySE(data=E, measurevar = "TAMF", groupvars = c("Soil", "Hydrology"))
            sumTDSE = summarySE(data=E, measurevar = "TDSE", groupvars = c("Soil", "Hydrology"))
            sumTEndo = summarySE(data=E, measurevar = "Tendo", groupvars = c("Soil", "Hydrology"))
            
#### TAMF ######
            
            #* Plot
            TAMFplot <- ggplot(sumTAMF, aes(x= Soil, y= TAMF, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=TAMF-se, ymax=TAMF+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "AMF colonization (%)") 
            p8 <- TAMFplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                   axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                   axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                   axis.text.x = element_text(size=14, colour ="black"),
                                   axis.text.y = element_text(size=14, colour ="black"))
            p9 <- p8 + ylim(0,0.45)
            p10 <- p9 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pTAMF <- p10 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                 panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
#### TDSE
            
            TDSEplot <- ggplot(sumTDSE, aes(x= Soil, y= TDSE, fill= Hydrology)) +
              geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
              geom_errorbar(aes(ymin=TDSE-se, ymax=TDSE+se), width=0.3, colour="black", position=position_dodge(0.9)) +
              theme (axis.title.y = element_text(vjust=1.8),axis.title.x = element_text(vjust=-0.5),axis.title = element_text(face="bold"), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank()) + 
              scale_fill_manual (values = c("blue","light gray","dark blue")) +
              labs(x="Soil Inoculation Treatments", y= "DSE colonization (%)") 
            p11 <- TDSEplot + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                    axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                    axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                    axis.text.x = element_text(size=14, colour ="black"),
                                    axis.text.y = element_text(size=14, colour ="black"))
            p12 <- p11 + ylim(0,0.45)
            p13 <- p12 + theme(panel.background = element_rect(fill = "white", colour = "black"))
            pTDSE <- p13 + theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                 panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
            
            
    # **Plotting all graphs together with a shared legend  
            FIGS3 <- ggarrange(pTAMF, pTDSE, ncol=2, common.legend = TRUE, legend="top")
            #grid.arrange(pIH, pFH, pRTHO, ncol=3)
            ggsave("FIGS3.pdf", plot=FIGS3, device="pdf", scale = 1, width = 20, height = 10, units = "cm", dpi = 300)
            

