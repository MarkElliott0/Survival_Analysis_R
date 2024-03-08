## Survival Analysis in R ## 
## Presentation to UKHSA - 04/03/2024

## This code uses a randomly generated dataset based on chemotherapy treatments for FLOT and ECX

## Install Packages - if required
#install.packages("survival")
#install.packages("survminer")

## Load packages 
library(tidyverse) # Various data manipulation packages 
library(survival) # package for the survival analysis function 
library(survminer) # package for further survival functions 
library(readxl) # if required for data import from excel 


## Load data sets for ECX and FLOT  
ECX <- read_excel("UGI_Study_Dummy_Dataset.xlsx", sheet = 'ECX')
FLOT <- read_excel("UGI_Study_Dummy_Dataset.xlsx", sheet = 'FLOT')


## Basic analysis/ summary on the cohort 
ECX_Summary <- summary(ECX)
FLOT_Summary <- summary(FLOT)

## Quick Analysis on the distribution of those that survived/censored or died in the ECX study 
ECX_Survived <- ECX %>% 
  filter(ECX$Status == 0)

## Quick Analysis on the distribution of those that survived/censored or died in the ECX study 
FLOT_Survived <- FLOT %>% 
  filter(FLOT$Status == 0)

## ECX summary
ECX_Summary_Survived <- summary(ECX_Survived)

# FLOT summary
FLOT_Summary_Survived <- summary(FLOT_Survived)

## T-test for FLOT and ECX - OS Months 
T_Test_Total_Chemo <- t.test(ECX$OS_Months, FLOT$OS_Months, paired = TRUE)
T_Test_Total_Chemo

## T-test for FLOT_Survived vs ECX Survived
T_Test_Survived <- t.test(ECX_Survived$OS_Months, FLOT_Survived$OS_Months, var.equal = TRUE)
T_Test_Survived

## Welch t-statistic
t.test_Welch <- t.test(ECX$OS_Months, FLOT$OS_Months, var.equal=FALSE)

# FLOT variance 
FLOT_Min <- min(FLOT_Survived$OS_Months)

FLOT_Max <- max(FLOT_Survived$OS_Months)

## ECX variance, difference between the Max and Min in OS_Months of survival ECX participants
ECX_Var <- ECX_Max - ECX_Min

## FLOT variance, difference between the Max and Min in OS_Months of survival FLOT participants
FLOT_Var <- FLOT_Max - FLOT_Min

hist_ECX_Survived <- hist(ECX_Survived$OS_Months, main = 'Histogram of ECX Survived by Months', xlab = 'Months')
hist_FLOT_Survived <- hist(FLOT_Survived$OS_Months, main = 'Histogram of FLOT Survived by Months', xlab = 'Months')


## Rejoin datasets using the rbind function, view the summary of the datasets (ECX, FLOT and combined)
## Set time variable,in this particular case, its OS_months, run model and fit the model to the visual function ggsurvplot()
## Plot additional visual without confidence interval
Total_Data <- rbind(ECX, FLOT)
Survived_Data <- rbind(ECX_Survived, FLOT_Survived)
summary(FLOT)
summary(ECX)
summary(Total_Data)

##Boxplot visual for full data - OS_Months by Chemo 
ggplot(data = Data, aes(x = `Neoadj. Chemo.`, y = OS_Months)) + # Define the dataset used and the x and y aesthetics
  geom_boxplot() + # Create a boxplot
  theme_classic() + # Assign a theme to the visuals 
  labs(title = 'OS_Months by Chemo Type - All Data', y = 'Months', x = 'Chemotherapy Type') # label the title and x, y axis

##Boxplot visual for survived data - OS_Months by Chemo 
ggplot(data = Survived_Data, aes(x = `Neoadj. Chemo.`, y = OS_Months)) + # Define the dataset used and the x and y aesthetics 
  geom_boxplot() + # Create a boxplot
  theme_classic() + # Assign a theme to the visuals 
  labs(title = 'OS_Months by Chemo Type - Survived', y = 'Months', x = 'Chemotherapy Type') # label the title and x, y axis

##Survival Difference with deaths included 
surv_diff <- survdiff(Surv(OS_Months, Status) ~ Data$`Neoadj. Chemo.`, data = Data)

surv_diff

## Anova test for FLOT vs ECX OS_Months
fit_Chemo <- aov(OS_Months~`Neoadj. Chemo.`, data = Data)
anova(fit_Chemo)

## Anova test for survived FLOT vs ECX OS_Months
fit_Chemo_surv <- aov(OS_Months~`Neoadj. Chemo.`, data = Survived_Data)
anova(fit_Chemo_surv)

## Create a survival object, usually used as a response variable in a model formula. Argument matching is special for this function, see Details below.
surv_object <- Surv(time = Data$OS_Months, event = Data$Status) # Define the model through time and event arguments

fit_surv_Chemo <- survfit(surv_object ~ Data$`Neoadj. Chemo.`, data = Data) # This function creates survival curves from either a formula (e.g. the Kaplan-Meier), a previously fitted Cox model, or a previously fitted accelerated failure time model.

Chemo_plot <- ggsurvplot(fit_surv_Chemo, data = Data, conf.int = TRUE, pval = TRUE, # Use the above fitted model to create object j, include conf interval and pval as part of the output.
                pval.method = TRUE, # Whether to add a text with the test name used for calculating the pvalue, that corresponds to survival curves' comparison
                xlim=c(0,56), # limit the x-axis
                risk.table = FALSE, # Whether to iniclude the risk table as part of the visual
                tables.theme = theme_cleantable(),legend.labs=c("ECX", "FLOT"), # Creating a theme and defining the legend table
                legend.title="Chemotherapy Type", # legend title
                palette=c("blue2", "pink2"), # colours to use for the survival lines
                title="Kaplan-Meier Curve for Chemotherapy Survival", # title of the visual
                break.time.by = 10,risk.table.fontsize = 4 ,risk.table.height = 0.20, # add breaks and adjust height and font size
                surv.median.line = c("hv"), # add survival median line (hv = horizontal and vertical)
                censor.shape = c("|"), # character or numeric value specifying the point shape of censors. Default value is "+" (3), a sensible choice is "|" (124).
                censor.size = 4, # censor.size: numveric value specifying the point size of censors.
                axes.offset = FALSE, tables.y.text = FALSE) # whether to add offsets to the axis and table y.text if TRUE, the y tick labels of tables will be colored by strata

## Edit the plot object of j and add a scale to the y aesthetic
Chemo_plot$plot <- Chemo_plot$plot + 
  scale_y_continuous(breaks = sort(c(seq(0, 1, 0.1))))

## Plot output of J
Chemo_plot

##Same plot as above but without the confidence interval

Chemo_Plot_NCI <- ggsurvplot(fit_surv_Chemo, data = Data, 
                conf.int = FALSE, 
                pval = TRUE, 
                pval.method = TRUE, 
                xlim=c(0,56), 
                risk.table = FALSE, 
                tables.theme = theme_cleantable(),legend.labs=c("ECX", "FLOT"),
                legend.title="Chemotherapy Type", 
                palette=c("blue2", "pink2"), 
                title="Kaplan-Meier Curve for Chemotherapy Survival", 
                break.time.by = 10,risk.table.fontsize = 4 ,risk.table.height = 0.20, 
                surv.median.line = c("hv"), 
                censor.shape = c("|"), 
                censor.size = 4, 
                axes.offset = FALSE, tables.y.text = FALSE)
Chemo_Plot_NCI$plot <- Chemo_Plot_NCI$plot + 
  scale_y_continuous(breaks = sort(c(seq(0, 1, 0.1))))

#print plot - without confidence intervals
Chemo_Plot_NCI

## Print Survival Analysis Table - summary
print(fit_surv_Chemo)
summary(fit_surv_Chemo)

## Write summary to table output
write.table(sumstats, file = "sumstats.txt", sep = ",", quote = FALSE, row.names = FALSE)

