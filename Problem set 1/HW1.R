myPackages <- c("lme4", "rsq", "MuMIn", "optimx", "MASS", "dplyr", "ggplot2", 
                "HLMdiag", "lmtest", "lmerTest", "psy", "PowerUpR","nlme", "MVN", 
                "skimr", "tidyverse", "pacman", "lavaan", "semTools", "gsl")

# then check if all the needed packages are installed by checking against the list 
# of already installed packages in R:
installed <- sapply(myPackages, function(p) p %in% rownames(installed.packages()))

# and install any missing packages
if (any(!installed)) {
  install.packages(myPackages[!installed])
}

# then initiate the packages management library:
library(pacman)

# to load the rest of the needed libraries in this lab at once:
p_load(lme4, rsq, MuMIn, optimx, MASS, dplyr, ggplot2, HLMdiag, lmtest, lmerTest, 
       psy, PowerUpR, nlme, MVN, skimr, tidyverse, lavaan, semTools, gsl)

#====================Data loading===========================================#
PISA09.PathAnalysis <- 
  read.csv("https://raw.githubusercontent.com/alibahkali/EDPSY960/main/PISA09.PathAnalysis.csv")

#====================Clean missing records==================================#
# Loop through each column in the data frame to convert missing values into NAs
for(col in names(PISA09.PathAnalysis)) {
  PISA09.PathAnalysis[[col]][PISA09.PathAnalysis[[col]] == 9999] <- NA
}

# Remove records with any missing values from the data frame
PISA09.PathAnalysis <- na.omit(PISA09.PathAnalysis)

#=====================Normality tests=======================================#
# List of continuous variables in the dataframe
variables <- c("CSTRAT", "ELAB", "MEMO", "READING", "ESCS", "GENDER", "IMMIGR")  

results <- sapply(PISA09.PathAnalysis[variables], function(x) {
  shapiro.test(sample(x, min(length(x), 5000), replace = FALSE))
})

print(results)

#===========================================================================#

for (var_name in variables) {
  # Histogram
  hist(PISA09.PathAnalysis[[var_name]], main=paste("Histogram of", var_name), xlab=var_name, breaks=30)
  
  # QQ Plot
  qqnorm(PISA09.PathAnalysis[[var_name]], main=paste("QQ Plot of", var_name))
  qqline(PISA09.PathAnalysis[[var_name]], col="red")
}

#===========================================================================#
# Assuming data_clean is the dataset after listwise deletion
mvn(PISA09.PathAnalysis[-c(6,7)], subset = NULL, mvnTest = c("mardia"), 
    covariance = TRUE, tol = 1e-25, alpha = 0.5,
    scale = FALSE, desc = TRUE, transform = "none", R = 1000,
    univariateTest = c("Lillie"),
    univariatePlot = "qqplot", multivariatePlot = "qq",
    multivariateOutlierMethod = "quan", bc = FALSE, bcType = "rounded",
    showOutliers = FALSE, showNewData = FALSE)

#===========================================================================#
PISA09.PathAnalysis$READINGs <- scale(PISA09.PathAnalysis$READING, 
                                      center = TRUE, scale = TRUE)

#=====================FIRST MODEL==============================================#
PISA.09.FirstModel <- '

# Regressions
#===================#
READING ~ A*MEMO + B*ELAB + C*CSTRAT
MEMO ~ D*ESCS + E*GENDER
ELAB ~ F*ESCS + G*GENDER + H*IMMIGR
CSTRAT ~ I*ESCS + J*GENDER + K*IMMIGR

# Mediation Analysis
#===================#
# Indirect effect of ESCS on READING through MEMO
DA := D*A
# Indirect effect of ESCS on READING through ELAB
FB := F*B
# Indirect effect of ESCS on READING through CSTRAT
IC := I*C
# Total indirect effect of ESCS on READING
DA_FB_IC := DA + FB + IC
#===================#
# Indirect effect of GENDER on READING through MEMO
EA := E*A
# Indirect effect of GENDER on READING through ELAB
GB := G*B
# Indirect effect of GENDER on READING through CSTRAT
JC := J*C
# Total indirect effect of GENDER on READING
EA_GB_JC := EA + GB + JC
#===================#
# Indirect effect of IMMIGR on READING through ELAB
HB := H*B
# Indirect effect of IMMIGR on READING through CSTRAT
KC := K*C
# Total indirect effect of IMMIGR on READING
HB_KC := HB + KC
'
#==============================================================================#
# Generate summary of path analysis
PISA.09.FirstModel.fit <- sem(PISA.09.FirstModel, data = PISA09.PathAnalysis, estimator = "ml")

summary(PISA.09.FirstModel.fit, fit.measures = TRUE, rsq = T)
#==============================================================================#
modindices(PISA.09.FirstModel.fit, standardized=TRUE, power=TRUE, delta=0.1, 
           alpha=.05, high.power=.80) %>% 
  filter(op != "~~") %>% 
  arrange(desc(mi))


modindices(PISA.09.FirstModel.fit, standardized=TRUE, power=TRUE, delta=0.1, 
           alpha=.05, high.power=.80) %>% 
  arrange(desc(mi))

#=====================SECOND MODEL=============================================#
PISA.09.SecondModel <- '

# Regressions
#===================#
READING ~ A*MEMO + B*ELAB + C*CSTRAT
MEMO ~ D*ESCS + E*GENDER + L*CSTRAT
ELAB ~ F*ESCS + G*GENDER + H*IMMIGR
CSTRAT ~ I*ESCS + J*GENDER + K*IMMIGR

# Mediation Analysis
#===================#
# Indirect effect of ESCS on READING through MEMO
DA := D*A
# Indirect effect of ESCS on READING through ELAB
FB := F*B
# Indirect effect of ESCS on READING through CSTRAT
IC := I*C
# Indirect effect of ESCS on READING through CSTRAT & MEMOR (new to second model)
ILA:= I*L*A
# Total indirect effect of ESCS on READING
DA_FB_IC_ILA := DA + FB + IC + ILA
#===================#
# Indirect effect of GENDER on READING through MEMO
EA := E*A
# Indirect effect of GENDER on READING through ELAB
GB := G*B
# Indirect effect of GENDER on READING through CSTRAT
JC := J*C
# Indirect effect of GENDER on READING through CSTRAT & MEMO (new to second model)
JLA := J*L*A
# Total indirect effect of GENDER on READING
EA_GB_JC_JLA := EA + GB + JC + JLA
#===================#
# Indirect effect of IMMIGR on READING through ELAB
HB := H*B
# Indirect effect of IMMIGR on READING through CSTRAT
KC := K*C
# Indirect effect of IMMIGR on READING through CSTRAT & MEMO (new to second model)
KLA := K*L*A
# Total indirect effect of IMMIGR on READING
HB_KC_KLA:= HB + KC + KLA
'
#==============================================================================#
# Generate summary of path analysis
PISA.09.SecondModel.fit <- sem(PISA.09.SecondModel, data = PISA09.PathAnalysis, estimator = "ml")

summary(PISA.09.SecondModel.fit, fit.measures = TRUE, rsq = T)

#==============================================================================#
modindices(PISA.09.SecondModel.fit, standardized=TRUE, power=TRUE, delta=0.1, 
           alpha=.05, high.power=.80) %>% 
  filter(op != "~~") %>% 
  arrange(desc(mi))


modindices(PISA.09.SecondModel.fit, standardized=TRUE, power=TRUE, delta=0.1, 
           alpha=.05, high.power=.80) %>% 
  arrange(desc(mi))

#=============================THIRD MODEL======================================#
PISA.09.ThirdModel <- '

# Regressions
#===================#
READING ~ A*MEMO + B*ELAB + C*CSTRAT
MEMO ~ D*ESCS + E*GENDER + L*CSTRAT
ELAB ~ F*ESCS + G*GENDER + H*IMMIGR + M*CSTRAT
CSTRAT ~ I*ESCS + J*GENDER + K*IMMIGR

# Mediation Analysis
#===================#
# Indirect effect of ESCS on READING through MEMO
DA := D*A
# Indirect effect of ESCS on READING through ELAB
FB := F*B
# Indirect effect of ESCS on READING through CSTRAT
IC := I*C
# Indirect effect of ESCS on READING through CSTRAT & MEMOR 
ILA := I*L*A
# Indirect effect of ESCS on READING through CSTRAT & ELAB (new to third model)
IMB := I*M*B
# Total indirect effect of ESCS on READING
DA_FB_IC_ILA_IMB := DA + FB + IC + ILA + IMB
#===================#
# Indirect effect of GENDER on READING through MEMO
EA := E*A
# Indirect effect of GENDER on READING through ELAB
GB := G*B
# Indirect effect of GENDER on READING through CSTRAT
JC := J*C
# Indirect effect of GENDER on READING through CSTRAT & MEMO 
JLA := J*L*A
# Indirect effect of GENDER on READING through CSTRAT & ELAB (new to third model)
JMB := J*M*B
# Total indirect effect of GENDER on READING
EA_GB_JC_JLA_JMB := EA + GB + JC + JLA + JMB
#===================#
# Indirect effect of IMMIGR on READING through ELAB
HB := H*B
# Indirect effect of IMMIGR on READING through CSTRAT
KC := K*C
# Indirect effect of IMMIGR on READING through CSTRAT & MEMO 
KLA := K*L*A
# Indirect effect of IMMIGR on READING through CSTRAT & ELAB (new to third model)
KMB := K*M*B
# Total indirect effect of IMMIGR on READING
HB_KC_KLA_KMB := HB + KC + KLA + KMB
'
#==============================================================================#
# Generate summary of path analysis
PISA.09.ThirdModel.fit <- sem(PISA.09.ThirdModel, data = PISA09.PathAnalysis, estimator = "ml")

summary(PISA.09.ThirdModel.fit, fit.measures = TRUE, rsq = T)
#==============================================================================#
modindices(PISA.09.ThirdModel.fit, standardized=TRUE, power=TRUE, delta=0.1, 
           alpha=.05, high.power=.80) %>% 
  filter(op != "~~") %>% 
  arrange(desc(mi))


modindices(PISA.09.ThirdModel.fit, standardized=TRUE, power=TRUE, delta=0.1, 
           alpha=.05, high.power=.80) %>% 
  arrange(desc(mi))

#=============================FOURTH MODEL======================================#
PISA.09.FourthModel <- '

# Regressions
#===================#
READING ~ A*MEMO + B*ELAB + C*CSTRAT
MEMO ~ D*ESCS + E*GENDER + L*CSTRAT + N*ELAB
ELAB ~ F*ESCS + G*GENDER + H*IMMIGR + M*CSTRAT
CSTRAT ~ I*ESCS + J*GENDER + K*IMMIGR

# Mediation Analysis
#===================#
# Indirect effect of ESCS on READING through MEMO
DA := D*A
# Indirect effect of ESCS on READING through ELAB
FB := F*B
# Indirect effect of ESCS on READING through ELAB & MEMOR
FNA := F*N*A
# Indirect effect of ESCS on READING through CSTRAT
IC := I*C
# Indirect effect of ESCS on READING through CSTRAT & MEMOR 
ILA := I*L*A
# Indirect effect of ESCS on READING through CSTRAT & ELAB 
IMB := I*M*B
# Indirect effect of ESCS on READING through CSTRAT & ELAB & MEMOR (new to the fourth model)
IMNA := I*M*N*A
# Total indirect effect of ESCS on READING
DA_FB_IC_ILA_IMB := DA + FB + IC + ILA + IMB + FNA + IMNA
#===================#
# Indirect effect of GENDER on READING through MEMO
EA := E*A
# Indirect effect of GENDER on READING through ELAB
GB := G*B
# Indirect effect of GENDER on READING through CSTRAT
JC := J*C
# Indirect effect of GENDER on READING through CSTRAT & MEMO 
JLA := J*L*A
# Indirect effect of GENDER on READING through CSTRAT & ELAB 
JMB := J*M*B
# Indirect effect of GENDER on READING through CSTRAT & ELAB & MEMO (new to the fourth model)
JMNA := J*M*N*A
# Total indirect effect of GENDER on READING
EA_GB_JC_JLA_JMB := EA + GB + JC + JLA + JMB + JMNA
#===================#
# Indirect effect of IMMIGR on READING through ELAB
HB := H*B
# Indirect effect of IMMIGR on READING through CSTRAT
KC := K*C
# Indirect effect of IMMIGR on READING through CSTRAT & MEMO 
KLA := K*L*A
# Indirect effect of IMMIGR on READING through CSTRAT & ELAB
KMB := K*M*B
# Indirect effect of IMMIGR on READING through CSTRAT & ELAB & MEMO (new to the fourth model)
KMNA := K*M*N*A
# Total indirect effect of IMMIGR on READING
HB_KC_KLA_KMB := HB + KC + KLA + KMB + KMNA
'
#==============================================================================#
# Generate summary of path analysis
PISA.09.FourthModel.fit <- sem(PISA.09.FourthModel, data = PISA09.PathAnalysis, estimator = "ml")

summary(PISA.09.FourthModel.fit, fit.measures = TRUE, rsq = T)
#==============================================================================#
modindices(PISA.09.FourthModel.fit, standardized=TRUE, power=TRUE, delta=0.1, 
           alpha=.05, high.power=.80) %>% 
  filter(op != "~~") %>% 
  arrange(desc(mi))


modindices(PISA.09.FourthModel.fit, standardized=TRUE, power=TRUE, delta=0.1, 
           alpha=.05, high.power=.80) %>% 
  arrange(desc(mi))

#=============================FIFTH MODEL======================================#
PISA.09.FifthModel <- '

# Regressions
#===================#
READING ~ A*MEMO + B*ELAB + C*CSTRAT + O*ESCS + Q*IMMIGR + P*GENDER 
MEMO ~ D*ESCS + E*GENDER + L*CSTRAT + N*ELAB
ELAB ~ F*ESCS + G*GENDER + H*IMMIGR + M*CSTRAT
CSTRAT ~ I*ESCS + J*GENDER + K*IMMIGR

# Mediation Analysis
#===================#
# Indirect effect of ESCS on READING through MEMO
DA := D*A
# Indirect effect of ESCS on READING through ELAB
FB := F*B
# Indirect effect of ESCS on READING through ELAB & MEMOR
FNA := F*N*A
# Indirect effect of ESCS on READING through CSTRAT
IC := I*C
# Indirect effect of ESCS on READING through CSTRAT & MEMOR 
ILA := I*L*A
# Indirect effect of ESCS on READING through CSTRAT & ELAB 
IMB := I*M*B
# Indirect effect of ESCS on READING through CSTRAT & ELAB & MEMOR 
IMNA := I*M*N*A
# Total indirect effect of ESCS on READING
DA_FB_IC_ILA_IMB := DA + FB + IC + ILA + IMB + FNA + IMNA
#===================#
# Indirect effect of GENDER on READING through MEMO
EA := E*A
# Indirect effect of GENDER on READING through ELAB
GB := G*B
# Indirect effect of GENDER on READING through CSTRAT
JC := J*C
# Indirect effect of GENDER on READING through CSTRAT & MEMO 
JLA := J*L*A
# Indirect effect of GENDER on READING through CSTRAT & ELAB 
JMB := J*M*B
# Indirect effect of GENDER on READING through CSTRAT & ELAB & MEMO 
JMNA := J*M*N*A
# Total indirect effect of GENDER on READING
EA_GB_JC_JLA_JMB := EA + GB + JC + JLA + JMB + JMNA
#===================#
# Indirect effect of IMMIGR on READING through ELAB
HB := H*B
# Indirect effect of IMMIGR on READING through CSTRAT
KC := K*C
# Indirect effect of IMMIGR on READING through CSTRAT & MEMO 
KLA := K*L*A
# Indirect effect of IMMIGR on READING through CSTRAT & ELAB
KMB := K*M*B
# Indirect effect of IMMIGR on READING through CSTRAT & ELAB & MEMO 
KMNA := K*M*N*A
# Total indirect effect of IMMIGR on READING
HB_KC_KLA_KMB := HB + KC + KLA + KMB + KMNA
'
#==============================================================================#
# Generate summary of path analysis
PISA.09.FifthModel.fit <- sem(PISA.09.FifthModel, data = PISA09.PathAnalysis, estimator = "ml")

summary(PISA.09.FifthModel.fit, fit.measures = TRUE, rsq = T)
#==============================================================================#
modindices(PISA.09.FifthModel.fit, standardized=TRUE, power=TRUE, delta=0.1, 
           alpha=.05, high.power=.80) %>% 
  filter(op != "~~") %>% 
  arrange(desc(mi))


modindices(PISA.09.FifthModel.fit, standardized=TRUE, power=TRUE, delta=0.1, 
           alpha=.05, high.power=.80) %>% 
  arrange(desc(mi))

varTable(PISA.09.FifthModel.fit)

#==============================================================================#
new_prompt <- taskCallbackManager()

new_prompt$add(function(expr, value, ok, visible) {
  options("prompt" = format(Sys.time(), "[%H:%M:%S]> "));
  return(TRUE) },
  name = "promptHandler")
