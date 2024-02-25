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
