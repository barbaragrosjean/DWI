
# WARNING !! Need installation on the server cannot run the code as it is done for now 
#install.packages('readr', repos='http://cran.us.r-project.org')
#install.packages('dplyr', repos='http://cran.us.r-project.org')
#install.packages("ggplot2", repos = "http://cran.us.r-project.org")

# I run it locally -> don't forget to change the path to 
# the data which are in the derivative folder!

library(stats)
library(readr)
library(dplyr)
library(ggplot2)


#------------------------------------------------------------#
#### BEHAVIORAL DATA #### 
#------------------------------------------------------------#

# Load behavioral data
# file_name <- file.path(data_path, 'behav_output', 'gain_blinded.csv') for server
file_name <- file.path('behav.csv') 
df_gain_diff <- read.csv(file_name)
df_behav <- df_gain_diff$gain

#------------------------------------------------------------#
#### SEED BASED  #### 
#------------------------------------------------------------#
# 1. LOAD DATA
file_name <- file.path('seed_metric_df.csv')
seed_metric_df <- read_csv(file_name)

# 2. FIT THE GLM

# binomial (for logistic regression), gaussian (for linear regression), poisson (for count data)...
glm_model_seed <- glm(df_behav ~ v_d_Ca_L + v_d_Ca_R + vm_dl_PU_L + vm_dl_PU_R, data = cbind(df_behav,seed_metric_df), family = gaussian()) # check family 

# View the ANOVA of the model
anova_table <- anova(glm_model_seed, test="F")

# Save the ANOVA table to a .csv file
write.csv(anova_table, file = "seed_metric_anova.csv")

# Visualization
ggplot(seed_metric_df, aes(x=v_d_Ca_L + v_d_Ca_R + vm_dl_PU_L + vm_dl_PU_R, y=df_behav)) + geom_point() + geom_smooth(method="glm", method.args=list(family="gaussian"))

#------------------------------------------------------------#
#### ROI2ROI - PUTAMEN RIGHT ONLY #### 
#------------------------------------------------------------#
# Consider the tracts putamen - .. only
# 1. LOAD DATA
file_name <- file.path('Pu.csv')
connectivity_df <- read_csv(file_name)

# 2. FIT THE GLM
# binomial (for logistic regression), gaussian (for linear regression), poisson (for count data)...
#glm_model_tract_Pu <- glm(df_behav~., data = connectivity_df, family = gaussian())
glm_model_tract_Pu <- glm(df_behav~., data = cbind(df_behav, connectivity_df), family = gaussian())


# View the ANOVA of the model
anova_table <- anova(glm_model_seed, test="F")

# Save the ANOVA table to a .csv file
write.csv(anova_table, file = "Pu_anova.csv")

# Visualization
#ggplot(connectivity_df, aes(x=., y=df_behav)) + geom_point() + geom_smooth(method="glm", method.args=list(family="gaussian"))


#------------------------------------------------------------#
#### ROI2ROI - CAUDATE RIGHT ONLY #### 
#------------------------------------------------------------#
# Consider the tracs caudate - .. only
# 1. LOAD DATA
file_name <- file.path('Ca.csv')
connectivity_df <- read_csv(file_name)

# 2. FIT THE GLM
# binomial (for logistic regression), gaussian (for linear regression), poisson (for count data)...
glm_model_tract_Ca <- glm(df_behav ~., data = cbind(df_behav,connectivity_df), family = gaussian())

# View the ANOVA of the model
anova_table <- anova(glm_model_seed, test="F")

# Save the ANOVA table to a .csv file
write.csv(anova_table, file = "Ca_anova.csv")

# Visualization
#ggplot(connectivity_df, aes(x=., y=df_behav)) + geom_point() + geom_smooth(method="glm", method.args=list(family="gaussian"))


#------------------------------------------------------------#
#### ROI2ROI -  #### 
#------------------------------------------------------------#
# Consider the all network connected to putamen (with tracts in betweeen the network)


# 1. LOAD DATA
file_name <- file.path('Pu_network.csv')
connectivity_df <- read_csv(file_name)

# 2. FIT THE GLM
# binomial (for logistic regression), gaussian (for linear regression), poisson (for count data)...
glm_model_tract_Pu_network <- glm(df_behav ~ ., data = cbind(df_behav,connectivity_df), family = gaussian())

# View the ANOVA of the model
anova_table <- anova(glm_model_seed, test="F")

# Save the ANOVA table to a .csv file
write.csv(anova_table, file = "Pu_network_anova.csv")

# Visualization
#ggplot(connectivity_df, aes(x=., y=df_behav)) + geom_point() + geom_smooth(method="glm", method.args=list(family="gaussian"))
