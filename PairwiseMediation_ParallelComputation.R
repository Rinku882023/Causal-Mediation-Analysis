
# Load necessary libraries
library(stats)  
library(Hmisc)
library(mediation)
library(doParallel)

d1<-data6
d1<-na.omit(d1)
#ensure your data is complete with no missing observations
#order your dataset's columns by: exposures, mediators, followed by covariates and outcome
#transform your exposures and mediators as you see fit/necessary according normality assumptions

# Initialize cluster
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# Define variables
nexp <- 1  # Define number of exposure analytes in your dataset
nmed <- 26743  # Define number of mediators in your dataset
ncovars <- 14  # datasheet number of covariates and outcome variable in your dataset


# Define function for mediation analysis
perform_mediation <- function(i, j, d1) {
  d2 <- d1[, c(i, j, ((nexp + nmed + 1):(nexp + nmed + ncovars)))]
  set.seed(111)
  # Define your covariates as needed
  m <- lm(d2[, 2] ~ d2[, 1] + Age_at_collect + Gender.x + Race_Group.x + Smoking_flag_Biobank + CD4Tnv + CD8Tnv + Mono + Neu + NK + Education + BMI_most_recent + Gap_BMI.Rec_PlasmaCollection+Gap_VitD_25OH.Rec_PlasmaCollection, data = d2)
  # Define your outcome variable below
  y <- lm(d2$fev1pctpred ~ d2[, 1] + d2[, 2] + Age_at_collect + Gender.x + Race_Group.x + Smoking_flag_Biobank + CD4Tnv + CD8Tnv + Mono + Neu + NK + Education + BMI_most_recent + Gap_BMI.Rec_PlasmaCollection+Gap_VitD_25OH.Rec_PlasmaCollection, data = d2)
  med <- mediation::mediate(m, y, sims = 1000, treat = "d2[, 1]", mediator = "d2[, 2]")  # Ensure mediation:: prefix
  return(cbind(nobs(y), med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p, med$d1, med$d1.ci[1], 
               med$d1.ci[2], med$d1.p, med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p, med$z1, 
               med$z1.ci[1], med$z1.ci[2], med$z1.p, med$n0, med$n0.ci[1], med$n0.ci[2], 
               med$n0.p, med$n1, med$n1.ci[1], med$n1.ci[2], med$n1.p, med$tau.coef, med$tau.ci[1], 
               med$tau.ci[2], med$tau.p, med$d.avg, med$d.avg.ci[1], med$d.avg.ci[2], med$d.avg.p,
               med$z.avg, med$z.avg.ci[1], med$z.avg.ci[2], med$z.avg.p, med$n.avg, med$n.avg.ci[1], 
               med$n.avg.ci[2], med$n.avg.p))
}

# Loop to conduct pairwise mediation with multiple exposures and mediators
results <- foreach(i = 1:nexp, .combine = rbind) %:%
  foreach(j = (nexp + 1):(nexp + nmed), .combine = rbind) %dopar% {
    perform_mediation(i, j, d1)
  }

# Combine results into med.results_fev1 dataframe
med.results_vitD_fev1 <- as.data.frame(results)

colnames(med.results_vitD_fev1) <- c('nobs', 'ACME.C', 'ACME.C.lo', 'ACME.C.hi', 'ACME.C.Pval', 'ACME.T', 'ACME.T.lo',
                                     'ACME.T.hi', 'ACME.T.pval', 'ADE.C', 'ADE.C.lo', 'ADE.C.hi', 'ADE.C.Pval', 'ADE.T',
                                     'ADE.T.lo', 'ADE.T.hi', 'ADE.T.pval', 'PMed.C', 'PMed.C.lo', 'PMed.C.hi', 'PMed.C.pval',
                                     'PMed.T', 'PMed.T.lo', 'PMed.T.hi', 'PMed.T.pval', 'TE', 'TE.lo', 'TE.hi', 'TE.pval',
                                     'ACME.avg', 'ACME.avg.lo', 'ACME.avg.hi', 'ACME.avg.pval', 'ADE.avg', 'ADE.avg.lo',
                                     'ADE.avg.hi', 'ADE.avg.pval', 'PMed.avg', 'PMed.avg.lo', 'PMed.avg.hi', 'PMed.avg.pval')


# Assign rownames
med.results_vitD_fev1$Variable <- paste(colnames(d1)[1], colnames(d1[2:26744), sep = '.')

# Stop the cluster
stopCluster(cl)
