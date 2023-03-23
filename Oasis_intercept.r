##########################################################################################
## Project: Spatial OASIS 
## File name: Oasis_intercept.r
## Author: Ashika Mani
## Date Created: 3/2/2023
## Last Modified: 3/2/2023
## Description:
##
##  This R program runs the Mimosa model on Spatial Oasis data
##
##########################################################################################

##########################################################################################
## Load the necessary libraries
##########################################################################################

library(fslr)
library(neurobase)
library(dplyr)
library(oasis)
library(oro.nifti)
library(abind)
library(ANTsR)
library(extrantsr)
library(WhiteStripe)
library(pROC)
library(grDevices)
library(EveTemplate)
library(lme4)
.libPaths("/home/ashika/myRpackages")

  
##########################################################################################
## Extract subject ids
##########################################################################################
in.dir <- '/project/Spatial_OASIS/processed/' 
in.dir_T1_GS <- '/project/Spatial_OASIS/Images/'
brainmask <-  '/project/Spatial_OASIS/OASIS_v1/'
out.dir <- '/project/Spatial_OASIS/z_score/'
intercept <- '/project/Spatial_OASIS/Intercept/'
subjectDir <- dir(in.dir)
unique_subjects <- unique(substr(subjectDir, 1, 7))


##########################################################################################
## Randomly split the data into three sets for training, training, and testing
## We have 98 subjects. Put 65 in each of training/training, and 33 in testing.
##########################################################################################
load(file=paste0(brainmask,"filepaths.rda"))
for (i in seq(nrow(filepaths))) {
  # Load files
  visit_id = filepaths$visit_id[i]
  fname = file.path(in.dir, 
                    paste0(visit_id, '_labs_flip.nii.gz'))
    filepaths$EVE[i] = fname
  }

x <- 1:2000
seed = sample(x,1)
set.seed(seed)
sets <- split(unique_subjects, sample(rep(1:2,c(65,33))))

training <- sets$'1'
testing <- sets$'2'

filepaths_train <- filepaths[filepaths$visit_id %in% training,] 
filepaths_test <- filepaths[filepaths$visit_id %in% testing,] 

##########################################################################################
## First, we have to loop through the subjects in the training set and create
## a training dataframe for each subject. 
##########################################################################################

for (i in 1:length(training)){
  subject.id  <- training[i]
  
  print(paste("Beginning: ",subject.id))
  
  # Read in ith subject's images
  T1_training = readnii(filepaths_train$T1[i])
  T2_training = readnii(filepaths_train$T2[i])
  FLAIR_training = readnii(filepaths_train$FLAIR[i])
  PD_training = readnii(filepaths_train$PD[i])
  gold_standard = readnii(filepaths_train$GS[i])
  GS_flip <- neurobase::flip_img(gold_standard , x=TRUE, y = TRUE)
  
  # Create training dataframe for ith subject 
  name <- paste("training_", subject.id, "_dataframe", sep="")
  assign(name,oasis_train_dataframe(flair = FLAIR_training,t1 = T1_training,t2 = T2_training,pd = PD_training,gold_standard = GS_flip,brain_mask = readnii(filepaths$brainmask[i]),preproc=FALSE,normalize=TRUE,cores=4))
  #save(list = name,file=paste0(intercept,"training_", subject.id, "_dataframe", ".rda"))

  print(paste("Finished: ",subject.id))
  
}

##########################################################################################
## Loop through the trained dataframes for all subjects in the training dataset,
## and add a column of the brain region labels to the trained dataframe
##########################################################################################
for (i in 1:length(training)){
  subject.id  <- training[i]
  
  print(paste("Beginning: ",subject.id))
  #load(file=paste0(brainmask,"training_", subject.id, "_dataframe", ".rda"))
  
  # Retrieve ith subject's trained dataframe
  trained_dataframe <- get(ls(pattern=paste("training_", subject.id, "_dataframe", sep="")))
  
  # Read in ith subject's brain labels
  EVE = readNIfTI( paste0(in.dir, subject.id, '_labs_flip.nii.gz') , reorient=FALSE)
 
  # Add column of brain labels to oasis dataframe for ith subject 
  select <- as.vector(EVE)[as.vector(trained_dataframe$voxel_selection)>0]
  trained_dataframe$oasis_dataframe <- cbind(trained_dataframe$oasis_dataframe,region = select)
  
  # Output new dataframe
  out_name <- paste("training_", subject.id, "_dataframe", sep="")
  assign(out_name,trained_dataframe)
  #save(list = out_name,file=paste0(intercept,"training_", subject.id, "_dataframe", ".rda"))
  
  print(paste("Finished: ",subject.id))
  
}

##########################################################################################
## Train model for subjects in training set
##########################################################################################
# for (i in 1:length(training)) {
#   subject.id  <- training[i]
#   load(file= paste0("/project/Spatial_OASIS/Intercept/training_", subject.id, "_dataframe.rda"))
# }

datlist <- mget(ls(pattern="training_"))
# train_oasis_model_whole_brain = oasis_training(datlist[[1]][[1]],datlist[[2]][[1]],datlist[[3]][[1]],datlist[[4]][[1]],datlist[[5]][[1]],datlist[[6]][[1]],datlist[[7]][[1]],datlist[[8]][[1]],datlist[[9]][[1]],datlist[[10]][[1]],datlist[[11]][[1]], datlist[[12]][[1]],datlist[[13]][[1]],datlist[[14]][[1]],datlist[[15]][[1]],datlist[[16]][[1]],datlist[[17]][[1]],datlist[[18]][[1]],datlist[[19]][[1]],datlist[[20]][[1]],datlist[[21]][[1]],datlist[[22]][[1]],datlist[[23]][[1]],datlist[[24]][[1]],datlist[[25]][[1]],datlist[[26]][[1]],datlist[[27]][[1]],datlist[[28]][[1]],datlist[[29]][[1]],datlist[[30]][[1]],datlist[[31]][[1]], datlist[[32]][[1]],datlist[[33]][[1]],datlist[[34]][[1]],datlist[[35]][[1]],datlist[[36]][[1]],datlist[[37]][[1]],datlist[[38]][[1]],datlist[[39]][[1]],datlist[[40]][[1]],datlist[[41]][[1]],datlist[[42]][[1]],datlist[[43]][[1]],datlist[[44]][[1]],datlist[[45]][[1]],datlist[[46]][[1]],datlist[[47]][[1]],datlist[[48]][[1]],datlist[[49]][[1]],datlist[[50]][[1]],datlist[[51]][[1]], datlist[[52]][[1]],datlist[[53]][[1]],datlist[[54]][[1]],datlist[[55]][[1]],datlist[[56]][[1]],datlist[[57]][[1]],datlist[[58]][[1]],datlist[[59]][[1]],datlist[[60]][[1]],datlist[[61]][[1]],datlist[[62]][[1]],datlist[[63]][[1]],datlist[[64]][[1]],datlist[[65]][[1]])
# save(file=paste0(brainmask,"train_oasis_model_whole_brain.rda"),train_oasis_model_whole_brain)


##########################################################################################
## Fit spatially-varying models
##########################################################################################


# Create list of dataframes
list_of_train_dataframes <- list(datlist[[1]][[1]],datlist[[2]][[1]],datlist[[3]][[1]],datlist[[4]][[1]],datlist[[5]][[1]],datlist[[6]][[1]],datlist[[7]][[1]],datlist[[8]][[1]],datlist[[9]][[1]],datlist[[10]][[1]],datlist[[11]][[1]], datlist[[12]][[1]],datlist[[13]][[1]],datlist[[14]][[1]],datlist[[15]][[1]],datlist[[16]][[1]],datlist[[17]][[1]],datlist[[18]][[1]],datlist[[19]][[1]],datlist[[20]][[1]],datlist[[21]][[1]],datlist[[22]][[1]],datlist[[23]][[1]],datlist[[24]][[1]],datlist[[25]][[1]],datlist[[26]][[1]],datlist[[27]][[1]],datlist[[28]][[1]],datlist[[29]][[1]],datlist[[30]][[1]],datlist[[31]][[1]], datlist[[32]][[1]],datlist[[33]][[1]],datlist[[34]][[1]],datlist[[35]][[1]],datlist[[36]][[1]],datlist[[37]][[1]],datlist[[38]][[1]],datlist[[39]][[1]],datlist[[40]][[1]],datlist[[41]][[1]],datlist[[42]][[1]],datlist[[43]][[1]],datlist[[44]][[1]],datlist[[45]][[1]],datlist[[46]][[1]],datlist[[47]][[1]],datlist[[48]][[1]],datlist[[49]][[1]],datlist[[50]][[1]],datlist[[51]][[1]], datlist[[52]][[1]],datlist[[53]][[1]],datlist[[54]][[1]],datlist[[55]][[1]],datlist[[56]][[1]],datlist[[57]][[1]],datlist[[58]][[1]],datlist[[59]][[1]],datlist[[60]][[1]],datlist[[61]][[1]],datlist[[62]][[1]],datlist[[63]][[1]],datlist[[64]][[1]],datlist[[65]][[1]])

# Convert to dataframe
train_vectors_multi <- do.call(rbind, list_of_train_dataframes)
train_vectors_multi <- as.data.frame(train_vectors_multi)
train_vectors_multi <- train_vectors_multi[train_vectors_multi$GoldStandard<=1,]

# Sort by region, and exclude 'background' (region=0)
train_vectors_multi <- train_vectors_multi[order(train_vectors_multi$region),]
train_vectors_multi <- train_vectors_multi[train_vectors_multi$region!=0,]

# Fit the OASIS model separately for each brain region 
# formula1 = GoldStandard ~ FLAIR_10 *FLAIR + FLAIR_20*FLAIR + PD_10 *PD + PD_20 *PD + T2_10 *T2 + T2_20 *T2 + T1_10 *T1 + T1_20 *T1
# oasis_region <- by(train_vectors_multi, train_vectors_multi$region, function(x) glm(formula = formula1, data=x, family=binomial))
# save(file=paste0(brainmask,"train_oasis_model_by_region.rda"),oasis_region)

# # Loop through oasis region and check for convergence errors. #
# region_errors = c()
# for (i in 1:dim(oasis_region)){
#    region_errors[i] <- anyNA(as.vector(oasis_region[[i]]$coefficients))
# }
# names <- getEveMapLabels(type="II")$structure
# names <- names[2:131] # Removing background
# names <- names[sort(unique(train_vectors_multi$region))]
# names[region_errors]

# Fit the OASIS model with random intercepts on brain region
formula2 = GoldStandard ~ FLAIR_10 *FLAIR + FLAIR_20*FLAIR + PD_10 *PD + PD_20 *PD + T2_10 *T2 + T2_20 *T2 + T1_10 *T1 + T1_20 *T1 + (1|region)
oasis_glmer <- glmer(formula = formula2,data = train_vectors_multi,family = binomial,nAGQ=0,control=glmerControl(optimizer = "nloptwrap"))

#save(file=paste0(intercept,"train_oasis_model_randomint.rda"),oasis_glmer)

# Fit the OASIS model with random slopes on brain region
# formula3 = GoldStandard ~ FLAIR_10 *FLAIR + FLAIR_20*FLAIR + PD_10 *PD + PD_20 *PD + T2_10 *T2 + T2_20 *T2 + T1_10 *T1 + T1_20 *T1 + (1 + FLAIR_10 *FLAIR + FLAIR_20*FLAIR + PD_10 *PD + PD_20 *PD + T2_10 *T2 + T2_20 *T2 + T1_10 *T1 + T1_20 *T1 |region)
# oasis_glmer2 <- glmer(formula = formula3,data = train_vectors_multi,family = binomial)
# save(file=paste0(brainmask,"train_oasis_model_random_slope.rda"),oasis_glmer2)


##########################################################################################
## Random intercepts OASIS predictions on testing dataset
##########################################################################################
source('/project/Spatial_OASIS/predict_functions.R') # Load functions for GLM 
#load(file=paste0(brainmask,"train_oasis_model_randomint.rda"))
for (i in 1:length(testing)){
  print(paste("Beginning: ",subject.id))
  subject.id  <- testing[i]
  
  # Read in ith subject's images
  brain_mask = readnii(filepaths_test$brainmask[i])
  FLAIR = readnii(filepaths_test$FLAIR[i])
  T1 = readnii(filepaths_test$T1[i])
  T2 = readnii(filepaths_test$T2[i])
  PD = readnii(filepaths_test$PD[i])
  EVE = readnii(filepaths_test$EVE[i])

  # Create training dataframe for ith subject 
  out <- paste("oasis_glmer_pred_", subject.id, sep="")
  assign(out,oasis_predict_glmer(flair = FLAIR,t1 = T1,t2 = T2,pd = PD,brain_mask = brain_mask, EVE = EVE, preproc=FALSE,normalize=TRUE,model=oasis_glmer,binary=TRUE,cores=4))
  #save(list = out, file=paste0(intercept,"oasis_glmer_pred_", subject.id, ".rda"))
  
  print(paste("Finished: ",subject.id))
  
}
#ls(pattern="oasis_glmer")


##########################################################################################
## Evaluate performance ROC - Overall
##########################################################################################


auc_glmer  <- data.frame(subject=rep(NA,length(testing)))

for (i in 1:length(testing)){
subject.id  <- testing[i]
auc_glmer[i,"subject"] <- subject.id
#load(paste0(intercept,"oasis_glmer_pred_", subject.id, ".rda"))
gold_standard = readnii(filepaths_test$GS[i])
GS_flip <- neurobase::flip_img(gold_standard , x=TRUE, y = TRUE)
brain_mask = readnii(filepaths_test$brainmask[i])
prob_true_GS <- GS_flip[brain_mask==1]
L <- get(paste("oasis_glmer_pred_", subject.id, sep=""))
probability_map <- L[["oasis_map"]]
#writenii(probability_map, file = paste0(brainmask,"pm_" ,subject.id))
#probability_map <- readnii(paste0('/project/Spatial_OASIS/z_score/pm_', subject.id), reorient = TRUE)
prob_true_predict <- probability_map[brain_mask==1]
ROC_obj = roc(as.vector(prob_true_GS), as.vector(prob_true_predict))
auc_glmer[i,"AUC"] <- auc(ROC_obj)[1]
auc_glmer[i,"pAUC"] <- auc(ROC_obj,partial.auc=c(1,0.99))[1]
thresh = ifelse(prob_true_predict>.25,1,0)
  dice_tab = as.data.frame(table(prob_true_GS,thresh))
  TP = dice_tab[4,3]
  FP = dice_tab[3,3]
  FN = dice_tab[2,3]
  DSC <- (2*TP)/((2*TP)+FP + FN)
  auc_glmer[i,"Dice"] <- DSC
  print(subject.id)
  print(auc_glmer[i,"AUC"])
  print(auc_glmer[i,"pAUC"])
  print(DSC)
}

# Calculate AUC and partial AUC
mean(auc_glmer$AUC)
mean(auc_glmer$pAUC)/0.01
mean(auc_glmer$Dice)
df <- data.frame(seed,mean(auc_glmer$AUC),mean(auc_glmer$pAUC)/0.01,mean(auc_glmer$Dice))
write.table(df, file = paste0(intercept,"results_intercept.csv"), sep = ',' , append = TRUE, quote = FALSE,col.names = FALSE, row.names = FALSE)

