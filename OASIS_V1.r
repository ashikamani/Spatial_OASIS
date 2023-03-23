##########################################################################################
## Project: Spatial OASIS 
## File name: Oasis_V1.r
## Author: Ashika Mani
## Date Created: 10/28/2022
## Last Modified: 3/20/2023
## Description:
##
##  This R program runs the Oasis model on Spatial Oasis data
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
.libPaths("/home/ashika/myRpackages")

  
##########################################################################################
## Extract subject ids
##########################################################################################
in.dir <- '/project/Spatial_OASIS/processed/' 
in.dir_T1_GS <- '/project/Spatial_OASIS/Images/'
brainmask <-  '/project/Spatial_OASIS/OASIS_v1/'
out.dir <- '/project/Spatial_OASIS/z_score/'
subjectDir <- dir(in.dir)
unique_subjects <- unique(substr(subjectDir, 1, 7))


##########################################################################################
## Creating a filename matrix
##########################################################################################


T1_files = list.files(path =
                        file.path(in.dir_T1_GS),
                      pattern = "T1[.]nii[.]gz",
                      full.names = TRUE)
T2_files = list.files(path =
                        file.path(in.dir),
                      pattern = "T2_N4_reg_stripped[.]nii[.]gz",
                      full.names = TRUE)
FLAIR_files = list.files(path =
                           file.path(in.dir),
                         pattern = "FLAIR_N4_reg_stripped[.]nii[.]gz",
                         full.names = TRUE)
PD_files = list.files(path =
                        file.path(in.dir),
                      pattern = "PD_N4_reg_stripped[.]nii[.]gz",
                      full.names = TRUE)
GS_files = list.files(path =
                        file.path(in.dir_T1_GS),
                      pattern = "_lesion_mask[.]nii[.]gz",
                      full.names = TRUE)
filepaths = data.frame(T1 = T1_files, T2 = T2_files,
                       FLAIR = FLAIR_files, PD = PD_files,GS = GS_files,
                       stringsAsFactors = FALSE)

# Obtain visit_id
have_data = nrow(filepaths) > 0
if (have_data) {
  ss = gsub(".*[/]([^.]+)[_].*", "\\1", nii.stub(filepaths$T1))
  filepaths$visit_id = sapply(ss, function(x) x)
  filepaths
}

create_brain_mask = function(...) {
  x = list(...)
  x = check_nifti(x)
  x = lapply(x, function(img) {
    img > 0
  })
  mask = Reduce("|", x)
  mask = datatyper(mask)
  mask
}

#Creating Brain mask for all subjects
filepaths$brainmask = NA
for (i in seq(nrow(filepaths))) {
  # Load files
  visit_id = filepaths$visit_id[i]
  fname = file.path(brainmask, 
                    paste0("brainmask_",
                           visit_id, ".nii.gz"))
    # T1 = readnii(filepaths$T1[i])
    # T2 = readnii(filepaths$T2[i])
    # FLAIR = readnii(filepaths$FLAIR[i])
    # PD = readnii(filepaths$PD[i])
   # brain_mask = T1 > 0 
  #Save brain mask to local working directory
   #writenii(brain_mask, 
        #    filename = fname)
    filepaths$brainmask[i] = fname
  }

  #save(file=paste0(brainmask,"filepaths.rda"),filepaths)

  

# Checking orientation 
#  T1 = readnii(filepaths$T1[3])
#     T2 = readnii(filepaths$T2[3])
#     FLAIR = readnii(filepaths$FLAIR[3])
#     PD = readnii(filepaths$PD[3])

#     brain = readnii(filepaths$brainmask[3])

# jpeg(paste0(brainmask,"flairandt1.jpg"))
# double_ortho(T2,FLAIR)
# dev.off()


##########################################################################################
## Randomly split the data into three sets for training, training, and testing
## We have 98 subjects. Put 65 in each of training/training, and 33 in testing.
##########################################################################################
load(file=paste0(brainmask,"filepaths.rda"))
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
  #save(list = name,file=paste0(brainmask,"training_", subject.id, "_dataframe", ".rda"))

  print(paste("Finished: ",subject.id))

  
}

##########################################################################################
## Train model for subjects in training set
##########################################################################################
# for (i in 1:length(training)) {
#   subject.id  <- training[i]
#   load(file= paste0("/project/Spatial_OASIS/z_score/training_", subject.id, "_dataframe.rda"))
# }

datlist <- mget(ls(pattern="training_"))
train_oasis_model_whole_brain = oasis_training(datlist[[1]][[1]],datlist[[2]][[1]],datlist[[3]][[1]],datlist[[4]][[1]],datlist[[5]][[1]],datlist[[6]][[1]],datlist[[7]][[1]],datlist[[8]][[1]],datlist[[9]][[1]],datlist[[10]][[1]],datlist[[11]][[1]], datlist[[12]][[1]],datlist[[13]][[1]],datlist[[14]][[1]],datlist[[15]][[1]],datlist[[16]][[1]],datlist[[17]][[1]],datlist[[18]][[1]],datlist[[19]][[1]],datlist[[20]][[1]],datlist[[21]][[1]],datlist[[22]][[1]],datlist[[23]][[1]],datlist[[24]][[1]],datlist[[25]][[1]],datlist[[26]][[1]],datlist[[27]][[1]],datlist[[28]][[1]],datlist[[29]][[1]],datlist[[30]][[1]],datlist[[31]][[1]], datlist[[32]][[1]],datlist[[33]][[1]],datlist[[34]][[1]],datlist[[35]][[1]],datlist[[36]][[1]],datlist[[37]][[1]],datlist[[38]][[1]],datlist[[39]][[1]],datlist[[40]][[1]],datlist[[41]][[1]],datlist[[42]][[1]],datlist[[43]][[1]],datlist[[44]][[1]],datlist[[45]][[1]],datlist[[46]][[1]],datlist[[47]][[1]],datlist[[48]][[1]],datlist[[49]][[1]],datlist[[50]][[1]],datlist[[51]][[1]], datlist[[52]][[1]],datlist[[53]][[1]],datlist[[54]][[1]],datlist[[55]][[1]],datlist[[56]][[1]],datlist[[57]][[1]],datlist[[58]][[1]],datlist[[59]][[1]],datlist[[60]][[1]],datlist[[61]][[1]],datlist[[62]][[1]],datlist[[63]][[1]],datlist[[64]][[1]],datlist[[65]][[1]])
#save(file=paste0(brainmask,"train_oasis_model_whole_brain.rda"),train_oasis_model_whole_brain)

##########################################################################################
## Whole-brain OASIS predictions on testing dataset
##########################################################################################
for (i in 1:length(testing)){
  subject.id  <- testing[i]
  
  print(paste("Beginning: ",subject.id))
  
  # Read in ith subject's images
  brain_mask = readnii(filepaths_test$brainmask[i])
  FLAIR = readnii(filepaths_test$FLAIR[i])
  T1 = readnii(filepaths_test$T1[i])
  T2 = readnii(filepaths_test$T2[i])
  PD = readnii(filepaths_test$PD[i])
  
  # Create training dataframe for ith subject 
  out <- paste("oasis_testingpred_subject_", subject.id, sep="")
  assign(out,oasis_predict(flair = FLAIR,t1 = T1,t2 = T2,pd = PD,brain_mask = brain_mask,preproc=FALSE,normalize=TRUE,model=train_oasis_model_whole_brain,binary=TRUE,cores=4))
  #save(list = out, file=paste0(brainmask,"oasis_testing_subject_", subject.id,".rda"))
  print(paste("Finished: ",subject.id))

}

##########################################################################################
## Evaluate performance ROC - Overall
##########################################################################################

# subject.id  <- testing[1]
# #load(file= paste0(out.dir,"oasis_testingpred_subject_",subject.id,".rda"))
# gold_standard = readnii(filepaths_test$GS[1])
# GS_flip <- neurobase::flip_img(gold_standard , x=TRUE, y = TRUE)
# brain_mask = readnii(filepaths_test$brainmask[1])
# prob_true_GS <- GS_flip[brain_mask==1]
# L <- get(paste("oasis_testingpred_subject_", subject.id, sep=""))
# probability_map <- L[["oasis_map"]]
# #writenii(probability_map, file = paste0(brainmask,"pm_" ,subject.id))
# #probability_map <- readnii(paste0('/project/Spatial_OASIS/z_score/pm_', subject.id), reorient = TRUE)
# prob_true_predict <- probability_map[brain_mask==1]
# roc(as.vector(prob_true_GS), as.vector(prob_true_predict))

auc_full  <- data.frame(subject=rep(NA,length(testing)))

for (i in 1:length(testing)){
subject.id  <- testing[i]
auc_full[i,"subject"] <- subject.id
#load(file= paste0(out.dir,"oasis_testingpred_subject_",subject.id,".rda"))
gold_standard = readnii(filepaths_test$GS[i])
GS_flip <- neurobase::flip_img(gold_standard , x=TRUE, y = TRUE)
brain_mask = readnii(filepaths_test$brainmask[i])
prob_true_GS <- GS_flip[brain_mask==1]
L <- get(paste("oasis_testingpred_subject_", subject.id, sep=""))
probability_map <- L[["oasis_map"]]
#writenii(probability_map, file = paste0(brainmask,"pm_" ,subject.id))
#probability_map <- readnii(paste0('/project/Spatial_OASIS/z_score/pm_', subject.id), reorient = TRUE)
prob_true_predict <- probability_map[brain_mask==1]
ROC_obj = roc(as.vector(prob_true_GS), as.vector(prob_true_predict))
auc_full[i,"AUC"] <- auc(ROC_obj)[1]
auc_full[i,"pAUC"] <- auc(ROC_obj,partial.auc=c(1,0.99))[1]
thresh = ifelse(prob_true_predict>.25,1,0)
  dice_tab = as.data.frame(table(prob_true_GS,thresh))
  TP = dice_tab[4,3]
  FP = dice_tab[3,3]
  FN = dice_tab[2,3]
  DSC <- (2*TP)/((2*TP)+FP + FN)
  auc_full[i,"Dice"] <- DSC
  print(subject.id)
  print(auc_full[i,"AUC"])
  print(auc_full[i,"pAUC"])
}

# Calculate AUC and partial AUC
mean(auc_full$AUC)
mean(auc_full$pAUC)/0.01
mean(auc_full$Dice)
df <- data.frame(seed,mean(auc_full$AUC),mean(auc_full$pAUC)/0.01,mean(auc_full$Dice))
write.table(df, file = paste0(brainmask,"results_whole.csv"), sep = ',' , append = TRUE, quote = FALSE,col.names = FALSE, row.names = FALSE)
