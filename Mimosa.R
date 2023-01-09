##########################################################################################
## Project: Spatial OASIS 
## File name: Mimosa.R
## Author: Ashika Mani
## Date Created: 10/28/2022
## Last Modified: 12/24/2022
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
library(mimosa)
library(dplyr)
library(oasis)
library(oro.nifti)
library(abind)
library(ANTsR)
library(extrantsr)
library(WhiteStripe)
library(pROC)
.libPaths("/home/ashika/myRpackages")

  
##########################################################################################
## Extract subject ids
##########################################################################################
in.dir <- '/project/Spatial_OASIS/processed/' 
in.dir_T1_GS <- '/project/Spatial_OASIS/Images/' 
out.dir <- '/project/Spatial_OASIS/Mimosa_data/'
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
                      pattern = "T2_stripped_reg[.]nii[.]gz",
                      full.names = TRUE)
FLAIR_files = list.files(path =
                           file.path(in.dir),
                         pattern = "FLAIR_stripped_reg[.]nii[.]gz",
                         full.names = TRUE)
PD_files = list.files(path =
                        file.path(in.dir),
                      pattern = "PD_stripped_reg[.]nii[.]gz",
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

# Create predictors

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

# Creating Brain mask for all subjects
filepaths$brainmask = NA
for (i in seq(nrow(filepaths))) {
  # Load files
  visit_id = filepaths$visit_id[i]
  fname = file.path(out.dir, 
                    paste0("brainmask_",
                           visit_id, ".nii.gz"))
  if (file.exists(fname)) {
    T1 = readnii(filepaths$T1[i])
    T2 = readnii(filepaths$T2[i])
    FLAIR = readnii(filepaths$FLAIR[i])
    PD = readnii(filepaths$PD[i])
   brain_mask = create_brain_mask(
      T1,
      T2,
      FLAIR,
      PD
    )
  #   Save brain mask to local working directory
    writenii(brain_mask, 
             filename = fname)
  }
  filepaths$brainmask[i] = fname
}

save(file=paste0(out.dir,"filepaths.rda"),filepaths)

##########################################################################################
## Manual WhiteStripe
##########################################################################################
load(file='/project/Spatial_OASIS/Mimosa_data/filepaths.rda')
for (i in 1:nrow(filepaths)){
  visit.id <- filepaths$visit_id[i]
  # Manual White Stripe
  FLAIR_1 = readnii(filepaths$FLAIR[i])
    T1_1 = readnii(filepaths$T1[i])
    T2_1 = readnii(filepaths$T2[i])
    PD_1 = readnii(filepaths$PD[i])
  WS_T1 = WhiteStripe::whitestripe(T1_1, type="T1", stripped = TRUE, verbose = T)
  WS_T2 = WhiteStripe::whitestripe(FLAIR_1, type="T2", stripped = TRUE, verbose = T)
  
  T1 = WhiteStripe::whitestripe_norm(T1_1, WS_T1$whitestripe.ind)
  FLAIR = WhiteStripe::whitestripe_norm(FLAIR_1, WS_T2$whitestripe.ind)
  T2 = WhiteStripe::whitestripe_norm(T2_1, WS_T2$whitestripe.ind)
  PD = WhiteStripe::whitestripe_norm(PD_1, WS_T2$whitestripe.ind)

 writeNIfTI(nim = T1,filename =paste0("/project/Spatial_OASIS/WhiteStripe/T1_", visit.id))
 writeNIfTI(nim = FLAIR,filename =paste0("/project/Spatial_OASIS/WhiteStripe/FLAIR_", visit.id))
 writeNIfTI(nim = T2,filename =paste0("/project/Spatial_OASIS/WhiteStripe/T2_", visit.id))
 writeNIfTI(nim = PD,filename =paste0("/project/Spatial_OASIS/WhiteStripe/PD_", visit.id))
}



##########################################################################################
## Randomly split the data into three sets for training, validation, and testing
## We have 98 subjects. Put 65 in each of training/validation, and 33 in testing.
##########################################################################################


set.seed(395)
sets <- split(unique_subjects, sample(rep(1:2,c(65,33))))

training <- sets$'1'
testing <- sets$'2'

filepaths_train <- filepaths[filepaths$visit_id %in% training,] 
filepaths_test <- filepaths[filepaths$visit_id %in% testing,] 


##########################################################################################
# Export subjectIDs for each dataset into a .txt file. Even though we set the seed to ensure
# we get the same randomization each time the code is run, let's have an external file with
# the IDs for safety. 
##########################################################################################
# write(training, "/project/Spatial_OASIS/Mimosa_data/train_IDs.txt", sep="\n")
# write(testing, "/project/Spatial_OASIS/Mimosa_data/test_IDs.txt", sep="\n")

##########################################################################################
## First, we have to loop through the subjects in the training set and train the model - DONE
##########################################################################################

# Create mimosa data - for ALL subjects

mimosa_df_list = vector(mode = "list",
                        length = nrow(filepaths))
names(mimosa_df_list) = filepaths$visit_id

for (i in seq(nrow(filepaths))) {
  # Load files
  T1_training = readnii(filepaths$T1[i])
  T2_training = readnii(filepaths$T2[i])
  FLAIR_training = readnii(filepaths$FLAIR[i])
  PD_training = readnii(filepaths$PD[i])
  gold_standard = readnii(filepaths$GS[i])
  brain_mask = readnii(filepaths$brainmask[i])
  
  # Manual White Stripe
  WS_T1 = WhiteStripe::whitestripe(T1_training, type="T1", stripped = TRUE, verbose = T)
  WS_T2 = WhiteStripe::whitestripe(FLAIR_training, type="T2", stripped = TRUE, verbose = T)
  
  T1 = WhiteStripe::whitestripe_norm(T1_training, WS_T1$whitestripe.ind)
  FLAIR = WhiteStripe::whitestripe_norm(FLAIR_training, WS_T2$whitestripe.ind)
  T2 = WhiteStripe::whitestripe_norm(T2_training, WS_T2$whitestripe.ind)
  PD = WhiteStripe::whitestripe_norm(PD_training, WS_T2$whitestripe.ind)
  # Obtain the mimosa predictor data.frame
  
  mimosa_df_list[[i]] = mimosa_data(
    brain_mask = brain_mask, 
    FLAIR = FLAIR, 
    T1 = T1,
    T2 = T2, 
    PD = PD, 
    tissue = FALSE, 
    gold_standard = gold_standard, 
    normalize = 'no', 
    cand_mask = NULL, 
    slices = NULL, 
    orientation = c("axial", "coronal", "sagittal"), 
    cores = 1, 
    verbose = TRUE)$mimosa_dataframe
}

# Turn list into a single data.frame which has all subjects predictor data.frames
mimosa_df = dplyr::bind_rows(mimosa_df_list, .id = "visit_id")

save(file=paste0(out.dir,"mimosa_df.rda"),mimosa_df)


##########################################################################################
## Train model for subjects in training set
##########################################################################################
                       
load(file='/project/Spatial_OASIS/Mimosa_data/mimosa_df.rda')

formula = gold_standard ~ FLAIR_10 * FLAIR + FLAIR_20 * FLAIR + 
  PD_10 * PD + PD_20 * PD + 
  T2_10 * T2 + T2_20 * T2 + 
  T1_10 * T1 + T1_20 * T1 +
  FLAIRonT1_intercepts + FLAIRonT2_intercepts + FLAIRonPD_intercepts +
  T1onT2_intercepts + T1onPD_intercepts + T2onPD_intercepts +
  T1onFLAIR_intercepts + T2onFLAIR_intercepts + PDonFLAIR_intercepts + 
  T2onT1_intercepts + PDonT1_intercepts + PDonT2_intercepts +
  FLAIRonT1_slopes + FLAIRonT2_slopes + FLAIRonPD_slopes +
  T1onT2_slopes + T1onPD_slopes + T2onPD_slopes +
  T1onFLAIR_slopes + T2onFLAIR_slopes + PDonFLAIR_slopes +
  T2onT1_slopes + PDonT1_slopes + PDonT2_slopes

mimosa_df_train <- mimosa_df[mimosa_df$visit_id %in% training,]

mimosa_model = mimosa_fit(mimosa_df_train, formula = formula)

save(file=paste0(out.dir,"mimosa_model.rda"), mimosa_model)

load(file='/project/Spatial_OASIS/Mimosa_data/mimosa_model.rda')

##########################################################################################
## Apply to test subjects, generate mimosa training + probability maps for each subject in testing set
##########################################################################################

filepaths_WhiteStripe = vector(mode = "list",
                        length = nrow(filepaths_test))
names(filepaths_WhiteStripe) = filepaths_test$visit_id

for (i in 1:length(testing)){
visit.id  <- testing[i]
  mimosa_testdata <- mimosa_data(
    brain_mask = readnii(filepaths_test$brainmask[i]),
    FLAIR = readnii(paste0("/project/Spatial_OASIS/WhiteStripe/FLAIR_",visit.id,".nii.gz")),
    T1 = readnii(paste0("/project/Spatial_OASIS/WhiteStripe/T1_",visit.id,".nii.gz")),
    T2 = readnii(paste0("/project/Spatial_OASIS/WhiteStripe/T2_",visit.id,".nii.gz")),
    PD = readnii(paste0("/project/Spatial_OASIS/WhiteStripe/PD_",visit.id,".nii.gz")),
    tissue = FALSE,
    gold_standard = NULL,
    normalize = 'no',
    cand_mask = NULL,
    slices = NULL,
    orientation = c("axial", "coronal", "sagittal"),
    cores = 1,
    verbose = TRUE)
    mimosa_testdata_df = mimosa_testdata$mimosa_dataframe
  mimosa_candidate_mask = mimosa_testdata$top_voxels
  predictions = predict(mimosa_model,
                        newdata = mimosa_testdata_df,
                        type = 'response')
  brain_mask = readnii(filepaths_test$brainmask[i])
  probability_map = niftiarr(brain_mask, 0)
  probability_map[mimosa_candidate_mask == 1] = predictions
  
  probability_map = fslsmooth(probability_map, 
                              sigma = 1.25,
                              mask = brain_mask, 
                              retimg = TRUE,
                              smooth_mask = TRUE,
                              reorient = TRUE)
  writeNIfTI(nim = probability_map,filename =paste0("/project/Spatial_OASIS/Mimosa_data/probability_map_", visit.id))
  save(probability_map,file ="probability_map.rda")
  # Generate Predicted lesion masks
  T1_files = list.files(path =
                      file.path('/project/Spatial_OASIS/WhiteStripe/'),
                      pattern = paste0("T1_",visit.id),
                      full.names = TRUE)
T2_files = list.files(path =
                      file.path('/project/Spatial_OASIS/WhiteStripe/'),
                      pattern = paste0("T2_",visit.id),
                      full.names = TRUE)
FLAIR_files = list.files(path =
                          file.path('/project/Spatial_OASIS/WhiteStripe/'),
                         pattern = paste0("FLAIR_",visit.id),
                         full.names = TRUE)
PD_files = list.files(path =
                      file.path('/project/Spatial_OASIS/WhiteStripe/'),
                      pattern = paste0("PD_",visit.id),
                      full.names = TRUE)
GS_files = list.files(path =
                        file.path(in.dir_T1_GS),
                      pattern = paste0(visit.id,"_lesion_mask[.]nii[.]gz"),
                      full.names = TRUE)                     

filepaths_WhiteStripe[[i]] = data.frame(T1 = T1_files, T2 = T2_files,
                       FLAIR = FLAIR_files, PD = PD_files, GS = GS_files,
                       stringsAsFactors = FALSE)

  mimosa_training = mimosa_training(
    brain_mask = filepaths_test$brainmask[i],
    FLAIR = filepaths_WhiteStripe[[i]]$FLAIR,
    T1 = filepaths_WhiteStripe[[i]]$T1,
    T2 = filepaths_WhiteStripe[[i]]$T2,
    PD = filepaths_WhiteStripe[[i]]$PD,
    tissue = FALSE, 
    gold_standard = filepaths_test$GS[i],
    normalize = 'no', 
    slices = NULL, 
    orientation = c("axial", "coronal", "sagittal"), 
    cores = 1, 
    verbose = TRUE, 
    outdir = NULL, 
    optimal_threshold = seq(0.25, 0.35, 0.01))
  
  threshold = mimosa_training$estimated_optimal_threshold
  segmentation_mask = probability_map > threshold

  
  L = list(top_voxels = mimosa_candidate_mask,
           predictions = predictions,
           probability_map = probability_map,
           mimosa_training = mimosa_training,
           threshold = threshold,
           segmentation_mask = segmentation_mask
  )
  out <- paste("test_", visit.id, "_df", sep="")
  assign(out,L)
  save(list = out,file=paste0(out.dir,"mimosa_test_pred_", visit.id, ".rda"))
}


##########################################################################################
## Evaluate performance 
##########################################################################################
for (i in 1:length(testing)){
  visit.id  <- testing[i]
  load(file= paste0("/project/Spatial_OASIS/Mimosa_data/mimosa_test_pred_", visit.id, ".rda"))
  gold_standard = readnii(filepaths_test$GS[i])
  brain_mask = readnii(filepaths_test$brainmask[i])
  prob_true_GS <- gold_standard[brain_mask==1]
  L <- get(paste("test_", visit.id, "_df", sep=""))
  probability_map <- L[["probability_map"]]
  writeNIfTI(nim = probability_map,filename =paste0("/project/Spatial_OASIS/Mimosa_data/probability_map_", visit.id))
  prob_true_predict <- probability_map[brain_mask==1]
  roc <- roc(prob_true_GS, prob_true_predict)
  out <- paste("test_", visit.id, "_roc", sep="")
  assign(out,roc)
  save(list = out,file=paste0("/project/Spatial_OASIS/Mimosa_data/mimosa_test_roc_", visit.id, ".rda"))
  
  #Plot
  png(file=paste0(out.dir,"roc_",visit.id,"_mimosa.png"))
  plot(roc)
  mtext(paste0("ROC Curve for Subject ", visit.id))
  dev.off()
}

