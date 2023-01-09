##########################################################################################
## Project: Spatial OASIS
## File name: Preprocessing.R
## Author: Ashika Mani
## Date Created: 09/26/2022
## Last Modified: 10/28/2022
## Description:
##
##  This R program Preprocesses the data for the Spatial OASIS project
##
##########################################################################################

##########################################################################################
## Load the necessary libraries
##########################################################################################
library(fslr)
library(oasis)
library(oro.nifti)
library(abind)
library(ANTsR)
library(extrantsr)

##########################################################################################
## Set directory paths
##########################################################################################

in.dir <- '/project/Spatial_OASIS/Images/' 
out.dir <- '/project/Spatial_OASIS/processed/' 
setwd(out.dir)

##########################################################################################
## Extract subject ids
##########################################################################################
subjectDir <- dir(in.dir)
unique_subjects <- unique(substr(subjectDir, 1, 7))

##########################################################################################
## Skull strip T2 and PD, then register them to T1
##########################################################################################


for (i in 1:length(unique_subjects)) {
  subject.id <- unique_subjects[i]
  
  print(paste("Beginning: ",subject.id))
  
  # Read in ith subject's images
  T2 = readNIfTI( paste0(in.dir, subject.id, '_T2.nii.gz') , reorient=T)
  T1 = readNIfTI( paste0(in.dir, subject.id, '_T1.nii.gz') , reorient=T)
  PD = readNIfTI( paste0(in.dir, subject.id, '_PD.nii.gz') , reorient=T)
  FLAIR = readNIfTI( paste0(in.dir,subject.id, '_FLAIR.nii.gz') , reorient=T)
  
  # Skull strip T2  then register them to T1, create file 
  T2_stripped = fslbet(T2, reorient = T)
  T2_stripped_reg = registration(filename = T2_stripped, template.file = T1, correct = T, correction = "N4", typeofTransform = "Rigid", interpolator = "welchWindowedSinc")
  writeNIfTI(T2_stripped_reg$outfile, filename = paste0(subject.id, '_T2_stripped_reg'), gzipped=T)
  # 
  # # Skull strip PD  then register them to T1, create file 
  PD_stripped = fslbet(PD, reorient = T)
  PD_stripped_reg = registration(filename = PD_stripped, template.file = T1, correct = T, correction = "N4", typeofTransform = "Rigid", interpolator = "welchWindowedSinc")
  writeNIfTI(PD_stripped_reg$outfile, filename = paste0(subject.id, '_PD_stripped_reg'), gzipped =T)
  
  # Skull strip FLAIR  then register them to T1, create file 
 # Don't strip again here FLAIR_stripped = fslbet(FLAIR, reorient = T)
  FLAIR_stripped_reg = registration(filename = FLAIR, template.file = T1, correct = T, correction = "N4", typeofTransform = "Rigid", interpolator = "welchWindowedSinc")
  writeNIfTI(FLAIR_stripped_reg$outfile, filename = paste0(subject.id, '_FLAIR_stripped_reg'), gzipped =T)
}
