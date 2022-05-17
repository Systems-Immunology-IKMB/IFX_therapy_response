library(caret)
library(dplyr)
library(pROC)

setwd("/Users/nehamishra//Projects/SYSCID/IFX_therapy_response_random_forest/")

pat_data <- read.csv("Pat_status.csv")
cd_pat <- subset(pat_data, pat_data$Diagnosis == "CD")$Patient
uc_pat <- subset(pat_data, pat_data$Diagnosis == "UC")$Patient
ibd_pat <- pat_data$Patient

#expression_data <- read.csv("output/output_IBD/methylation_data/imputed_methylation_data.txt", sep = '\t')
expression_data <- read.csv("output/imputed_vst_data_emed_bfu_combined_02.txt", sep = '\t')
rownames(expression_data) <- expression_data$Patient
expression_data$Patient <- NULL

cd_data <- expression_data[as.character(cd_pat), ]
uc_data <- expression_data[as.character(uc_pat), ]
ibd_data <- expression_data[as.character(ibd_pat), ]

#selected_features <- read.csv("output/output_IBD/methylation_data/Feature_selection_janitza/Selected_features_from_DMP_janitza_120421.txt", header = TRUE, sep = '\t')
#selected_features <- read.csv("output/output_IBD/Feature_selection_janitza_ECG/Selected_features_from_ECG_janitza_120421.txt", header = TRUE, sep = '\t')
selected_features <- read.csv("output/output_UC/Baseline_module_genes/Selected_features_from_baseline_modules_separately_janitza_250920.txt", header = TRUE, sep = '\t')
selected_features <- subset(selected_features, selected_features$Module %in% c("plum","navajowhite2"))
selected_features <- subset(selected_features, !is.na(selected_features$Timepoint))
selected_features <- selected_features$Features

model_data <- uc_data[, c(as.character(selected_features), "Status")]
output_folder <- "output/output_UC/Model_training_testing_using_SCS_plum_navajowhite2///"

training <- model_data[!is.na(model_data$Status), ]

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE)

set.seed(248)

rfFit <- train(Status ~ ., data=training, 
               method="rf", preProc=c("center", "scale"), 
               trControl=fitControl, tuneLength=10)

rfFit
saveRDS(rfFit, file.path(output_folder, "expression_data_rf_model.rds"))

