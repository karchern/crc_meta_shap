library(mlr3extralearners)
library(mlr3learners)
library(mlr3)
library(tidyverse)
library(here)

profiles <- read.table(here('data/Profiles.all.samples.tsv'), sep = "\t", check.names = F)
meta <- read.table(here('data/Metadata.all.samples.tsv'), sep = "\t", check.names = F)

meta <- meta %>% filter(Cohort == "Zeller")

profiles <- profiles[, colnames(profiles) %in% rownames(meta)]
profiles <- profiles[, match(rownames(meta), colnames(profiles))]
stopifnot(all(colnames(profiles) == rownames(meta)))

profiles <- profiles[apply(profiles, 1, \(x) max(x) > 0.001), ]
profiles <- profiles[apply(profiles, 1, \(x) mean(x > 0) > 0.1), ]

profiles <- log10(profiles + 1E-4)

profiles <- profiles %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sampleID") %>%
    inner_join(meta %>%
        rownames_to_column("sampleID") %>%
        select(sampleID, Condition) %>%
        inner_join(data.frame(Condition = c("CRC", "CTR"))), by = "sampleID") %>%
    column_to_rownames("sampleID")

profiles$Condition <- as.factor(profiles$Condition)
colnames(profiles) <- make.names(colnames(profiles))

task_data <- as_task_classif(profiles, target = "Condition")

rf <- lrn("classif.randomForest", predict_type = "prob", importance = "accuracy")
lasso <- lrn("classif.cv_glmnet", alpha = 1, predict_type = 'prob')

set.seed(12321)
numFolds <- 5
cv <-  rsmp("repeated_cv", folds = numFolds, repeats = 5)
cv$instantiate(task_data)

rf_eval <- resample(task_data, rf, cv)
rf_eval$aggregate(msr("classif.auc"))

lasso_eval <- resample(task_data, lasso, cv)
lasso_eval$aggregate(msr("classif.auc"))

rf_learners <- rf_eval$learners
lasso_learners <- lasso_eval$learners

for (foldIndex in 1:numFolds){
    training_data <- profiles[cv$train_set(foldIndex), ] %>% select(-Condition)
    test_data <- profiles[cv$test_set(foldIndex), ]  %>% select(-Condition)
    write_tsv(training_data, here('data', 'fold_info', str_c(str_c("training_data_fold, ", foldIndex, sep = ""), ".tsv")))
    write_tsv(test_data, here('data', 'fold_info', str_c(str_c("test_data_fold, ", foldIndex, sep = ""), ".tsv")))
    write_rds(rf_learners, here('data/models', str_c("model__fold_id_", foldIndex, "__RF.rds")))
    write_rds(lasso_learners, here('data/models', str_c("model__fold_id_", foldIndex, "__lasso_learners.rds")))
}


