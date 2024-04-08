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

rf <- lrn("classif.ranger", predict_type = "prob")
lasso <- lrn("classif.cv_glmnet", alpha = 1, predict_type = 'prob')

set.seed(12321)
numFolds <- 5
numRepeats <- 5

### This was initially working but kernelshap later fails since no task is saved in the objects if we fit them this wise
### Hence, manually extract training set and fit models that way.

# cv <-  rsmp("repeated_cv", folds = numFolds, repeats = numRepeats)
# cv$instantiate(task_data)

# rf_eval <- resample(task_data, rf, cv, store_models = TRUE, store_task = TRUE)
# rf_eval$aggregate(msr("classif.auc"))

# lasso_eval <- resample(task_data, lasso, cv, store_models = TRUE)
# lasso_eval$aggregate(msr("classif.auc"))

# rf_learners <- rf_eval$learners
# lasso_learners <- lasso_eval$learners

cv <-  rsmp("repeated_cv", folds = numFolds, repeats = numRepeats)
cv$instantiate(as_task_classif(profiles, target = "Condition"))

for (foldRepeatIndex in 1:(numFolds*numRepeats)){
    print(foldRepeatIndex)
    foldIndex <- cv$folds(foldRepeatIndex)
    repeatIndex <- cv$repeats(foldRepeatIndex)
    training_data <- profiles[cv$train_set(foldIndex), ]
    test_data <- profiles[cv$test_set(foldIndex), ]

    tasky <- as_task_classif(training_data, target = "Condition")
    rf$train(tasky)
    lasso$train(tasky)

    write_tsv(training_data %>% 
        rownames_to_column('sampleID') %>%
        relocate('sampleID', "Condition"), here('data', 'fold_info', str_c(str_c("training_data_fold", foldIndex, "__repeat_", repeatIndex, sep = ""), ".tsv")))
    write_tsv(test_data  %>% 
        rownames_to_column('sampleID') %>%
        relocate('sampleID', "Condition"), here('data', 'fold_info', str_c(str_c("test_data_fold", foldIndex, "__repeat_", repeatIndex, sep = ""), ".tsv")))
    write_rds(rf, here('data/models', str_c("model__fold_id_", foldIndex, "__repeat_", repeatIndex, "__RF.rds")))
    write_rds(lasso, here('data/models', str_c("model__fold_id_", foldIndex, "__repeat_", repeatIndex, "__lasso.rds")))
}


