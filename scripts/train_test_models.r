library(mlr3extralearners)
library(mlr3learners)
library(mlr3)
library(tidyverse)
library(here)
library(ggembl)

profiles <- read.table(here('data/Profiles.all.samples.tsv'), sep = "\t", check.names = F)
meta <- read.table(here('data/Metadata.all.samples.tsv'), sep = "\t", check.names = F)

# dataset <- "Zeller"
dataset <- "Feng"

meta <- meta %>% filter(Cohort == dataset)
# meta <- meta %>% filter(Cohort == "Feng")

profiles <- profiles[, colnames(profiles) %in% rownames(meta)]
profiles <- profiles[, match(rownames(meta), colnames(profiles))]
stopifnot(all(colnames(profiles) == rownames(meta)))

profiles <- profiles[apply(profiles, 1, \(x) max(x) > 0.0005), ]
profiles <- profiles[apply(profiles, 1, \(x) mean(x > 0) > 0.05), ]

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

# profiles$Condition <- as.factor(profiles$Condition)
profiles$Condition <- factor(profiles$Condition, levels = c("CTR", "CRC"))
colnames(profiles) <- make.names(colnames(profiles))

rf <- lrn("classif.ranger", predict_type = "prob", importance = "impurity")
lasso <- lrn("classif.cv_glmnet", alpha = 1, predict_type = 'prob')


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

set.seed(12321)
cv <- rsmp("repeated_cv", folds = numFolds, repeats = numRepeats)
set.seed(12321)
cv$instantiate(as_task_classif(profiles, target = "Condition"))

rf_aucs <- c()
lasso_aucs <- c()
rfs <- list()
lassos <- list()



for (foldRepeatIndex in 1:(numFolds * numRepeats)) {
    print(foldRepeatIndex)
    foldIndex <- cv$folds(foldRepeatIndex)
    repeatIndex <- cv$repeats(foldRepeatIndex)
    training_data <- profiles[cv$train_set(foldIndex), ]
    test_data <- profiles[cv$test_set(foldIndex), ]

    tasky <- as_task_classif(training_data, target = "Condition")
    rf$train(tasky)
    lasso$train(tasky)

    rf_aucs[length(rf_aucs) + 1] <- rf$predict_newdata(test_data)$score(msr('classif.auc'))
    lasso_aucs[length(lasso_aucs) + 1] <- lasso$predict_newdata(test_data)$score(msr('classif.auc'))

    rfs[[length(rfs) + 1]] <- rf
    lassos[[length(lassos) + 1]] <- lasso


    write_tsv(training_data %>%
        rownames_to_column('sampleID') %>%
        relocate('sampleID', "Condition"), here('data', 'fold_info', str_c(str_c("training_data_fold", foldIndex, "__repeat_", repeatIndex, "__", dataset, sep = ""), ".tsv")))
    write_tsv(test_data %>%
        rownames_to_column('sampleID') %>%
        relocate('sampleID', "Condition"), here('data', 'fold_info', str_c(str_c("test_data_fold", foldIndex, "__repeat_", repeatIndex, "__", dataset, sep = ""), ".tsv")))
    write_rds(rf, here('data/models', str_c("model__fold_id_", foldIndex, "__repeat_", repeatIndex, "__", dataset, "__RF.rds")))
    write_rds(lasso, here('data/models', str_c("model__fold_id_", foldIndex, "__repeat_", repeatIndex, "__", dataset, "__lasso.rds")))
}

plot <- ggplot(data = data.frame(aucs = c(rf_aucs, lasso_aucs), model = rep(c("RF", "Lasso"), each = length(rf_aucs)))) +
    geom_boxplot(aes(x = model, y = aucs)) +
    theme_publication()

ggsave(plot = plot, filename = here('plots', str_c('aucs_', dataset, '.pdf')), width = 3, height = 3, dpi = 300)
