# Try from /g/scb2/zeller/karcher/mambaforge/envs/siamcat
library(mlr3)
library(SIAMCAT)
library(tidyverse)
library(here)

dataset <- "Quinten"

extractFold <- function(siamcatObject, what_fold = NULL) {
    if (what_fold == "training") {
        sc_folds <- siamcatObject@data_split$training.folds
    } else if (what_fold == "testing") {
        sc_folds <- siamcatObject@data_split$test.folds
    } else {
        stop("what_fold must be either 'training' or 'testing'")
    }
    all_folds <- map2(sc_folds, 1:length(sc_folds), \(resamp, resamp_id) {
        folds <- map2(resamp, 1:length(resamp), \(fold, fold_id) {
            fold <- data.frame(sampleID = fold)
            fold$fold_id <- fold_id
            return(fold)
        })
        folds <- do.call('rbind', folds)
        folds$resamp_id <- resamp_id
        return(folds)
    })
    all_folds <- do.call('rbind', all_folds)
    return(all_folds)
}

data_raw <- readRDS('/g/scb/zeller/karcher/SHAP/data/sc.obj_native_ML_blocked_cazy.RDS')
meta <- as.data.frame(data_raw@label$label)
colnames(meta)[1] <- 'Condition'
tmp <- as.data.frame(data_raw@label$info)
colnames(tmp)[1] <- "code"
tmp$trans <- rownames(tmp)
meta <- meta %>%
    rownames_to_column('sampleID') %>%
    inner_join(tmp, by = c("Condition" = "code")) %>%
    mutate(Condition = trans) %>%
    select(-trans)
norm_feat <- data_raw@norm_feat$norm.feat %>% t() %>% as.data.frame()
models_list <- data_raw@model_list$models
n <- names(models_list)
models_list <- map(models_list, \(x) x$model)
names(models_list) <- n

fold_ids <- extractFold(data_raw, what_fold = "testing")

for (foldIndex in 1:10) {
    for (repeatIndex in 1:10) {
        print(repeatIndex)
        # training_data <- profiles[cv$train_set(foldIndex), ]
        # test_data <- profiles[cv$test_set(foldIndex), ]
        test_data <- fold_ids %>%
            filter(resamp_id == repeatIndex) %>%
            filter(fold_id == foldIndex) %>%
            select(-fold_id, -resamp_id) %>%
            inner_join(
                norm_feat %>%
                    rownames_to_column('sampleID'),
                by = 'sampleID'
            ) %>%
            inner_join(meta, by = 'sampleID') %>%
            relocate(sampleID, Condition)

        write_tsv(test_data, here('data', 'fold_info', str_c(str_c("test_data_fold", foldIndex, "__repeat_", repeatIndex, "__", dataset, sep = ""), ".tsv")))
        # write_rds(rf, here('data/models', str_c("model__fold_id_", foldIndex, "__repeat_", repeatIndex, "__", dataset, "__RF.rds")))
        # write_rds(lasso, here('data/models', str_c("model__fold_id_", foldIndex, "__repeat_", repeatIndex, "__", dataset, "__lasso.rds")))
        write_rds(models_list[[str_c('cv_fold', foldIndex, "_rep", repeatIndex)]], here('data/models', str_c("model__fold_id_", foldIndex, "__repeat_", repeatIndex, "__", dataset, "__lasso.rds")))
    }
}
