# Try from /g/scb2/zeller/karcher/mambaforge/envs/siamcat
library(mlr3)
library(SIAMCAT)
library(tidyverse)
library(here)

#dataset <- "Selin20240604AllData"
#dataset <- "Selin20240604Balanced"
#dataset <- "Selin20240604BalancedAndBlocked"

# For Selin's 3 types of models, this runs for around 10 minutes

model_type <- 'RF'
for (dataset in c("Selin20240604AllData", "Selin20240604Balanced", "Selin20240604BalancedAndBlocked")) {

    dataset_name_path_map <- list(
        "Selin20240604AllData" = '/g/scb/zeller/pekel/mi-eocrc/Results/Training.all.data.rf.Rdata',
        "Selin20240604Balanced" = '/g/scb/zeller/pekel/mi-eocrc/Results/Model.all.data.rf.balanced.Rdata',
        "Selin20240604BalancedAndBlocked" = '/g/scb/zeller/pekel/mi-eocrc/Results/Model.all.data.rf.balanced.blocked.Rdata'
    )

    dataset_name_model_object_map <- list(
        "Selin20240604AllData" = 'models.all.data.rf.v5',
        "Selin20240604Balanced" = 'models.all.data.rf.balanced',
        # TODO: The 2 lines below are remnant of the fact that Selin by mistake sent me the wrong model. Make sure this is correct for final iteration
        "Selin20240604BalancedAndBlocked" = 'models.all.data.rf.balanced.blocked'
        #"Selin20240604BalancedAndBlocked" = 'models.all.data.rf.balanced'
    )   

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

    load(dataset_name_path_map[[dataset]])
    assign("data_raw", get(dataset_name_model_object_map[[dataset]]), envir = .GlobalEnv)
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
    fold_ids_training <- extractFold(data_raw, what_fold = "training")

    for (foldIndex in 1:data_raw@data_split$num.folds) {
        for (repeatIndex in 1:data_raw@data_split$num.resample) {
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

            train_data <- fold_ids_training %>%
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
            write_tsv(train_data, here('data', 'fold_info', str_c(str_c("training_data_fold", foldIndex, "__repeat_", repeatIndex, "__", dataset, sep = ""), ".tsv")))
            write_rds(models_list[[str_c('cv_fold', foldIndex, "_rep", repeatIndex)]], here('data/models', str_c("model__fold_id_", foldIndex, "__repeat_", repeatIndex, "__", dataset, "__", model_type, ".rds")))
            }
        }
}

