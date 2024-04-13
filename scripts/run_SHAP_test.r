# Based on https://mlr3book.mlr-org.com/chapters/chapter12/model_interpretation.html#sec-shapley
library(mlr3verse)
library(mlr3learners)
# library(iml) # provides Predictor object, providing different model-agnostic interpretation methods (SHAP but also others)
library(kernelshap) # provides fast approximation of shapley values via kernelshap function
library(shapviz)
library(tidyverse)
library(here)


if (interactive()) {
    args <- list()
    r_id <- 1
    f_id <- 5
    on_what <- "testing"
    which_model <- "lasso"
    dataset <- 'Quinten'
} else {
    args <- commandArgs(trailingOnly = TRUE)
    r_id <- as.integer(args[1])
    f_id <- as.integer(args[2])
    on_what <- args[3]
    which_model <- args[4]
    dataset <- args[5]
}

set.seed(11323)

# Load model
lasso_model <- readRDS(here(str_c("data/models/model__fold_id_", f_id, "__repeat_", r_id, '__', dataset, "__lasso.rds")))
# RF_model <- readRDS(here(str_c("data/models/model__fold_id_", f_id, "__repeat_", r_id, '__', dataset, "__RF.rds")))

# models <- list("RF" = RF_model, "lasso" = lasso_model)
models <- list("lasso" = lasso_model)

# training_data_and_labels <- read_tsv(here('data', 'fold_info', str_c("training_data_fold", f_id, "__repeat_", r_id, '__', dataset, ".tsv")))
testing_data_and_labels <- read_tsv(here('data', 'fold_info', str_c("test_data_fold", f_id, "__repeat_", r_id, '__', dataset, ".tsv")))

# profiles_training <- training_data_and_labels %>% select(-Condition, -sampleID) %>% as.data.frame()
profiles_testing <- testing_data_and_labels %>% select(-Condition, -sampleID) %>% as.data.frame()
# training_labels <- training_data_and_labels %>% select(sampleID, Condition)
testing_labels <- testing_data_and_labels %>% select(sampleID, Condition)

if (on_what == "training") {
    ps <- kernelshap(
        models[[which_model]], X = profiles_training, bg_X = profiles_training,
    )
} else if (on_what == "testing") {
    ps <- kernelshap(
        models[[which_model]], X = profiles_testing, bg_X = profiles_testing,
    )
} else {
    stop("on_what must be either 'training' or 'testing'")
}

write_rds(ps, here(str_c("results/kernelshap_objects/resamp_id_", r_id, "__fold_id_", f_id, "__on_", on_what, "__which_model_", which_model, '__', dataset, ".rds")))

# ps_shapviz <- shapviz(ps)
# ps_shapviz_importance <- sv_importance(ps_shapviz, kind = 'both', max_display = 30)
# sv_waterfall(ps_shapviz, row_id = 1)

# library(randomForest)
# rfOOBImpoortance <- randomForest(
#     as.formula(str_c("Condition ~ ", str_c(colnames(profileWithMeta)[colnames(profileWithMeta) != "Condition"], collapse = " + "))),
#     # as.formula(str_c("effect_size ~ ", str_c(colnames(trainData)[1:10000], collapse = " + "))),
#     data = profileWithMeta %>% mutate(Condition = factor(Condition, levels = c("CTR", "CRC"))),
#     importance = TRUE
# )$importance
