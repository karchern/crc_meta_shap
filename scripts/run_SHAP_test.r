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
    r_id <- 3
    f_id <- 5
} else {
    args <- commandArgs(trailingOnly = TRUE)
    r_id <- as.integer(args[1])
    f_id <- as.integer(args[2])
}

set.seed(11323)

# profile <- read.table(here("data/Profiles.all.samples.tsv"), check.names = F)
# meta <- read.table(here('data/Metadata.all.samples.tsv'), check.names = F) %>%
#     rownames_to_column('sampleID') %>%
#     filter(str_detect(Cohort, "Zeller")) %>%
#     filter(Condition %in% c("CRC", "CTR"))

model <- readRDS(here(str_c("data/models/model__resamp_id_", resamp_id, "__fold_id_", fold_id, ".rds")))
profiles_all <- read.table(here("data/otu_table.tsv"))
training_fold <- read_tsv(here('data/training_folds.tsv')) %>%
    filter(fold_id == f_id, resamp_id == r_id) %>%
    select(-fold_id, -resamp_id)
testing_fold <- read_tsv(here('data/testing_folds.tsv')) %>%
    filter(fold_id == f_id, resamp_id == r_id) %>%
    select(-fold_id, -resamp_id)
profiles_testing <- profiles_all %>%
    # filter(sampleID %in% testing_fold$sampleID)
    inner_join()


# profile <- profile[, colnames(profile) %in% meta$sampleID]
# stopifnot(all(colnames(profile) == meta$sampleID))
# profile <- profile[apply(profile, 1, function(x) max(x) > 0.01 && mean(x > 0) > 0.2), ]
# profile <- log10(profile + 1E-5)
# profileWithMeta <- profile %>%
#     t() %>%
#     as.data.frame() %>%
#     rownames_to_column('sampleID') %>%
#     left_join(meta[, c("sampleID", "Condition")], by = 'sampleID') %>%
#     as.data.frame() %>%
#     column_to_rownames('sampleID')
# colnames(profileWithMeta) <- make.names(colnames(profileWithMeta))
# tsk <- as_task_classif(profileWithMeta, target = "Condition", positive = "CRC")

# split = partition(tsk, ratio = 0.8)
# learner = lrn("classif.ranger", predict_type = "prob")
# learner$train(tsk, row_ids = split$train)

# This runs for a bit
ps <- kernelshap(
    model, X = training_fold, bg_X = testing_fold,
)
ps_shapviz <- shapviz(ps)
ps_shapviz_importance <- sv_importance(ps_shapviz, kind = 'both', max_display = 30)
sv_waterfall(ps_shapviz, row_id = 1)

library(randomForest)
rfOOBImpoortance <- randomForest(
    as.formula(str_c("Condition ~ ", str_c(colnames(profileWithMeta)[colnames(profileWithMeta) != "Condition"], collapse = " + "))),
    # as.formula(str_c("effect_size ~ ", str_c(colnames(trainData)[1:10000], collapse = " + "))),
    data = profileWithMeta %>% mutate(Condition = factor(Condition, levels = c("CTR", "CRC"))),
    importance = TRUE
)$importance
