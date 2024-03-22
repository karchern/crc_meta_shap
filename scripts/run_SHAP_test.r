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
    on_what <- "training"
} else {
    args <- commandArgs(trailingOnly = TRUE)
    r_id <- as.integer(args[1])
    f_id <- as.integer(args[2])
    on_what <- args[3]
}

set.seed(11323)

# Load model
model <- readRDS(here(str_c("data/models/model__resamp_id_", r_id, "__fold_id_", f_id, ".rds")))
# Load profile
profiles_all <- read.table(here("data/otu_table.tsv"))
rownames(profiles_all) <- make.names(rownames(profiles_all))
# subset to training and testing fold
training_fold <- read_tsv(here('data/training_folds.tsv')) %>%
    filter(fold_id == f_id, resamp_id == r_id) %>%
    select(-fold_id, -resamp_id)
testing_fold <- read_tsv(here('data/testing_folds.tsv')) %>%
    filter(fold_id == f_id, resamp_id == r_id) %>%
    select(-fold_id, -resamp_id)
profiles_training <- profiles_all %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('sampleID') %>%
    inner_join(training_fold) %>%
    column_to_rownames('sampleID') %>%
    # t() %>%
    as.data.frame()
profiles_testing <- profiles_all %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('sampleID') %>%
    inner_join(testing_fold) %>%
    column_to_rownames('sampleID') %>%
    # t() %>%
    as.data.frame()

if (on_what == "training") {
    ps <- kernelshap(
        model, X = profiles_training, bg_X = profiles_training[sample(nrow(profiles_training), 250), ],
    )
} else if (on_what == "testing") {
    ps <- kernelshap(
        model, X = profiles_testing, bg_X = profiles_testing[sample(nrow(profiles_testing), 250), ],
    )
} else {
    stop("on_what must be either 'training' or 'testing'")
}

write_rds(ps, here(str_c("results/kernelshap_objects/resamp_id_", r_id, "__fold_id_", f_id, "__on_", on_what, ".rds")))

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
