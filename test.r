# Based on https://mlr3book.mlr-org.com/chapters/chapter12/model_interpretation.html#sec-shapley
library(mlr3verse)
library(mlr3learners)
# library(iml) # provides Predictor object, providing different model-agnostic interpretation methods (SHAP but also others)
library(kernelshap) # provides fast approximation of shapley values via kernelshap function
library(shapviz)
library(tidyverse)

# Load Zeller_2014 profiles + metadata
profile <- read.table("/g/scb/zeller/karcher/SHAP_test/Profiles.all.samples.tsv", check.names = F)
meta <- read.table('/g/scb/zeller/karcher/SHAP_test/Metadata.all.samples.tsv', check.names = F) %>%
    rownames_to_column('sampleID') %>%
    filter(str_detect(Cohort, "Zeller")) %>%
    filter(Condition %in% c("CRC", "CTR"))
profile <- profile[, colnames(profile) %in% meta$sampleID]
stopifnot(all(colnames(profile) == meta$sampleID))
profile <- profile[apply(profile, 1, function(x) max(x) > 0.01 && mean(x > 0) > 0.2), ]
profile <- log10(profile + 1E-5)
profileWithMeta <- profile %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('sampleID') %>%
    left_join(meta[, c("sampleID", "Condition")], by = 'sampleID') %>%
    as.data.frame() %>%
    column_to_rownames('sampleID')
colnames(profileWithMeta) <- make.names(colnames(profileWithMeta))
tsk <- as_task_classif(profileWithMeta, target = "Condition", positive = "CRC")

split = partition(tsk, ratio = 0.8)
learner = lrn("classif.ranger", predict_type = "prob")
learner$train(tsk, row_ids = split$train)

# features in test data
credit_x = tsk$data(rows = split$test,
    cols = tsk$feature_names)
# target in test data
credit_y = tsk$data(rows = split$test,
    cols = tsk$target_names)

# This runs for a bit
ps <- kernelshap(
    learner, X = credit_x, bg_X = credit_x,
)
ps_shapviz <- shapviz(ps)
ps_shapviz_importance <- sv_importance(ps_shapviz, kind = 'both', max_display = 30)
sv_waterfall(ps_shapviz, row_id = 1)

library(randomForest)
rfOOBImpoortance <- randomForest(
    as.formula(str_c("Condition ~ ", str_c(colnames(profileWithMeta)[colnames(profileWithMeta) != "Condition"], collapse = " + "))),
    # as.formula(str_c("effect_size ~ ", str_c(colnames(trainData)[1:10000], collapse = " + "))),
    data = profileWithMeta  %>% mutate(Condition = factor(Condition, levels = c("CTR", "CRC"))),
    importance = TRUE
)$importance

shap <- apply(ps_shapviz$CRC$S, 2, mean) %>% as.data.frame()
colnames(shap)[1] <- "rf_oob"
shap <- shap %>%
    rownames_to_column('feature') 

rfOOBImpoortance <- rfOOBImpoortance %>%
    as.data.frame() %>%
    rownames_to_column('feature') %>%
    .[, c("feature", "MeanDecreaseAccuracy")]

inner_join(oob, rfOOBImpoortance) %>%
    ggplot(aes(x = rf_oob, y = MeanDecreaseAccuracy)) +
    geom_point() +
    theme_classic()
