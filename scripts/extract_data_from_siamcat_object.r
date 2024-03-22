library(SIAMCAT)
library(tidyverse)
library(here)

# args <- commandArgs(trailingOnly = TRUE)
args <- c('/g/scb/zeller/pekel/mi-eocrc/Results/Rdata/Training.all.data.rf.v5.Rdata')
modelPath <- args[1]

load(modelPath) # loads models.all.data.rf.v5


# x@data_split$training.folds does not contain expliciut information regarding fold. I assume it is simply correspond to 1:numFolds
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

# Extract model objects, input data and training/testing folds
load(modelPath)

mlr3_models <- map(models.all.data.rf.v5@model_list$models, \(model) model$model)
training_folds <- extractFold(models.all.data.rf.v5, "training")
testing_folds <- extractFold(models.all.data.rf.v5, "testing")
stopifnot(dim(inner_join(training_folds, testing_folds, by = c('sampleID', 'fold_id', 'resamp_id')))[1] == 0)

# Save extracted data
write_rds(mlr3_models, here('data/mlr3_models.rds'))
write_tsv(models.all.data.rf.v5@phyloseq@otu_table@.Data %>% as.data.frame(), here('data/otu_table.tsv'))
write_tsv(models.all.data.rf.v5@phyloseq@sam_data, here('data/metadata.tsv'))
write_tsv(training_folds, here('data/training_folds.tsv'))
write_tsv(testing_folds, here('data/testing_folds.tsv'))
