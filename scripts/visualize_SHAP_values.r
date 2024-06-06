# load packages
library(tidyverse)
library(here)
library(ggembl)
# Loading the Matrix package seems crucial to avoid some weird error
library(Matrix)

source('/g/scb/zeller/karcher/SHAP/scripts/visualize_SHAP_values_utility.r')

for (dataset in c(
    "Selin20240604AllData"
    "Selin20240604Balanced",
    "Selin20240604BalancedAndBlocked"
    )) {
        vis_all(dataset = dataset, label_case = "1", model_types_to_evaluate = c("RF"))
    }

