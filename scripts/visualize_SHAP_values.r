# load packages
library(tidyverse)
library(here)
library(ggembl)
library(readxl)
library(patchwork)



dataset <- "Quinten"
# pc <- -4

substrate_annotations <- read_xlsx("/g/scb/zeller/karcher/cazy_gut_microbiome/scripts/../data/20230607_glycan_annotations_cleaned.xlsx")

# load models and clean up
modelPaths <- list.files(here('data', 'models'), pattern = ".rds", full.names = TRUE)
modelPaths <- modelPaths[str_detect(modelPaths, dataset)]
models <- map(modelPaths, \(x) {
    readRDS(x)
})
names(models) <- modelPaths
models <- enframe(models, name = 'raw_path', value = "model_object")
models <- models %>%
    mutate(path = str_split_fixed(raw_path, "/", n = 20)[, 9]) %>%
    mutate(fold = str_split_fixed(path, "__", n = 3)[, 2]) %>%
    mutate(fold = str_replace(fold, "fold_id_", "")) %>%
    mutate(fold = as.numeric(fold)) %>%
    mutate(resampling = str_split_fixed(path, "__", n = 4)[, 3]) %>%
    mutate(resampling = str_replace(resampling, "repeat_", "")) %>%
    mutate(resampling = as.numeric(resampling)) %>%
    mutate(model_type = str_split_fixed(path, "__", n = 6)[, 5]) %>%
    mutate(model_type = str_replace(model_type, ".rds", "")) %>%
    select(-raw_path) %>%
    mutate(beta_values = map(model_object, \(x) {
        lam_min_index <- x$model$index['1se', ]
        betas <- x$model$glmnet.fit$beta[, lam_min_index]
        return(as.data.frame(betas) %>% rownames_to_column('feature'))
    }))


# load shap values and clean up
shapPaths <- list.files(here("results", 'kernelshap_objects'), pattern = ".rds", full.names = TRUE)
shapPaths <- shapPaths[str_detect(shapPaths, dataset)]
shap <- map(shapPaths, \(x) {
    readRDS(x)
})
names(shap) <- shapPaths
shap <- enframe(shap, name = 'raw_path', value = "shap_values")
shap <- shap %>%
    mutate(path = str_split_fixed(raw_path, "/", n = 20)[, 9]) %>%
    mutate(fold = str_split_fixed(path, "__", n = 4)[, 2]) %>%
    mutate(fold = str_replace(fold, "fold_id_", "")) %>%
    mutate(fold = as.numeric(fold)) %>%
    mutate(resampling = str_split_fixed(path, "__", n = 4)[, 1]) %>%
    mutate(resampling = str_replace(resampling, "resamp_id_", "")) %>%
    mutate(resampling = as.numeric(resampling)) %>%
    mutate(on = str_split_fixed(path, "__", n = 5)[, 3]) %>%
    mutate(on = str_replace(on, "on_", "")) %>%
    mutate(model_type = str_split_fixed(path, "__", n = 5)[, 4]) %>%
    mutate(model_type = str_replace(model_type, ".rds", "")) %>%
    mutate(model_type = str_replace(model_type, "which_model_", "")) %>%
    select(-raw_path) %>%
    identity()

# Loada training/testing data
profilePaths <- list.files(here("data", 'fold_info'), pattern = ".tsv", full.names = TRUE)
profilePaths <- profilePaths[str_detect(profilePaths, dataset)]
profiles <- map(profilePaths, \(x) {
    read_tsv(x)
})
names(profiles) <- profilePaths
profiles <- enframe(profiles, name = 'raw_path', value = "profile")
profiles <- profiles %>%
    mutate(path = str_split_fixed(raw_path, "/", n = 20)[, 9]) %>%
    mutate(fold = str_split_fixed(path, "__", n = 4)[, 1]) %>%
    mutate(fold = str_replace(fold, ".*fold", "")) %>%
    mutate(fold = as.numeric(fold)) %>%
    mutate(resampling = str_split_fixed(path, "__", n = 4)[, 2]) %>%
    mutate(resampling = str_replace(resampling, "repeat_", "")) %>%
    mutate(resampling = str_replace(resampling, ".tsv", "")) %>%
    mutate(resampling = as.numeric(resampling)) %>%
    mutate(on = str_split_fixed(path, "__", n = 4)[, 1]) %>%
    mutate(on = str_replace(on, "_fold.*", "")) %>%
    mutate(on = case_when(
        on == "test_data" ~ "testing",
        on == "training_data" ~ "training"
    )) %>%
    select(-path, -raw_path)

shap <- shap %>%
    left_join(profiles, by = c("fold", "resampling", 'on'))


shap <- shap %>%
    mutate(shap_values_long = pmap(list(1:dim(shap)[1], shap_values, profile), \(i, x, p) {
        # tmp <- x$S$CRC %>% as.data.frame()
        tmp <- x$S$`1` %>% as.data.frame()
        tmp$sampleID <- p$sampleID
        tmp <- tmp %>%
            pivot_longer(-sampleID) %>%
            rename(
                feature = name,
                shap_value = value
            )
        return(tmp)
    }))

shap <- shap %>%
    mutate(on = ifelse(on == "training", "training\nset", "test\nset"))


# Boxplots of shap values
shap <- shap %>%
    select(fold, resampling, on, model_type, shap_values_long) %>%
    unnest(shap_values_long) %>%
    # I evalaute shap on training and testing
    # For testing folds, I get one shap value per model and resampling
    # for training folds, I get  4 shap values (in 5x cv) per model and resampling
    # In any case, take the median shap value for each sampleID
    group_by(sampleID, on, model_type, feature)

shap_tmp <- shap

# Understand mean/variance relationship of test-sample shap values over resampling rounds

lel <- shap_tmp %>%
    # filter(on == "test\nset") %>%
    filter(on == "training\nset") %>%
    group_by(sampleID, feature, model_type) %>%
    summarize(`mean(shap)\n(over resampled models)` = mean(shap_value), `var(shap)\n(over resampled models)` = var(shap_value))

ggsave(plot = ggplot(data = lel,
    aes(x = log10(abs(`mean(shap)\n(over resampled models)`)), y = log10(`var(shap)\n(over resampled models)`))) +
    # geom_point() +
    # geom_hex() +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
    facet_grid(. ~ model_type) +
    xlim(c(-8, -2)) +
    ylim(c(-17, -2)) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = 'solid') +
    geom_abline(intercept = -1, slope = 1, color = "black", linetype = 'dashed', alpha = 0.8) +
    geom_abline(intercept = -2, slope = 1, color = "black", linetype = 'dashed', alpha = 0.6) +
    geom_abline(intercept = -3, slope = 1, color = "black", linetype = 'dashed', alpha = 0.4) +
    geom_abline(intercept = -4, slope = 1, color = "black", linetype = 'dashed', alpha = 0.2) +
    geom_abline(intercept = -5, slope = 1, color = "black", linetype = 'dashed', alpha = 0.1) +
    ggtitle(str_c("Dataset: ", dataset, "\nSHAP values mean vs variance\nover resampling rounds")) +
    xlab("log10(abs(mean(shap)))\nover resampled models") +
    ylab("log10(var(shap))\nover resampled models") +
    theme_publication(), filename = here('plots', str_c(dataset, "_shap_mean_vs_var.pdf")), width = 5, height = 3)



# for (mt in c("RF", "lasso")) {
for (mt in c("lasso")) {
    shap_tmp2 <- shap_tmp %>% filter(model_type == mt)
    shap_tmp2 <- shap_tmp2 %>%
        inner_join(shap_tmp2 %>%
            group_by(feature, model_type) %>%
            summarize(n = mean(abs(shap_value))) %>%
            arrange(desc(n)) %>%
            head(10), by = c("feature", "model_type"))
    shap_tmp2$feature <- factor(shap_tmp2$feature, levels = shap_tmp2 %>%
        group_by(feature) %>%
        summarize(n = mean(abs(shap_value))) %>%
        arrange(desc(n)) %>%
        pull(feature))
    plot <- ggplot(data = shap_tmp2, aes(x = sampleID, y = shap_value, fill = resampling)) +
        theme_presentation() +
        facet_wrap(~feature, ncol = 2) +
        geom_boxplot() +
        geom_abline(intercept = 0, slope = 0, color = "red", linetype = 'dashed') +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        ylab("SHAP value") +
        xlab("Samples")
    ggsave(plot = plot, filename = here('plots', str_c(dataset, "__", mt, "_shap_values_boxplot_by_resampling.pdf")), width = 12, height = 8)
}


shap_tmp <- shap_tmp %>% summarize(shap_value = median(shap_value))
# for (mt in c("lasso", "RF")) {
for (mt in c("lasso")) {

    shap_tmp2 <- shap_tmp %>%
        filter(model_type == mt)
    shap_tmp2$feature <- factor(shap_tmp2$feature, levels = shap_tmp2 %>%
        group_by(feature) %>%
        summarize(n = mean(abs(shap_value))) %>%
        arrange(desc(n)) %>%
        pull(feature))

    shap_tmp2 <- shap_tmp2 %>%
        inner_join(shap_tmp2 %>%
            group_by(feature) %>%
            summarize(n = mean(abs(shap_value))) %>%
            arrange(desc(n)) %>%
            head(10), by = "feature")


    plot <- ggplot() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        # geom_boxplot(data = shap_tmp2, aes(x = feature, y = shap_value, fill = on), outlier.color = NA) +
        geom_point(data = shap_tmp2, aes(x = feature, y = shap_value, color = on), position = position_jitterdodge(), alpha = 0.3, size = 1) +
        geom_point(data = shap_tmp2 %>%
            filter(on == 'training\nset') %>%
            group_by(on, model_type, feature) %>%
            summarize(
                mea = mean(abs(shap_value))
            ), aes(x = feature, y = mea, group = on), position = position_nudge(0.25), color = 'black', shape = 3, show.legend = FALSE) +
        geom_point(data = shap_tmp2 %>%
            group_by(on, model_type, feature) %>%
            filter(on == 'test\nset') %>%
            summarize(
                mea = mean(abs(shap_value))
            ), aes(x = feature, y = mea, group = on), position = position_nudge(-0.25), color = 'black', shape = 3, show.legend = FALSE) +
        # geom_point(data = shap_tmp2 %>%
        #     filter(on == 'training\nset') %>%
        #     group_by(on, model_type, feature) %>%
        #     summarize(
        #         mea = median(abs(shap_value))
        #     ), aes(x = feature, y = mea, group = on), position = position_nudge(0.25), color = 'black', shape = 16, show.legend = FALSE) +
        # geom_point(data = shap_tmp2 %>%
        #     group_by(on, model_type, feature) %>%
        #     filter(on == 'test\nset') %>%
        #     summarize(
        #         mea = median(abs(shap_value))
        #     ), aes(x = feature, y = mea, group = on), position = position_nudge(-0.25), color = 'black', shape = 16, show.legend = FALSE) +
        theme_presentation() +
        coord_flip() +
        ggtitle(str_c("Dataset: ", dataset, "\nModel: ", mt, "\neach dot is a sample\nblack cross mean(abs(shap))\nsorted by mean(abs(shap)) shap")) +
        NULL

    shap_tmp3 <- shap_tmp %>%
        group_by(feature) %>%
        summarize(frac_pos = mean(shap_value > 0)) %>%
        arrange(desc(frac_pos))

    ggsave(plot = plot, filename = here('plots', str_c(mt, "__", dataset, '__shap_values_boxplot_', '.pdf')), width = 5, height = 5)
}


# For simplicity, let's move on with RF on testing data

# mt <- "RF"
mt <- "lasso"
more_plot_data <- shap_tmp %>%
    # filter(model_type == mt, on == "test\nset")
    filter(model_type == mt, on == "training\nset")

more_plot_data$feature <- factor(more_plot_data$feature, levels = more_plot_data %>%
    group_by(feature) %>%
    summarize(n = mean(abs(shap_value))) %>%
    arrange(desc(n)) %>%
    pull(feature))

more_plot_data <- more_plot_data %>%
    left_join(
        profiles %>%
            select(profile) %>%
            unnest() %>%
            distinct() %>%
            pivot_longer(-c(sampleID, Condition)) %>%
            rename(feature = name, feature_value = value), by = c('sampleID', "feature"))

more_plot_data <- more_plot_data %>%
    mutate(Condition = case_when(
        Condition == "After fasting" ~ "After fasting",
        Condition == "Before" ~ "Before fasting"
    ))

# more_plot_data$feature <- factor(more_plot_data$feature, levels = l)

# Get spearman cors between genus abundance and shap to pimp the mean(abs(shap)) summary metric
more_plot_data <- more_plot_data %>%
    left_join(more_plot_data %>%
        group_by(feature) %>%
        summarize(
            spearman = cor(shap_value, feature_value, method = "spearman")
        ), by = 'feature') %>%
    mutate(spearman_sign = ifelse(spearman > 0, 1, -1))


l1 <- more_plot_data %>%
    group_by(feature) %>%
    summarize(n = mean(abs(shap_value) * spearman_sign)) %>%
    arrange(desc(n)) %>%
    head(25) %>%
    pull(feature)
l2 <- more_plot_data %>%
    group_by(feature) %>%
    summarize(n = mean(abs(shap_value) * spearman_sign)) %>%
    arrange(desc(n)) %>%
    tail(25) %>%
    pull(feature)
more_plot_data <- more_plot_data %>%
    inner_join(data.frame(feature = c(l1, l2)), by = "feature")
more_plot_data$feature <- factor(more_plot_data$feature, levels = c(l1, l2))

plot <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    # geom_point(data = more_plot_data %>% filter(feature_value == pc), aes(x = feature, y = shap_value, shape = Condition), position = position_jitter(height = 0, width = 0.25), color = 'black', alpha = 0.35) +
    # geom_point(data = more_plot_data %>% filter(feature_value > pc), aes(x = feature, y = shap_value, shape = Condition, color = feature_value), position = position_jitter(height = 0, width = 0.25), alpha = 0.5) +
    geom_point(data = more_plot_data, aes(x = feature, y = shap_value, shape = Condition, color = feature_value), position = position_jitter(height = 0, width = 0.25), alpha = 0.5) +
    geom_point(data = more_plot_data %>%
        group_by(feature) %>%
        summarize(n = mean(abs(shap_value)) * spearman_sign), aes(x = feature, y = n), shape = 18, color = 'orange', size = 2.5, inherit.aes = F) +
    geom_point(data = more_plot_data %>%
        group_by(feature) %>%
        summarize(n = mean(abs(shap_value)) * spearman_sign), aes(x = feature, y = n), shape = 5, color = 'black', size = 2.5, inherit.aes = F) +
    theme_presentation() +
    coord_flip() +
    # ggtitle(str_c("Dataset: ", dataset, "\nModel: ", mt, "\neach dot is a sample\nblack dots samples\nwith feature == 0")) +
    scale_color_continuous(low = "blue", high = "yellow") +
    scale_shape_manual(values = c("After fasting" = 16, "Before fasting" = 1)) +
    NULL +
    ylab("SHAP value") +
    xlab("Genus")

ggsave(plot = plot, filename = here('plots', str_c(dataset, "_SHAP_vs_relAb.pdf")), width = 5, height = 5)

# for (f in c("Fusobacterium", "Parvimonas", "Peptostreptococcus", "Porphyromonas", "CAG.41", "Anaerostipes", "Eubacterium_G")) {
#     print(str_c("Processing feature: ", f))

#     tmp <- more_plot_data %>%
#         filter(feature == f)

#     plot <- ggplot() +
#         geom_hline(yintercept = 0) +
#         geom_point(data = tmp, aes(x = feature_value, y = shap_value, color = Condition), alpha = 0.3) +
#         theme_presentation() +
#         ggtitle(f)


#     ggsave(plot = plot, filename = here('plots', str_c(dataset, "__", mt, "_SHAP_vs_relAb_scatter_", f, ".pdf")), width = 3.75, height = 2.8)

# }

plot1 <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    # geom_boxplot(data = more_plot_data, aes(x = feature, y = shap_value, fill = Condition), alpha = 0.3) +
    # geom_point(data = more_plot_data %>% filter(feature_value == pc), aes(x = feature, y = shap_value, fill = Condition, shape = Condition), position = position_jitterdodge(jitter.width = 0.25), color = 'black', alpha = 0.35) +
    # geom_point(data = more_plot_data %>% filter(feature_value > pc), aes(x = feature, y = shap_value, fill = Condition, shape = Condition, color = feature_value), position = position_jitterdodge(jitter.width = 0.25), alpha = 0.5) +
    geom_point(data = more_plot_data %>%
        rename(`CAZy copy number` = feature_value) %>%
        mutate(`CAZy copy number` = `CAZy copy number` + 1), aes(x = feature, y = shap_value, fill = Condition, shape = Condition, color = `CAZy copy number`), position = position_jitterdodge(jitter.width = 0.2), alpha = 0.3) +
    geom_point(data = more_plot_data %>%
        group_by(feature) %>%
        summarize(n = mean(abs(shap_value)) * spearman_sign), aes(x = feature, y = n), shape = 18, color = 'orange', size = 2.5, inherit.aes = F) +
    geom_point(data = more_plot_data %>%
        group_by(feature) %>%
        summarize(n = mean(abs(shap_value)) * spearman_sign), aes(x = feature, y = n), shape = 5, color = 'black', size = 2.5, inherit.aes = F) +
    theme_presentation() +
    coord_flip() +
    scale_color_continuous(low = "blue", high = "yellow", trans = 'log10') +
    scale_shape_manual(values = c("After fasting" = 16, "Before fasting" = 1)) +
    ylab("SHAP value") +
    xlab("CAZy family")

# ggsave(plot = plot1, filename = here('plots', str_c(dataset, "_SHAP_vs_relAb_by_case_control.pdf")), width = 5, height = 8.5)

plot2 <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_boxplot(data = more_plot_data, aes(x = feature, y = shap_value, fill = Condition), alpha = 1) +
    theme_presentation() +
    coord_flip() +
    scale_color_continuous(low = "blue", high = "yellow", trans = 'log10') +
    scale_shape_manual(values = c("After fasting" = 16, "Before fasting" = 1)) +
    ylab("SHAP value") +
    xlab("CAZy family")

subtrate_info_for_plot_3 <- shap_tmp %>%
    group_by(feature) %>%
    summarize(n = mean(abs(shap_value))) %>%
    arrange(desc(n)) %>%
    left_join(substrate_annotations, by = c("feature" = "Subfamily"))

expand_subtrate_classes_into_matrix <- function(v, families) {
    # Prep
    v <- map(v, \(x) str_split(x, ",")[[1]])
    # Get all unique classes
    u <- unique(unlist(v))
    v <- map(v, \(x) data.frame(Seen = as.integer(u %in% x)))
    v <- do.call('cbind', v)
    u[is.na(u)] <- "NA"
    rownames(v) <- u
    colnames(v) <- families
    return(v)
}

data_for_plot_3 <- expand_subtrate_classes_into_matrix(subtrate_info_for_plot_3$FUNCTION_AT_DESTINATION_1, subtrate_info_for_plot_3$feature) %>%
    rownames_to_column('substrate_class') %>%
    pivot_longer(-substrate_class) %>%
    rename() %>%
    rename(
        family = name,
        Substrate = value) %>%
    mutate(Substrate = Substrate == 1) %>%
    inner_join(data.frame(feature = unique(more_plot_data$feature)), by = c('family' = 'feature')) %>%
    mutate(family = factor(family, levels = levels(more_plot_data$feature))) %>%
    mutate(substrate_class = factor(
        substrate_class,
        levels = c(unique(substrate_class)[!unique(substrate_class) %in% c("NA", "Other", "Unknown")], c("NA", "Other", "Unknown"))))

plot3 <- ggplot(data = data_for_plot_3) +
    geom_tile(aes(x = substrate_class, y = family, fill = Substrate)) +
    theme_presentation() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylab("CAZy family") +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
    NULL

ggsave(plot = plot1 + plot2 + plot3 + plot_layout(width = c(2, 1.75, 1.5), guides = 'collect'), filename = here('plots', str_c(dataset, "_SHAP_vs_relAb_by_case_control.pdf")), width = 10, height = 11)

# Adding code to compare SHAP values to single-feature wilcox values
data <- profiles %>%
    filter(resampling == 1) %>%
    unnest() %>%
    select(-fold, -resampling, -on) %>%
    select(-sampleID)

cond <- data$Condition
data$Condition <- NULL

wilcox_results <- list()
gfcs <- list()
for (fam in colnames(data)) {
    wilcox_results[[fam]] <- wilcox.test(data[[fam]] ~ cond)
    gfcs[[fam]] <- CalcGFC(data[cond == "Before", fam, drop = T], data[cond == "After fasting", fam, drop = T])
}

CalcGFC <- function(x.pos, x.neg, probs.fc = seq(.05, .95, .05)) {
    q.p <- quantile(x.pos, probs = probs.fc)
    q.n <- quantile(x.neg, probs = probs.fc)
    return(sum(q.p - q.n) / length(q.p))
}

enframe(wilcox_results) %>%
    rename(family = name, wilcox_test = value) %>%
    mutate(p = map_dbl(wilcox_test, \(x) x$p.value)) %>%
    mutate(gfcs = unlist(gfcs)) %>%
    arrange(p) %>%
    head(50) %>%
    select(family) %>%
    mutate(source = 'wilcox') %>%
    full_join(data.frame(family = levels(more_plot_data$feature)) %>% mutate(source = 'shap'), by = 'family', suffix = c("_wilcox", "_shap")) %>%
    mutate(overlap = case_when(
        !is.na(source_wilcox) & !is.na(source_shap) ~ "both",
        !is.na(source_wilcox) ~ "wilcox",
        !is.na(source_shap) ~ "shap"
    )) %>%
    group_by(overlap) %>%
    tally()

wilcox_results_fin <- enframe(wilcox_results) %>%
    rename(family = name, wilcox_test = value) %>%
    mutate(p = map_dbl(wilcox_test, \(x) x$p.value)) %>%
    mutate(gfcs = unlist(gfcs)) %>%
    arrange(p) %>%
    mutate(p.adj = p.adjust(p, method = "BH")) %>%
    filter(p.adj < 0.001) %>%
    mutate(positive = gfcs < 0) %>%
    group_by(positive) %>%
    nest() %>%
    mutate(data = map2(data, positive, function(da, po) {
        if (po) {
            return(da %>% arrange(gfcs) %>% head(25))
        } else {
            return(da %>% arrange(gfcs) %>% tail(25))
        }
    })) %>%
    unnest()
# .. get wilcox-based ranking.
data_for_plot_3 <- expand_subtrate_classes_into_matrix(subtrate_info_for_plot_3$FUNCTION_AT_DESTINATION_1, subtrate_info_for_plot_3$feature) %>%
    rownames_to_column('substrate_class') %>%
    pivot_longer(-substrate_class) %>%
    rename() %>%
    rename(
        family = name,
        Substrate = value) %>%
    mutate(Substrate = Substrate == 1) %>%
    # inner_join(data.frame(feature = unique(more_plot_data$feature)), by = c('family' = 'feature')) %>%
    inner_join(data.frame(feature = wilcox_results_fin$family), by = c('family' = 'feature')) %>%
    # mutate(family = factor(family, levels = levels(more_plot_data$feature))) %>%
    mutate(family = factor(family, levels = wilcox_results_fin$family)) %>%
    mutate(substrate_class = factor(
        substrate_class,
        levels = c(unique(substrate_class)[!unique(substrate_class) %in% c("NA", "Other", "Unknown")], c("NA", "Other", "Unknown"))))

plot3 <- ggplot(data = data_for_plot_3) +
    geom_tile(aes(x = substrate_class, y = family, fill = Substrate)) +
    theme_presentation() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ylab("CAZy family") +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey")) +
    NULL

ggsave(plot = plot3 + plot_layout(width = c(2, 1.75, 1.5), guides = 'collect'), filename = here('plots', str_c(dataset, "_SHAP_vs_relAb_by_case_control_wilcox_ranking_only_substrate_ranking_plot.pdf")), width = 5.5, height = 11)
