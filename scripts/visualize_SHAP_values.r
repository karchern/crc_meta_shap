# load packages
library(tidyverse)
library(here)
library(ggembl)

dataset <- "Zeller"
pc <- -4

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
    mutate(model_type = str_split_fixed(path, "__", n = 5)[, 4]) %>%
    mutate(model_type = str_replace(model_type, ".rds", "")) %>%
    select(-raw_path)


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
        tmp <- x$S$CRC %>% as.data.frame()
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
shap_tmp <- shap %>%
    select(fold, resampling, on, model_type, shap_values_long) %>%
    unnest(shap_values_long) %>%
    # I evalaute shap on training and testing
    # For testing folds, I get one shap value per model and resampling
    # for training folds, I get  4 shap values (in 5x cv) per model and resampling
    # In any case, take the median shap value for each sampleID
    group_by(sampleID, on, model_type, feature)
shap_tmp <- shap_tmp %>% summarize(shap_value = median(shap_value))
for (mt in c("lasso", "RF")) {

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


mt <- "RF"
more_plot_data <- shap_tmp %>%
    filter(model_type == mt, on == "training\nset")

more_plot_data$feature <- factor(more_plot_data$feature, levels = more_plot_data %>%
    group_by(feature) %>%
    summarize(n = mean(abs(shap_value))) %>%
    arrange(desc(n)) %>%
    pull(feature))

more_plot_data <- more_plot_data %>%
    inner_join(more_plot_data %>%
        group_by(feature) %>%
        summarize(n = mean(abs(shap_value))) %>%
        arrange(desc(n)) %>%
        head(10), by = "feature")
l <- levels(more_plot_data$feature)

more_plot_data <- more_plot_data %>%
    left_join(
        profiles %>%
            select(profile) %>%
            unnest() %>%
            distinct() %>%
            pivot_longer(-c(sampleID, Condition)) %>%
            rename(feature = name, log10relAb = value), by = c('sampleID', "feature"))

more_plot_data$feature <- factor(more_plot_data$feature, levels = l)

plot <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_point(data = more_plot_data %>% filter(log10relAb == pc), aes(x = feature, y = shap_value, shape = Condition), position = position_jitter(height = 0, width = 0.25), color = 'black', alpha = 0.35) +
    geom_point(data = more_plot_data %>% filter(log10relAb > pc), aes(x = feature, y = shap_value, shape = Condition, color = log10relAb), position = position_jitter(height = 0, width = 0.25)) +
    geom_point(data = more_plot_data %>%
        group_by(feature) %>%
        summarize(n = mean(abs(shap_value))), aes(x = feature, y = n), shape = 18, color = 'orange', size = 2.5, inherit.aes = F) +
    geom_point(data = more_plot_data %>%
        group_by(feature) %>%
        summarize(n = mean(abs(shap_value))), aes(x = feature, y = n), shape = 5, color = 'black', size = 2.5, inherit.aes = F) +
    theme_presentation() +
    coord_flip() +
    ggtitle(str_c("Dataset: ", dataset, "\nModel: ", mt, "\neach dot is a sample\nblack dots samples\nwith feature == 0")) +
    scale_color_continuous(low = "blue", high = "red") +
    scale_shape_manual(values = c("CRC" = 16, "CTR" = 1)) +
    NULL +
    ylab("SHAP value") +
    xlab("Genus")

ggsave(plot = plot, filename = here('plots', str_c(dataset, "_SHAP_vs_relAb.pdf")), width = 5, height = 5)

for (f in c("Fusobacterium", "Parvimonas", "Peptostreptococcus", "Porphyromonas", "CAG.41", "Anaerostipes", "Eubacterium_G")) {
    print(str_c("Processing feature: ", f))

    tmp <- more_plot_data %>%
        filter(feature == f)

    plot <- ggplot() +
        geom_hline(yintercept = 0) +
        geom_point(data = tmp, aes(x = log10relAb, y = shap_value)) +
        theme_presentation() +
        ggtitle(f)


    ggsave(plot = plot, filename = here('plots', str_c(dataset, "__", mt, "_SHAP_vs_relAb_scatter_", f, ".pdf")), width = 2.5, height = 2.8)

}
