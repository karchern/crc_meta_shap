# load packages
library(tidyverse)
library(patchwork)
library(here)
library(ggembl)
# Loading the Matrix package seems crucial to avoid some weird error
library(Matrix)

vis_all <- function(dataset, label_case, model_types_to_evaluate) {
    #dataset <- "Selin20240604Balanced"
    ###
    #pc <- -4 # This needs to be manually set based on the dataset
    # pc (pseudo count) used to be mandatory and helpful for vis downstream if your profiles contain relative abundances
    # Selin's data seems to  be z-scored, so I'm getting rid of this functiooality
    ###
    #label_case <- '1' # This needs to be manually set based on the dataset
    #model_types_to_evaluate <- c("RF") # This needs to be manually set based on the dataset to produce plots in the beginning of this script

    # load models and clean up
    modelPaths <- list.files(here('data', 'models'), pattern = ".rds", full.names = TRUE)
    modelPaths <- modelPaths[str_detect(modelPaths, str_c(dataset, "__"))]
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
    shapPaths <- shapPaths[str_detect(shapPaths, str_c(dataset, ".rds"))]
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
    profilePaths <- profilePaths[str_detect(profilePaths, str_c(dataset, ".tsv"))]
    profiles <- map(profilePaths, \(x) {
        prof <- read.table(x, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>% as_tibble()
        colnames(prof) <- map_chr(colnames(prof), \(x) str_replace(x, "-", '.'))
        colnames(prof) <- ifelse(colnames(prof) == "51.20", "X51.20", colnames(prof))        
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
            tmp <- x$S[[label_case]] %>% as.data.frame()
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

    # Understand mean/variance relationship of test-sample shap values over resampling rounds

    lel <- shap_tmp %>%
        filter(on == "test\nset") %>%
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

    #This seems to take several minutes to run, which is a bummer.
    for (mt in model_types_to_evaluate) {
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
    for (mt in model_types_to_evaluate) {

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
    more_plot_data_all <- shap_tmp %>%
        filter(model_type == mt, on == "test\nset")

    more_plot_data_all$feature <- factor(more_plot_data_all$feature, levels = more_plot_data_all %>%
        group_by(feature) %>%
        summarize(n = mean(abs(shap_value))) %>%
        arrange(desc(n)) %>%
        pull(feature))
    l <- levels(more_plot_data_all$feature)        

    more_plot_data_all <- more_plot_data_all %>%
        left_join(
            profiles %>%
                select(profile) %>%
                unnest() %>%
                distinct() %>%
                pivot_longer(-c(sampleID, Condition)) %>%
                rename(feature = name, feature_value = value), by = c('sampleID', "feature"))

    # Get spearman cors between genus abundance and shap to pimp the mean(abs(shap)) summary metric
    more_plot_data_all <- more_plot_data_all %>%
        left_join(more_plot_data_all %>%
            group_by(feature) %>%
            summarize(
                spearman = cor(shap_value, feature_value, method = "spearman")
            ), by = 'feature') %>%
        mutate(spearman_sign = ifelse(spearman > 0, 1, -1))

    more_plot_data_all$feature <- factor(more_plot_data_all$feature, levels = rev(l))

    more_plot_data <- more_plot_data_all %>%
        inner_join(more_plot_data %>%
            group_by(feature) %>%
            summarize(n = mean(abs(shap_value))) %>%
            arrange(desc(n)) %>%
            head(10), by = "feature")


    plot <- ggplot() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        # geom_point(data = more_plot_data %>% filter(feature_value == pc), aes(x = feature, y = shap_value), position = position_jitter(height = 0, width = 0.25), color = 'black', alpha = 0.35) +
        # geom_point(data = more_plot_data %>% filter(feature_value > pc), aes(x = feature, y = shap_value, color = feature_value), position = position_jitter(height = 0, width = 0.25), alpha = 0.5) +
        geom_point(data = more_plot_data, aes(x = feature, y = shap_value, color = feature_value), position = position_jitter(height = 0, width = 0.45), alpha = 0.2, size = 1) +    
        geom_point(data = more_plot_data %>%
            group_by(feature) %>%
            summarize(n = mean(abs(shap_value)) * spearman_sign) %>%
            distinct(), aes(x = feature, y = n), shape = 18, color = 'orange', size = 2.5, inherit.aes = F) +
        geom_point(data = more_plot_data %>%
            group_by(feature) %>%
            summarize(n = mean(abs(shap_value)) * spearman_sign) %>%
            distinct(), aes(x = feature, y = n), shape = 5, color = 'black', size = 2.5, inherit.aes = F) +
        theme_presentation() +
        coord_flip() +
        # ggtitle(str_c("Dataset: ", dataset, "\nModel: ", mt, "\neach dot is a sample\nblack dots samples\nwith feature == 0")) +
        #scale_color_continuous(low = "blue", high = "yellow") +
        scale_color_gradientn(
            colors = c('#008000', "white",'#964B00'), 
            #values = c(-3,0,3),
            limits = c(-3, 3),
            oob = scales::squish) +
        # scale_color_gradient2() +
        scale_shape_manual(values = c("CRC" = 16, "CTR" = 1)) +
        NULL +
        ylab("SHAP value") +
        xlab("Genus")
    plot_right <- ggplot() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        geom_bar(data = more_plot_data %>%
            group_by(feature) %>%
            summarize(n = mean(abs(shap_value)) * spearman_sign) %>%
            distinct(), aes(y = n, x = feature), stat = 'identity') +
            theme_presentation() +
            coord_flip()   +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank()) +
            ylab("mean(abs(SHAP value))\n* spearman sign")
    ggsave(plot = plot + plot_right + plot_layout(widths = c(1, 0.8), guides = 'collect'), filename = here('plots', str_c(dataset, "_SHAP_vs_relAb.pdf")), width = 6, height = 5)

    for (f in c("Fusobacterium", "Parvimonas", "Peptostreptococcus", "Porphyromonas", "CAG.41", "Anaerostipes", "Eubacterium_G")) {
        print(str_c("Processing feature: ", f))

        tmp <- more_plot_data %>%
            filter(feature == f)

        plot <- ggplot() +
            geom_hline(yintercept = 0) +
            geom_point(data = tmp, aes(x = feature_value, y = shap_value, color = Condition), alpha = 0.3) +
            theme_presentation() +
            ggtitle(f)


        ggsave(plot = plot, filename = here('plots', str_c(dataset, "__", mt, "_SHAP_vs_relAb_scatter_", f, ".pdf")), width = 3.75, height = 2.8)

    }

    # Compare global shap values with single-feature wilcox test values
    shap <- more_plot_data_all %>%
            group_by(feature) %>%
            summarize(n = mean(abs(shap_value)) * spearman_sign) %>%
            distinct()
    wilcox <- profiles %>%
                select(profile) %>%
                unnest() %>%
                distinct() %>%
                pivot_longer(-c(sampleID, Condition)) %>%
                rename(feature = name, feature_value = value) %>%
                group_by(feature) %>%
                nest() %>%
                mutate(wilcox = map(data, \(x) {
                    wilcox.test(x %>% filter(Condition == "CRC") %>% pull(feature_value), x %>% filter(Condition == "CTR") %>% pull(feature_value))
                })) %>%
                mutate(wilcox_p = map_dbl(wilcox, "p.value"))
    shap_wilcox <- full_join(
        shap %>%
            rename(shap = n), 
        wilcox %>%
            select(feature, wilcox_p))



    # plot <- ggplot() +
    #     geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    #     #    geom_boxplot(data = more_plot_data, aes(x = feature, y = shap_value, fill = Condition), alpha = 0.3) +
    #     # geom_point(data = more_plot_data %>% filter(feature_value == pc), aes(x = feature, y = shap_value, fill = Condition, shape = Condition), position = position_jitterdodge(jitter.width = 0.25), color = 'black', alpha = 0.35) +
    #     # geom_point(data = more_plot_data %>% filter(feature_value > pc), aes(x = feature, y = shap_value, fill = Condition, shape = Condition, color = feature_value), position = position_jitterdodge(jitter.width = 0.25), alpha = 0.5) +
    #     geom_point(data = more_plot_data, aes(x = feature, y = shap_value, fill = Condition, shape = Condition, color = feature_value), position = position_jitterdodge(jitter.width = 0.25), alpha = 0.5) +
    #     geom_point(data = more_plot_data %>%
    #         group_by(feature) %>%
    #         summarize(n = mean(abs(shap_value)) * spearman_sign), aes(x = feature, y = n), shape = 18, color = 'orange', size = 2.5, inherit.aes = F) +
    #     geom_point(data = more_plot_data %>%
    #         group_by(feature) %>%
    #         summarize(n = mean(abs(shap_value)) * spearman_sign), aes(x = feature, y = n), shape = 5, color = 'black', size = 2.5, inherit.aes = F) +
    #     theme_presentation() +
    #     coord_flip() +
    #     scale_color_continuous(low = "blue", high = "yellow") +
    #     scale_shape_manual(values = c("CRC" = 16, "CTR" = 1)) +
    #     NULL +
    #     ylab("SHAP value") +
    #     xlab("Genus")

    # ggsave(plot = plot, filename = here('plots', str_c(dataset, "_SHAP_vs_relAb_by_CRC_CTR.pdf")), width = 5, height = 4)


    # models2 <- models %>%
    #     filter(model_type == "lasso") %>%
    #     select(fold, resampling, beta_values) %>%
    #     unnest() %>%
    #     rename(beta_value = betas) %>%
    #     mutate(beta_value = -1 * beta_value) %>%
    #     inner_join(
    #         shap %>%
    #             filter(model_type == 'lasso', on == "test\nset") %>%
    #             select(fold, resampling, on, model_type, shap_values_long) %>%
    #             unnest(shap_values_long) %>%
    #             # I evalaute shap on training and testing
    #             # For testing folds, I get one shap value per model and resampling
    #             # for training folds, I get  4 shap values (in 5x cv) per model and resampling
    #             # In any case, take the median shap value for each sampleID
    #             group_by(sampleID, on, model_type, feature) %>%
    #             group_by(fold, resampling, feature, on, model_type) %>%
    #             summarize(shap_value = median(shap_value)),
    #         by = c("fold", "resampling", "feature")) %>%
    #     arrange(desc(beta_value))

    # models2$feature <- factor(models2$feature, levels = models2 %>%
    #     group_by(feature) %>%
    #     summarize(n = median(beta_value)) %>%
    #     arrange(desc(n)) %>%
    #     pull(feature))

    # ggsave(plot = ggplot(
    #     data = models2 %>%
    #         filter(feature %in% levels(feature)[1:10])
    # ) +
    #     geom_boxplot(aes(x = feature, y = beta_value)) +
    #     geom_jitter(aes(x = feature, y = beta_value), alpha = 0.3) +
    #     theme_publication() +
    #     coord_flip() +
    #     xlab("Genus") +
    #     ylab("Lasso beta value")
    # , filename = here('plots', str_c(dataset, "_lasso_beta_vs_shap.pdf")), width = 3, height = 3)
    
}
