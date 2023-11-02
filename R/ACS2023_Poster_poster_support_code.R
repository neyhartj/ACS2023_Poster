## Poster support code
##
## Generate tables and figures for the poster
##

# Packages to load
library(tidyverse)
library(kableExtra)
library(neyhart)

# Set colors
heat_colors <- neyhart_palette("fall")[c(1,5)]

trait_abbr <- c("FruitYield" = "FY", "FruitWeight" = "FW", "Tacy" = "TAC", "TotalRottenFruitPercent" = "PFR")
trial_abbr <- c("ayt_2003" = "AYT03", "ayt_2013" = "AYT13")



main_dir <- neyhart::find_dir("ARS|CranberryLab")

# Project directories
proj_dir <- file.path(main_dir, "/Projects/CranberryPhenotypicBLUP")

# Directories
pheno_dir <- file.path(proj_dir, "data")
result_dir <- file.path(proj_dir, "output")
fig_dir <- file.path(getwd(), "figures/")


## Dimensions to help with producting figures
paper_width <- 34 # (inches)
paper_height <- 46

one_col <- 0.3013333 * paper_width
two_col <- 0.6266666 * paper_width
three_col <- 0.9279999 * paper_width







# Load results ------------------------------------------------------------

results_files <- list.files(result_dir, pattern = ".RData", full.names = TRUE)
results_objs <- sapply(results_files, function(x) load(file = x, envir = .GlobalEnv))





# Table 3. Stage 1 results -------------------------------------------------


# summarize stage 1 heritability and AIC comparing the best serial correlation model versus null
stage1_summ <- stage1b_model_fitting %>%
  imap_dfr(~{
    trl <- .y
    summ <- imap_dfr(.x, ~{
      trt <- .y
      map(.x[grepl("stage1", names(.x))], "fit") %>%
        imap_dfr(~mutate(.x, model = .y)) %>%
        mutate(trait = .y)
    }) %>%
      mutate(trial = trl)
  })

# Rescale AIC
stage1_summ2 <- stage1_summ %>%
  arrange(trial, trait, desc(model)) %>%
  split(list(.$trial, .$trait)) %>%
  map_df(~mutate(.x, AIC = AIC - AIC[1])) %>%
  mutate(model = ifelse(grepl("null", model), "null", "correlation")) %>%
  arrange(trial, trait, desc(model))

# Format for table creation
stage1_table <- stage1_summ2 %>%
  mutate(trial = toupper(trial),
         trial1 = paste0(str_sub(trial, 1, 3), str_sub(trial, -2, -1))) %>%
  select(Trial = trial1, Trait = trait, model, AIC) %>%
  spread(model, AIC) %>%
  select(Trial, Trait, Null = null, Correlation = correlation) %>%
  subset(grepl("AYT", Trial))


stage1_table %>%
  kable(x = ., format = "latex", align = "lccc", booktabs = TRUE) %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(c(" ", " ", "Model AIC" = 2))
  # column_spec(1, width = paste0(round(one_col * 0.3, 2), units)) %>%
  # column_spec(2, width = paste0(round(one_col * 0.3, 2), units)) %>%
  # column_spec(3, width = paste0(round(one_col * 0.3, 2), units))




# Figure 3. Stage 2 genetic covariance matrices ---------------------------


# Only do this for the ayt trials
stage2_model_fitting_use <- stage2_model_fitting_df %>%
  subset(grepl("ayt", trial)) %>%
  subset(!sapply(stage2_out, is.null)) %>%
  mutate(aic = map_dbl(stage2_out, "aic"))

# Remove null models for stage 1; then select the model with the lowest AIC
# per trial and trait
stage2_model_fitting_use1 <- stage2_model_fitting_use %>%
  subset(spatial_serial_structure == "best") %>%
  group_by(trial, trait) %>%
  top_n(n = 1, wt = -aic) %>%
  ungroup()

# Extract the G_h matrices
pheno_modeling_stage2_df <- stage2_model_fitting_use1 %>%
  mutate(G_h = map(stage2_out, ~as.matrix(.x$vars@geno1)))


# Convert the G matrices to data.frame
blup_model_Gmats <- pheno_modeling_stage2_df %>%
  select(trial, trait, model = gen_har_serial_structure, G_h) %>%
  # mutate(G_h = modify_if(G_h, trait == "FruitYield", ~fy_unit_convert(sqrt(.x), "Mg.ha")^2)) %>%
  mutate(G_h = imap(G_h, ~{
    # Convert to half cov - half cor
    covmat <- .x
    cormat <- cov2cor(covmat)
    # Make the upper triangle the correlation
    covmat[lower.tri(covmat)] <- cormat[lower.tri(cormat)]

    # Tidy up
    covmat %>%
      as.data.frame() %>%
      rownames_to_column("harvest_yap1") %>%
      gather(harvest_yap2, varcomp, -harvest_yap1)

  }))

# Create a plot for each G_h matrix
blup_model_Gmats_plots <- blup_model_Gmats %>%
  group_by(trial, trait, model) %>%
  do(plot = {
    row <- .

    dat <- row$G_h[[1]] %>%
      mutate_at(vars(contains("harvest_yap")), as.factor) %>%
      mutate(harvest_yap2 = fct_rev(harvest_yap2)) %>%
      # mutate(varcomp_text = format_numbers(varcomp)
      mutate(varcomp_text = formatC(varcomp, digits = 3))

    # Add indicator for upper triangle
    dummy_mat <- matrix(0, nlevels(dat$harvest_yap1), nlevels(dat$harvest_yap1),
                        dimnames = list(levels(dat$harvest_yap1), levels(dat$harvest_yap1)))
    dummy_mat[upper.tri(dummy_mat)] <- 1
    dummy_mat[lower.tri(dummy_mat)] <- 3
    diag(dummy_mat) <- 2

    dummy_mat_df <- dummy_mat %>%
      as.data.frame() %>%
      rownames_to_column("harvest_yap1") %>%
      gather(harvest_yap2, varcomp, -harvest_yap1) %>%
      mutate(position = c("upper", "diag", "lower")[varcomp]) %>%
      select(-varcomp)

    dat1 <- merge(x = dat, y = dummy_mat_df) %>%
      mutate(tile_fill = case_when(
        position == "lower" ~ varcomp,
        position == "diag" ~ as.numeric(NA),
        position == "upper" ~ 0
      ))

    trl <- row$trial
    trt <- row$trait
    mod <- row$model

    # Create the plot
    ggplot(dat1, aes(x = harvest_yap1, y = harvest_yap2, label = varcomp_text)) +
      geom_tile(aes(fill = tile_fill), color = "black") +
      geom_text(size = 2) +
      scale_x_discrete(position = "top") +
      scale_fill_gradient2(na.value = "grey85", low = heat_colors[1], high = heat_colors[2],
                           midpoint = 0, limits = c(-1, 1), guide = "none") +
      labs(subtitle = paste0("Trial: ", trial_abbr[trl], "; Trait: ", trait_abbr[trt], "; Model: ", toupper(mod))) +
      # theme_presentation2() +
      theme(axis.title = element_blank())

  }) %>% ungroup()



# Subset the matrix for the best model for the ayt2013 trial
gmats_ayt2003 <- blup_model_Gmats_plots %>%
  subset(trial == "ayt_2003")

# Create a combined imaged
# g_combined <- wrap_plots(gmats_ayt2013$plot)

g_combined <- gmats_ayt2003 %>%
  mutate(data = map(plot, "data")) %>%
  select(-plot) %>%
  unnest(data) %>%
  mutate(trait = factor(trait_abbr[trait], levels = trait_abbr),
         model = toupper(model)) %>%
  ggplot(aes(x = harvest_yap1, y = harvest_yap2, label = varcomp_text)) +
  geom_tile(aes(fill = tile_fill), color = "black") +
  geom_text(size = 2) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(na.value = "grey85", low = heat_colors[1], high = heat_colors[2], midpoint = 0,
                       limits = c(-1, 1), guide = "none") +
  facet_wrap(~ trait + model, ncol = 2, scales = "free", labeller = labeller(.multi_line = FALSE)) +
  # theme_presentation2() +
  theme(axis.title = element_blank(), strip.placement = "outside")


# Save the plot
ggsave(filename = "figure3.png", plot = g_combined, path = fig_dir,
       height = 4, width = 5, dpi = 300)




# Table of covariance structures ------------------------------------------

\begin{table}
\centering
\begin{tabular}{lc}
\toprule
\cmidrule(l{3pt}r{3pt})
Cov. Str. & Description \\
\midrule
ID & Constant genetic variance across years\\
CS & Constant genetic variance across years and constant correlation between years\\
CSH & Heterogenous genetic variance across years and constant correlation between years\\
AR1 & Constant genetic variance across years and decaying correlation between years\\
AR1H & Heterogenous genetic variance across years and decaying correlation between years\\
\bottomrule
\end{tabular}
\end{table}


cov_str_desc <- tribble(
  ~`Cov. Str.`, ~Description,
  "ID", "Constant genetic variance across years",
  "CS", "Constant genetic variance across years and constant correlation between years",
  "CSH", "Heterogenous genetic variance across years and constant correlation between years",
  "AR1", "Constant genetic variance across years and decaying correlation between years",
  "AR1H", "Heterogenous genetic variance across years and decaying correlation between years",
)

cov_str_desc %>%
  kable(x = ., format = "latex", align = paste0(c("l", rep("c", ncol(cov_str_desc) - 1)), collapse = ""),
        booktabs = TRUE) %>%
  kable_styling(full_width = FALSE)







# Figure 4. Cross-validation results --------------------------------------

# Display RMSE as a percentage of the worst model
future_year_crossval_df1 <- future_year_crossval_df %>%
  subset(spatial_serial_structure == "best") %>%
  unnest(crossv_out) %>%
  mutate(model = factor(gen_har_serial_structure, levels = c("idt", "corv", "corh", "ar1", "ar1h"))) %>%
  group_by(trial, trait, model) %>%
  summarize_at(vars(corr, rmse, avg_r2), mean, na.rm = TRUE)

future_year_crossval_df2 <- future_year_crossval_df1 %>%
  group_by(trial, trait) %>%
  mutate(rmse = rmse / max(rmse))


# Quick plot
g_crossval_figure <- future_year_crossval_df2 %>%
  mutate(trial = trial_abbr[trial],
         model = str_replace(model, "corh", "csh") %>% str_replace("corv", "cs"),
         model = fct_inorder(model)) %>%
  # subset(trial == "ayt_2003") %>%
  ggplot(aes(x = model, y = rmse * 100)) +
  geom_col() +
  scale_y_continuous(name = "RMSEP (% of the worst model)", breaks = pretty) +
  scale_x_discrete(name = expression(bold(G)~matrix~covariance~structure), labels = toupper) +
  facet_grid(trait ~ trial, scales = "free_y", switch = "y",
             labeller = labeller(trait = trait_abbr, trial = label_both)) +
  theme_poster()

ggsave(filename = "figure5.png", plot = g_crossval_figure, path = fig_dir, height = 7, width = 6)


# Create a table
crossval_table <- future_year_crossval_df1 %>%
  ungroup() %>%
  # subset(trial == "ayt_2003") %>%
  mutate(trait_abbr = trait_abbr[trait],
         trial_abbr = trial_abbr[trial],
         rmse_text = format_numbers(rmse)) %>%
  select(Trial = trial_abbr, Trait = trait_abbr, model, rmse_text) %>%
  spread(model, rmse_text)


crossval_table %>%
  kable(x = ., format = "latex", align = paste0(c("l", rep("c", ncol(crossval_table) - 1)), collapse = ""),
        booktabs = TRUE) %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(c(" ", " ", "RMSEP" = ncol(crossval_table) - 2))







# Figure list -------------------------------------------------------------

#
# 1. Figure 1. diagram of correlations that are expected / overhead example of cranberry
# trial
# 2. Figure 2. phenotypic correlations in different trials
# 3. Figure 3. stage 1 heritability / AIC
# 4. Figure 4. BLUP reliability / cross-validation accuracy
# 5. Figure 5 or table. Example of genotype rank change





# Figure X - BLUP reliabilities -------------------------------------------






# Plot BLUP reliabilities with different covariance structures

# Plot reliabilities of overall BLUPs
gen_blup <- blup_output_df %>%
  select(trial, trait, spatial_serial_structure, gen_har_serial_structure, gen_blup) %>%
  unnest(gen_blup)


gen_blup %>%
  ggplot(aes(x = gen_har_serial_structure, y = r2, fill = spatial_serial_structure)) +
  geom_boxplot() +
  facet_grid(trial ~ trait)













