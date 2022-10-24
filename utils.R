library(texreg)
library(sjstats)
library(comprehenr)
library(glue)

# CCDF
ccdf <- function(x) {
  ecdf_fun <- ecdf(x)
  df <- data.frame(x = unique(x), y = 1 - sapply(unique(x), ecdf_fun))
  return(df)
}

# Make formula
make_formula <- function(dv, ev, controls) {
  frmla <-
    as.formula(paste(
      dv,
      paste(c(ev, controls), collapse = " + "), 
      sep = " ~ "
    ))
  return(frmla)
}

# Log-like transformation
asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

# Plot correlation matrix
get_corhm_plot <- function(cor_matrix){
  cor_matrix[lower.tri(cor_matrix)] <- NA
  cor_matrix_melted <- melt(cor_matrix, na.rm = TRUE)
  # Heatmap
  corhm_plot <- ggplot(cor_matrix_melted, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#1976D2", high = "#E53935", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Correlation") +
    theme_bw()+ # minimal theme
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 12, hjust = 1,  color = "black"),
          axis.text.y = element_text(size = 12, color = "black")) +
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = ifelse(value > 0 & value < 1, paste0("  ", value), value)), color = "black", size = 3.5) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = grid::unit(c(2, 3, 1, 2), "mm"),
      legend.justification = c(1, 0),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
      guides(fill = guide_colorbar(barwidth = 1, barheight = 18,
                    title.position = "right", title.hjust = 0.5))
  
  return(corhm_plot)
}


# Pretty summary statistics
descriptiveStatistics <- function(x, digits = 3,
                                  subset = c("mean", "median", "min", "max", "sd", "skew", "kurtosis"),
                                  filename = NULL) {
  m <- psych::describeBy(x, group = rep("", nrow(x)), mat = TRUE, digits = digits, fast = F)[, subset]
  
  # fix the concatenated "1" for the grouping in describeBy
  rownames(m) <- colnames(x)
  
  # fix kurtosis (as psych::describeBy is actuall computing the excess kurtosis)
  colnames(m) <- plyr::mapvalues(colnames(m), c("kurtosis"), c("excess_kurtosis"))
  
  if (!is.null(filename)) {
    col_names <- plyr::mapvalues(colnames(m),
                                 c("mean", "median", "min", "max", "sd", "skew", "excess_kurtosis"),
                                 c("Mean", "Median", "Min.", "Max", "Std. dev.", "Skewness", "Excess kurtosis"))
    showColumns(col_names)
    
    print(xtable::xtable(m, digits = digits),
          only.contents = TRUE, include.colnames = FALSE, booktabs = TRUE,
          file = filename, type = "latex")
  }
  
  return(m)
}


# Regression model coefficient names to TeX
coef_names_tex <- function(model, model_variables, display_variables, cut = NULL) {
  
  comb_vars <- expand.grid(model_variables, model_variables, stringsAsFactors = FALSE) %>% filter(Var1 != Var2)
  
  dict <- apply(comb_vars, 1, paste, collapse = ":") %>% enframe(., value = "value") %>% bind_cols(comb_vars) %>% 
    left_join(tibble(model_variables, display_variables) %>% 
               dplyr::select(model_variables, DV1 = display_variables), by = c("Var1" = "model_variables")) %>%
    left_join(tibble(model_variables, display_variables) %>% 
                dplyr::select(model_variables, DV2 = display_variables), by = c("Var2" = "model_variables")) %>%
    mutate(display_variables = paste0(DV1, " $\\times$ ", DV2)) %>%
    dplyr::select(model_variables = value, display_variables = display_variables) %>% 
    bind_rows(tibble(model_variables, display_variables), .)
  
  custom_names_display <- lapply(model, function(x) rownames(summary(x)$coefficients)) %>% 
    unlist() %>% unique() %>% enframe(.) %>% dplyr::select(-name)
  
  if(!is.null(cut)) {
    custom_names_display <- custom_names_display %>% mutate(value = gsub(cut, "", value)) 
  }
  
  custom_names_display <- custom_names_display %>%
    left_join(dict, by = c("value" = "model_variables")) %>% 
    mutate(display_variables = ifelse(is.na(display_variables), value, display_variables)) %>%
    pull(display_variables)
  
  return(custom_names_display)
}


marginal_effects_plot <- function(
  model,
  ev,
  dv_label="",
  ev_label="",
  plot_legend=TRUE,
  xlim_min=-3,
  xlim_max=7,
  x_breaks=seq(-3, 7, 1),
  x_labels=c("", expression("-2" ~ sigma), "", expression(mu), "", expression("2" ~ sigma), "", expression("4" ~ sigma), "", expression("6" ~ sigma), ""),
  ylim_min=NA,
  ylim_max=10
  ) {

  df_meplot <- ggpredict(model, c(ev, "hate_dummy")) %>%
    mutate(group = if_else(group==1, "Hateful", "Normal"))

  if(plot_legend){
    legend = c(0.65, 0.9)
    show = TRUE
  } else {
    legend = "none"
    show = FALSE
  }

  vals <- seq(xlim_min, xlim_max, 1)
  label_vals <- to_list(for(i in vals) if(i != 0) expression(i ~ sigma))

  meplot <- ggplot(df_meplot, aes(x = x, y = predicted, colour = group)) +
    geom_line(aes(group = group, color = group, linetype = group), size = 1.5, show.legend=show) +
    scale_linetype_discrete(labels= c("Hateful", "Normal")) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), linetype = 0, alpha = 0.2, show.legend=FALSE) +
    scale_colour_manual(values=c("#CE1126", "#003893", "#56B4E9"), labels= c("Hateful", "Normal")) +
    labs(y = dv_label, x = ev_label) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = c(xlim_min, xlim_max)
    ) +
    ylim(ylim_min, ylim_max) +
    theme(
      plot.margin = unit(c(3,1,0,1), "lines"),
      legend.position = legend
    )
  return(meplot)
}


marginal_effects_plot_dummy <- function(model, ev, dv_label="", ev_label="", plot_legend=FALSE, ylim_min=NA, ylim_max=10) {
  df_meplot <- ggpredict(model, c(ev, "hate_dummy")) %>%
    mutate(
      group = if_else(group==1, "Hateful", "Normal"),
      x = if_else(x==1, TRUE, FALSE)
    )

  if(plot_legend){
    legend = c(0.3, 0.9)
  } else {
    legend = "none"
  }

  meplot <- ggplot(df_meplot, aes(x = x, y = predicted, colour = group)) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = -0.5), fatten = 10, size = 1, shape = 15) + #, alpha = 0.5
    # geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1) +
    scale_colour_manual(values=c("#CE1126", "#003893", "#56B4E9"), labels= c("Hateful", "Normal")) +
    labs(y = dv_label, x = ev_label) +
    ylim(ylim_min, ylim_max) +
    theme(
      plot.margin = unit(c(3,1,0,1), "lines"),
      legend.position = legend
    )
  return(meplot)
}