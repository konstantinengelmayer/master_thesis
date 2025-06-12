library(dplR)
library(sf)
source("scripts/r/functions/extract_TRW.R")          # tree-ring widths
source("scripts/r/functions/extract_climate_variables.R")  

TRW_df <- extract_TRW(
  xls_folder  = "data/dendro_data/MW3_Daten/Dendro-Data/d-XLS",
  points_path = "data/vector_data/trees_all_plots.gpkg"
)

chronologies_by_plot <- TRW_df |>
  group_by(plot) |>
  group_modify(~{
    wide  <- pivot_wider(.x, Year, names_from = tree_id, values_from = ring_width)
    rwi   <- detrend(as.rwl(wide[,-1]), method = "Spline", nyrs = 13, f = 0.6)         # detrend
    chron <- chron(as.data.frame(rwi), prefix = paste0("CH_", unique(.x$plot)))
    chron |>
      tibble::rownames_to_column("year") |>
      mutate(year = as.numeric(year) + 1649,     # calendar offset
             plot = unique(.x$plot))
  }) |>
  ungroup() |>
  rename(TRI = starts_with("std"))               # TRI column



# 1. point layer with plot IDs
cores <- sf::read_sf("data/vector_data/trees_all_plots.gpkg")  # has column 'plot'

# 2. raster files (named!)
files <- c(
  tmax = "data/climate/dwd_4_variables/air_temp_max_1991_2024_merged.tif",
  pet  = "data/climate/dwd_4_variables/evapo_p_1991_2024_merged.tif",
  ppt  = "data/climate/dwd_4_variables/precipitation_1991_2024_merged.tif",
  soi  = "data/climate/dwd_4_variables/soi_moist_1991_2024_merged.tif"
)

# 3. run
plot_monthly_clim <- extract_plot_climate(
  points_sf = cores,
  var_files = files,
  buffer_m   = 1,
  years      = 2000:2022,
  fun        = mean,
  na.rm      = TRUE,
  start_date = as.Date("1991-01-01"),
  plot_col   = "plot"
)


## ---- monthly climate -------------------------------------------------
clim <- plot_monthly_clim %>%                   # has plot • year • month • vars …
  mutate(date = as.Date(sprintf("%04d-%02d-15", year, month)))

## ---- annual growth ---------------------------------------------------
growth <- chronologies_by_plot %>%
  select(plot, year, TRI) %>%
  dplyr::filter(!is.na(TRI))               # drop NaNs created by detrending

## convert the climate tibble to data.table
clim_dt <- as.data.table(clim)        # plot • year • month • ppt • tmax ...
clim_dt[, date := as.IDate(sprintf("%04d-%02d-15", year, month))]
setorder(clim_dt, plot, date)

growth_dt <- as.data.table(growth)    # plot • year • TRI

vars          <- c("ppt", "tmax", "pet", "soi")
max_window    <- 24

## make one long table: plot • date • var • value
clim_long <- melt(
  clim_dt,
  id.vars        = c("plot", "date", "year", "month"),
  measure.vars   = vars,
  variable.name  = "var",
  value.name     = "value",
  variable.factor= FALSE
)

## pre-allocate the wide columns cum1 … cum24  (much faster than lapply inside [.])
clim_long[, paste0("cum", 1:max_window) := 
            lapply(1:max_window, function(w)
              frollsum(value, w, align = "right")),  # NA padded on the left
          by = .(plot, var)]

## melt cum1 … cum24 back to long, add window_length and start_month
cum_long <- melt(
  clim_long,
  id.vars        = c("plot", "var", "year", "month"),
  measure.vars   = patterns("^cum"),
  variable.name  = "window_length",
  value.name     = "clim_sum",
  variable.factor= FALSE
)[
  , `:=`(
    window_length = as.integer(sub("cum", "", window_length)),
    start_month   = month                                 # month column
  )
][
  !is.na(clim_sum)                                         # drop leading NAs
]

## trim to years present in growth_dt to avoid extra NA joins
cum_long <- cum_long[year %in% growth_dt$year]

setkey(cum_long,   plot, year)
setkey(growth_dt,  plot, year)

merged <- cum_long[growth_dt, nomatch = 0L]   # inner join

## correlation per plot × var × window × start-month
corr_dt <- merged[
  , .(r = cor(clim_sum, TRI, use = "pair")),
  by = .(plot, var, start_month, window_length)
]

best <- corr_dt[
  order(-abs(r)),
  .SD[1],
  by = .(plot, var)
]


library(ggplot2)
library(RColorBrewer)

library(ggplot2)
library(RColorBrewer)
library(dplyr)

draw_heatmap <- function(df, plot_id) {
  df_plot <- df %>% dplyr::filter(plot == plot_id)
  
  ## 1. fixed 0.1-wide breaks spanning this plot’s range
  r_min    <- floor(min(df_plot$r, na.rm = TRUE) * 10) / 10
  r_max    <- ceiling(max(df_plot$r, na.rm = TRUE) * 10) / 10
  r_breaks <- seq(r_min, r_max, by = 0.1)
  
  ## 2. human-readable labels, already in ascending order
  r_labels <- paste0(sprintf("%.1f", head(r_breaks, -1)),
                     " – ",
                     sprintf("%.1f", tail(r_breaks, -1)))
  
  df_plot <- df_plot %>%
    mutate(
      r_bin = cut(r,
                  breaks = r_breaks,
                  include.lowest = TRUE,
                  labels = r_labels,
                  right = FALSE)                  # keep left-closed intervals
    )
  
  ## 3. palette: reverse Spectral so blue (low) → red (high)
  n_bins <- length(r_labels)
  pal    <- rev(colorRampPalette(brewer.pal(min(11, n_bins), "Spectral"))(n_bins))
  names(pal) <- r_labels                          # lock colours to labels
  
  ## 4. draw
  ggplot(df_plot,
         aes(factor(start_month), factor(window_length), fill = r_bin)) +
    geom_tile(colour = "white") +
    scale_fill_manual(
      values = pal,
      drop   = FALSE,
      name   = "r",
      guide  = guide_legend(reverse = TRUE)       # red at top of legend
    ) +
    facet_wrap(~ var, nrow = 1) +
    labs(
      title = paste("Rolling-sum correlation with TRI –", plot_id),
      x = "Start month",
      y = "Window length (months)"
    ) +
    theme_minimal() +
    theme(
      panel.grid    = element_blank(),
      axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.key.height = unit(0.6, "cm"),
      legend.text   = element_text(size = 7)
    )
}


plots <- unique(corr_dt$plot)

walk(plots, function(p) {
  g <- draw_heatmap(corr_dt, p)
  filename <- paste0("figures/climate_TRI_heatmaps_by_plot/corr_heatmap_", p, ".png")
  ggsave(filename, g, width = 10, height = 4, dpi = 300, bg = "white")
  message("saved: ", filename)
})


# Step 1: extract full info for the max |r| row per plot
strength_max_detailed <- corr_dt %>%
  group_by(plot) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  summarise(
    clim_connect_max    = max(abs(r), na.rm = TRUE),
    var_max             = var,
    start_month_max     = start_month,
    window_length_max   = window_length,
    r_max               = r,
    .groups = "drop"
  )




## ─────────────────────────────────────────────────────────────
## 2. mean of the top-N strongest windows  (N = 5)
## ─────────────────────────────────────────────────────────────
N <- 20
strength_meanTopN <- corr_dt %>%
  group_by(plot) %>%
  slice_max(abs(r), n = N, with_ties = FALSE) %>%
  summarise(clim_connect_meanTopN = mean(abs(r)),
            .groups = "drop")


## ─────────────────────────────────────────────────────────────
## 3. combined R² of a multi-window linear model
##    (best one window per climate variable)
## ─────────────────────────────────────────────────────────────
best_windows <- corr_dt %>%                       # pick strongest window per var
  group_by(plot, var) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup()

# pull those windows’ time series and join TRI
tri_clim <- best_windows %>%
  select(plot, var, start_month, window_length) %>%
  left_join(cum_long,
            by = c("plot", "var", "start_month", "window_length")) %>%   # adds year & clim_sum
  left_join(growth, by = c("plot", "year")) %>%                         # adds TRI
  dplyr::filter(!is.na(TRI))

# one linear model per plot
strength_r2 <- tri_clim %>%
  nest(data = -plot) %>%
  mutate(
    mod = map(data, ~ lm(TRI ~ clim_sum + 0, data = .x)),          # +0: no intercept
    clim_connect_r2 = map_dbl(mod, ~ summary(.x)$r.squared)
  ) %>%
  select(plot, clim_connect_r2)


## ─────────────────────────────────────────────────────────────
## 4. fraction of the surface with |r| ≥ 0.3
## ─────────────────────────────────────────────────────────────
thr <- 0.3
strength_area <- corr_dt %>%
  group_by(plot) %>%
  summarise(clim_connect_area = mean(abs(r) >= thr, na.rm = TRUE),
            .groups = "drop")


## ─────────────────────────────────────────────────────────────
## 5. put all four metrics side-by-side  (optional)
## ─────────────────────────────────────────────────────────────
strength_all <- strength_max %>%
  left_join(strength_meanTopN, by = "plot") %>%
  left_join(strength_r2,       by = "plot") %>%
  left_join(strength_area,     by = "plot")

writexl::write_xlsx(strength_all, "data/vector_data/plot_climate_growth_reaction.xlsx")
