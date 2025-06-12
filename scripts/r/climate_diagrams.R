# ─────────────────────────────────────────────────────────────
#  POINT-BASED CLIMATE & ELEVATION SUMMARY  •  FIVE AREAS (DE)
# ─────────────────────────────────────────────────────────────

# 1 ── libraries ──────────────────────────────────────────────
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(lubridate)

# 2 ── paths (edit if needed) ─────────────────────────────────
shp_file    <- "data/vector_data/all_merged.shp"               # polygons (just for names)
points_file <- "data/vector_data/trees_all_plots.gpkg"         # sampling points
raster_dir  <- "data/climate/dwd_4_variables/"
temp_file   <- file.path(raster_dir, "air_temp_mean_1991_2024_merged.tif")
prec_file   <- file.path(raster_dir, "precipitation_1991_2024_merged.tif")
dem_file    <- "data/satellite_data/predictors/dem.tif"

# 3 ── vectors & rasters ──────────────────────────────────────
areas  <- st_read(shp_file)
areas$area <- areas$name <- c("Eifel", "Kellerwald", "Koenigsforst",
                              "Calderner Wald", "Lindenberger Wald") 
points <- st_read(points_file) %>% select(-area)
  

Tmean <- rast(temp_file) / 10        # 0.1 °C → °C
Pmm   <- rast(prec_file)             # mm
dem   <- rast(dem_file)              # m

# label raster layers "YYYY-MM"
start_date <- ymd("1991-01-01")
names(Tmean) <- names(Pmm) <- format(seq(start_date, by = "1 month", length.out = nlyr(Tmean)), "%Y-%m")

# harmonise CRS and attach area name to each point
crs_target <- crs(Tmean)
points <- st_transform(points, crs_target) |>
  st_join(st_transform(areas, crs_target), left = FALSE)

# 4 ── extract raster values at points ────────────────────────
pts_v <- vect(points)
temp_df <- terra::extract(Tmean, pts_v, ID = FALSE) |> as_tibble()
prec_df <- terra::extract(Pmm,  pts_v, ID = FALSE) |> as_tibble()
elev_df <- terra::extract(dem,  pts_v, ID = FALSE) |> as_tibble() |> rename(elev_m = 1)

# tidy to long format (one row per point⋅month) ----------------
pt_meta <- points |> st_drop_geometry() |> mutate(pt_id = row_number()) |> select(pt_id, area)

temp_long <- temp_df |> mutate(pt_id = row_number()) |>
  pivot_longer(-pt_id, names_to = "date", values_to = "temp")

prec_long <- prec_df |> mutate(pt_id = row_number()) |>
  pivot_longer(-pt_id, names_to = "date", values_to = "prec")

clim_long <- left_join(temp_long, prec_long, by = c("pt_id", "date")) |>
  mutate(date = ymd(paste0(date, "-01")),
         month = month(date)) |>
  left_join(pt_meta, by = "pt_id")

# 5 ── summarise: annual precip + mean Apr-Sep temp per point --
veg_months <- 4:9  # April-September
summary_points <- clim_long |>
  group_by(area, pt_id, year = year(date)) |>
  summarise(
    prec_annual = sum(prec, na.rm = TRUE),
    temp_veg    = mean(temp[month %in% veg_months], na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(area) |>
  summarise(
    prec_mm    = mean(prec_annual),
    temp_veg_C = mean(temp_veg),
    elev_m     = mean(elev_df$elev_m),   # mean over points in the area
    .groups = "drop"
  ) |>
  mutate(across(c(prec_mm, temp_veg_C, elev_m), round, 1)) |>
  arrange(area)

print(summary_points)




library(dplyr)      # data wrangling
library(tidyr)      # pivot_* helpers
library(lubridate)  # date-time helpers (year(), yday(), …)
library(zoo)        # rollapply() for rolling sums
library(purrr)      # map_dfc(), walk()
library(ggplot2)    # all plotting
library(dplR)
library(terra)
library(sf)
library(zoo)        # for rolling means
library(signal)     # for Savitzky-Golay smoothing
library(readxl)
library(tibble)
library(RColorBrewer) 
library(data.table)
library(signal)
library(minpack.lm)


# custom extractors (paths relative to your repo)
source("scripts/r/functions/extract_TRW.R")          # tree-ring widths
source("scripts/r/functions/extract_SDC_buffered.R") # Landsat reflectance
source("scripts/r/functions/extract_spi_by_plot.R")
source("scripts/r/functions/extract_predictors.R")

# ─────────────────────────────────────────────────────────────────────────────
# 1. SATELLITE DATA  – daily VIs per plot
# ─────────────────────────────────────────────────────────────────────────────
## 1A  raw reflectance pixels --------------------------------------------------
SDC_df <- extract_SDC_points(
  folder_paths = c("data/satellite_data/SDC/lindenberg_eifel_koenigsforst",
                   "data/satellite_data/SDC/kellerwald_lahntal"),
  points_path  = "data/vector_data/trees_all_plots.gpkg",
  method = "original",
)
names(SDC_df)[1] <- "plot" 

evi_ts <- SDC_df %>%
  mutate(
    EVI2 =  ((nir - red)/(nir + red))*nir     # classic NDVI formula
  ) %>%
  select(plot, date, EVI2, area) %>%      # keep only what you need
  arrange(plot, date)

# alle Jahre ermitteln
all_years <- evi_ts %>% distinct(year = year(date)) %>% pull()

# ── 1) einfach anzeigen  ─────────────────────────────────────────────
walk(all_years, function(yr) {
  p <- evi_ts %>%
    dplyr::filter(year(date) == yr) %>%
    ggplot(aes(date, EVI2)) +
    geom_line() +
    facet_wrap(~ plot, scales = "free_y") +
    labs(title = paste("EVI2-Zeitreihe – Jahr", yr),
         x = NULL, y = "EVI2") +
    theme_minimal()
  
  print(p)           # zeigt Plot im Plot-Pane
})

# ────────────────────────────────────────────────────────────────────
#  Hilfsfunktionen
kv_part <- function(t, a, b, c) {            # Krümmungs-Ableitung
  z  <- exp(a + b*t);  bc <- b*c
  num1 <-  b^3 * c * z * (3*z*(1-z)*(1+z)^3)
  num2 <- (1+z)^2 * (1 + 2*z - 5*z^2)
  den  <- ((1+z)^4 + (bc*z)^2)^(5/2)
  z * (num1 - (bc*z)^3 * num2) / den
}

extrema <- function(kv) which(diff(sign(diff(kv))) != 0) + 1

# ════════════════════════════════════════════════════════════════════════
#  one_fit()  –  getrennte Logistik‐Fits für Frühling & Herbst
#               + SOS / EOS nach Zhang aus Krümmungs-Ableitung
#               + sofortiger Kontroll-Plot (Rohwerte + Fits + SOS/EOS)
#               ↦ liefert ein tibble (plot, year, SOS_doy, EOS_doy, conv)
# ════════════════════════════════════════════════════════════════════════

# ───────────────────────────────────────────────────────────────────
#  one_fit()  –  Doppel­logistik mit Grid-Search-Starts + nlsLM
# ───────────────────────────────────────────────────────────────────
one_fit <- function(df_one,
                    sg_p = 2, sg_n = 11,      # Glättung
                    amp_min = 0.05){          # min. Jahres-Amplitude
  
  out <- tibble(plot = first(df_one$plot),
                year = year(first(df_one$date)),
                SOS_doy = NA_real_, EOS_doy = NA_real_, conv = FALSE)
  if (nrow(df_one) < 15) return(out)
  
  # 3.1  Glättung ------------------------------------------------------------
  df_one <- df_one %>%
    arrange(date) %>%
    mutate(doy  = yday(date),
           EVI2 = signal::sgolayfilt(EVI2, p = sg_p, n = sg_n)) %>%
    drop_na()
  
  if (diff(range(df_one$EVI2)) < amp_min) return(out)
  
  # 3.2  feste Parameter -----------------------------------------------------
  d_fix  <- min(df_one$EVI2)
  c1_fix <- max(df_one$EVI2) - d_fix           # Frühling (+)
  c2_fix <- -c1_fix                            # Herbst   (−)
  
  # 3.3  Zeitachse skalieren  (t ≈ −6 … +6) ---------------------------------
  df_one <- df_one %>%
    mutate(t = (doy - 182.5) / 30)             # 30 d ≈ 1 t-Einheit
  
  # 3.4  Brute-Force Raster  -------------------------------------------------
  grid <- expand.grid(
    a1 = seq(-6,  0, 1),                       # Frühling-Lage
    b1 = seq(-0.20, -0.03, 0.03),              # Frühling-Steilheit
    a2 = seq( 1,  6, 1),                       # Herbst-Lage  (später)
    b2 = seq( 0.03, 0.20, 0.03)                # Herbst-Steilheit
  )
  
  y <- df_one$EVI2;  t <- df_one$t
  rss <- apply(grid, 1, function(g){
    g <- as.numeric(g)
    y_hat <- d_fix +
      c1_fix/(1+exp(g[1] + g[2]*t)) +
      c2_fix/(1+exp(g[3] + g[4]*t))
    sum((y - y_hat)^2)
  })
  start_best <- as.list(grid[which.min(rss), ])
  
  # 3.5  Fein-Fit  mit nlsLM -------------------------------------------------
  lower <- c(a1=-Inf, b1=-0.25, a2=-Inf, b2= 0.01)
  upper <- c(a1= Inf, b1=-0.01, a2= Inf, b2= 0.25)
  
  form <- EVI2 ~ d_fix +
    c1_fix/(1+exp(a1 + b1*t)) +
    c2_fix/(1+exp(a2 + b2*t))
  
  fit <- try(
    nlsLM(form, data = df_one,
          start = start_best,
          control = nls.lm.control(maxiter = 800)),
    silent = TRUE)
  
  if (inherits(fit,"try-error")) return(out)
  out$conv <- TRUE
  pa <- coef(fit)
  
  # 3.6  SOS / EOS nach Zhang-Krümmung --------------------------------------
  t_grid <- seq(-6, 6, 0.1)
  kv_s   <- kv_part(t_grid, pa["a1"], pa["b1"],  c1_fix)
  kv_a   <- kv_part(t_grid, pa["a2"], pa["b2"],  c2_fix)
  sp_e   <- extrema(kv_s)
  au_e   <- extrema(kv_a)
  
  # → DOY für jedes Extrem
  doy_sp <- t_grid[sp_e] * 30 + 182.5
  doy_au <- t_grid[au_e] * 30 + 182.5
  
  if (!length(doy_sp) || !length(doy_au)) return(out)
  
  if (mean(doy_au) < mean(doy_sp)) {
    tmp_doy <- doy_sp;  doy_sp <- doy_au;  doy_au <- tmp_doy   # DOY tauschen
    tmp_e   <- sp_e;    sp_e   <- au_e;    au_e   <- tmp_e     # Index-Vekt. tauschen
  }
  
  sos_range <- c( 90, 160)      # Frühling
  eos_range <- c(240, 330)      # Herbst
  
  ## ----------  SPRING  ------------------------------------------------------
  inside_sp <- which(doy_sp >= sos_range[1] & doy_sp <= sos_range[2])
  
  if (length(inside_sp)) {                       # mindestens eins im Fenster
    SOS_t <- t_grid[ sp_e[ inside_sp[1] ] ]      # frühestes nehmen
  } else {                                       # alle draußen → kleinstes Δ
    dist_sp <- pmax(sos_range[1]-doy_sp, 0) +
      pmax(doy_sp-sos_range[2], 0)      # Distanz zu Fenster
    SOS_t <- t_grid[ sp_e[ which.min(dist_sp) ] ]
  }
  
  ## ----------  AUTUMN  ------------------------------------------------------
  inside_au <- which(doy_au >= eos_range[1] & doy_au <= eos_range[2])
  
  if (length(inside_au)) {                       # mindestens eins im Fenster
    EOS_t <- t_grid[ au_e[ tail(inside_au, 1) ] ]  # spätestes nehmen
  } else {                                       # alle draußen → kleinstes Δ
    dist_au <- pmax(eos_range[1]-doy_au, 0) +
      pmax(doy_au-eos_range[2], 0)
    EOS_t <- t_grid[ au_e[ which.min(dist_au) ] ]
  }
  
  # ------------------------------------------------------------------
  #  **Sicherheits-Swap**: immer das frühere Datum als SOS ausgeben
  # ------------------------------------------------------------------
  
  SOS <- SOS_t*30 + 182.5
  EOS <- EOS_t*30 + 182.5
  
  if (EOS >= SOS + 60) {
    out$SOS_doy <- round(SOS)
    out$EOS_doy <- round(EOS)
  }
  out
}



# ────────────────────────────────────────────────────────────────────
#  Haupt-Pipeline  –  Loop über alle Plot-Jahre
phenology_df <- evi_ts %>%
  mutate(doy  = yday(date),
         year = year(date)) %>%
  group_by(plot, year) %>%
  group_split() %>%
  map_dfr(one_fit)

df_one <- evi_ts %>%
  dplyr::filter(plot == "SF05", year(date)== 2019) %>%
  mutate(doy = yday(date))
  

# ────────────────────────────────────────────────────────────────────
#  Ergebnis prüfen
print(head(phenology_df, 12))
table(phenology_df$conv)
# optional: speichern

pred <- tibble(
  doy = t_grid,
  fit = d_fix +
    c1_fix / (1 + exp(coef(fit)["a1"] + coef(fit)["b1"]*t_grid)) +
    c2_fix / (1 + exp(coef(fit)["a2"] + coef(fit)["b2"]*t_grid))
)

## 2. Rohpunkte + Modell in einem Plot
ggplot() +
  geom_point(data = df_one, aes(((doy - 182.5) / 30) , EVI2), colour = "grey40") +   # gemessene EVI2
  geom_line (data = pred,   aes(doy, fit),  colour = "red", size = 1) +  # Logistic-Fit
  labs(title = paste0("Plot ", unique(df_one$plot),
                      " – Jahr ", unique(lubridate::year(df_one$date))),
       x = "Tag des Jahres", y = "EVI2") +
  theme_minimal()



## ------------------------------------------------------------
## 1. Daten ins Long-Format bringen  (SOS / EOS untereinander)
## ------------------------------------------------------------
phen_long <- phenology_df %>%          # Ergebnis-Tabelle aus dem Fit
  dplyr::filter(conv) %>%                     # nur erfolgreiche Fits
  pivot_longer(cols = c(SOS_doy, EOS_doy),
               names_to  = "phase",
               values_to = "doy")

# phase = "SOS_doy" oder "EOS_doy"

## ------------------------------------------------------------
## Variante A  –  Boxplot pro Phase über ALLE Plots & Jahre
## ------------------------------------------------------------
ggplot(phen_long, aes(x = phase, y = doy, fill = phase)) +
  geom_boxplot(width = 0.5, alpha = .6, outlier.shape = 21) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Verteilung von Start- und End-of-Season (alle Plots, alle Jahre)",
       x = NULL, y = "Tag des Jahres") +
  theme_minimal() +
  theme(legend.position = "none")

## ------------------------------------------------------------
## Variante B  –  Facetten pro Plot, farbcodiert nach Phase
## ------------------------------------------------------------
ggplot(phen_long, aes(x = phase, y = doy, fill = phase)) +
  geom_boxplot(width = 0.5, alpha = .6, outlier.shape = 21) +
  facet_wrap(~ plot, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "SOS / EOS – Boxplots getrennt nach Plot",
       x = NULL, y = "Tag des Jahres") +
  theme_minimal() +
  theme(legend.position = "none")

## ------------------------------------------------------------
## Variante C  –  Jahres-Entwicklung (Boxplot je Jahr, Farbe = Phase)
## ------------------------------------------------------------
ggplot(phen_long, aes(x = factor(year), y = doy, fill = phase)) +
  geom_boxplot(outlier.shape = 21, position = position_dodge(width = .7)) +
  scale_fill_brewer(palette = "Set2", name = NULL,
                    labels = c("EOS", "SOS")) +
  labs(title = "SOS und EOS pro Jahr (alle Plots zusammen)",
       x = "Jahr", y = "Tag des Jahres") +
  theme_minimal()



