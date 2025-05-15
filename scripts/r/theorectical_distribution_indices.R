# ===============================================================
#  Histogram grid for 19 vegetation / moisture / soil indices
#  using random reflectance samples in the 0–1 physical range
# ===============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------
# 1)  Generate random reflectance inputs ---------------------------
# -----------------------------------------------------------------
set.seed(20250424)                   # reproducibility
n_sample <- 5e5                      # try 1e6 if you have plenty of RAM/CPU

raw <- tibble(
  nir   = runif(n_sample, 0.01, 0.99),
  red   = runif(n_sample, 0.01, 0.99),
  green = runif(n_sample, 0.01, 0.99),
  blue  = runif(n_sample, 0.01, 0.99),
  swir1 = runif(n_sample, 0.01, 0.99),
  swir2 = runif(n_sample, 0.01, 0.99)
)

# -----------------------------------------------------------------
# 2)  Compute all indices ------------------------------------------
# -----------------------------------------------------------------
index_tbl <- raw %>% mutate(
  NDVI   = (nir - red)                  / (nir + red),
  EVI    = 2.5 * (nir - red)            / (nir + 6 * red - 7.5 * blue + 1),
  GNDVI  = (nir - green)                / (nir + green),
  NDMI   = (nir - swir1)                / (nir + swir1),
  NDWI   = (green - nir)                / (green + nir),
  NBR    = (nir - swir2)                / (nir + swir2),
  NBR2   = (swir1 - swir2)              / (swir1 + swir2),
  ARVI   = (nir - (2 * red - blue))     / (nir + (2 * red - blue)),
  OSAVI  = (nir - red)                  / (nir + red + 0.16),
  NIRv   = ((nir - red) / (nir + red))  *  nir,
  VMI    = (nir - swir1) * (nir - swir2) / (nir + swir1 + swir2),
  NMDI   = (nir - (swir1 - swir2))      / (nir + (swir1 - swir2)),
  MSI    = swir1 / nir,
  MSAVI2 = (nir * swir1 - swir2 * red)  / (nir * swir1 + swir2 * red + 0.001),
  VMSI   = ((nir + swir1) - (swir2 + red)) /
    ((nir + swir1) + (swir2 + red) + 0.001),
  GSDI   = (green - swir2)              / (green + swir2),
  CNSI   = (nir * swir1)                / (swir2 * red + 0.001),
  EVDI   = (nir - red - swir2)          / (nir + red + swir2),
  TBMI   = (swir1 - swir2 - red)        / (swir1 + swir2 + red),
  SWIR12 = nir / swir2
)%>%                                         # <- all indices now created
  select(-nir, -red, -green, -blue, -swir1, -swir2)   #  <-- drop raw bands

# -----------------------------------------------------------------
# 3)  Long format  →  histogram facet grid -------------------------
# -----------------------------------------------------------------
index_long <- index_tbl %>%
  pivot_longer(everything(), names_to = "index", values_to = "value")

ggplot(index_long, aes(x = value)) +
  geom_histogram(bins = 200, fill = "steelblue", alpha = 0.85) +
  facet_wrap(~ index, scales = "free", ncol = 4) +
  labs(
    title = "Theoretical distribution of 19 remote-sensing indices\n(0–1 random reflectance inputs)",
    x     = "Index value",
    y     = "Count (random samples)"
  ) +
  theme_minimal(base_size = 11)
