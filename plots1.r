# 1. Load required package
library(tidyverse)

# 2.
setwd("C:/Users/nisha/Downloads")

# 3. Load
url <- "https://raw.githubusercontent.com/octantbio/mc4r-dms/main/paper/mc4r-dms.tsv"
df <- read_tsv(url)

# 4. Create and format 'compound_dose'
df <- df %>%
  mutate(compound_dose = paste(compound, dose)) %>%
  mutate(compound_dose = trimws(compound_dose)) %>%
  # Convert it to a factor to enforce the row order
  mutate(compound_dose = factor(compound_dose, levels = c("aMSH high", "aMSH low", "Setmelanotide")))

# 5. Create the Heatmap 
heatmap_plot <- ggplot(df, aes(x = pos, y = aa, fill = statistic)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_grid(compound_dose ~ pathway) +
  theme_minimal() +
  labs(title = "MC4R DMS Heatmap", x = "Position", y = "Amino Acid", fill = "Statistic")

print(heatmap_plot)

# 6. Create the Volcano Plot
volcano_plot <- ggplot(df, aes(x = log2FoldChange, y = -log10(p.adj))) +
  geom_point(alpha = 0.5, color = "darkgray") +
  facet_grid(compound_dose ~ pathway) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(p-adj)")

print(volcano_plot)

# 7. Prepare the PyMOL input data
pymol_input <- df %>%
  group_by(pos) %>%
  slice_max(order_by = statistic, n = 1, with_ties = FALSE) %>% # Gets the best score per position
  ungroup() %>%
  mutate(
    chain_id = "R", # Using 'R' for receptor chain
    res_id = pos,
    score = statistic
  ) %>%
  select(chain_id, res_id, score)

# 8. Export the PyMOL file
write_tsv(pymol_input, "best_mutations_for_pymol.txt", col_names = FALSE)
