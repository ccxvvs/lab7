# Load required libraries
library(tidyverse)

# 1. 
df_exp <- read_tsv("mc4r-dms.tsv")
df_pred <- read_csv("predicted_ddg_results_part_1.csv")

# 2. Merge the data
df_merged <- inner_join(df_exp, df_pred, by = c("pos", "aa"))

# 3. 
df_merged <- df_merged %>%
  mutate(
    transformed_ddg = log10(pmax(predicted_ddg, 0) + 1),
    # Combine compound and dose for the facet grid rows
    compound_dose = paste(compound, dose, sep = "_"),
    # Explicitly calculate the negative log2 fold change to prevent formula errors
    neg_log2FoldChange = -log2FoldChange 
  )

# 4. 
correlations <- df_merged %>%
  group_by(pathway, compound_dose) %>%
  summarize(
    # Use the raw predicted_ddg for rank correlation, comparing against the negative L2FC
    rho = cor(predicted_ddg, neg_log2FoldChange, method = "spearman", use = "complete.obs"),
    p_val = cor.test(predicted_ddg, neg_log2FoldChange, method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  )

print("--- Spearman Rank Correlations ---")
print(correlations)

# 5. Plot the facet grid
ggplot(df_merged, aes(x = transformed_ddg, y = neg_log2FoldChange)) +
  geom_point(alpha = 0.2, size = 1) + 
  geom_smooth(method = "lm", color = "red") + 
  facet_grid(compound_dose ~ pathway) + 
  labs(
    title = "Experimental vs. Predicted ddG",
    x = "Rosetta predicted ddG (log10(x+1))",
    y = "-log2FoldChange (Howard2025)"
  ) +
  theme_minimal()
