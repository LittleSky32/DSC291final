setwd("D:/ucsd/dsc291")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read the results file
alz_results <- read.table("alz_baseline.results", header = TRUE, sep = "\t")

# Filter out rows with NaN or missing enrichment p-values
alz_results <- alz_results %>%
  filter(!is.na(Enrichment_p))


######### significant annotation ##########
# Define a significance threshold for enrichment p-values (e.g., 0.05)
# Filter significant results (p < 0.05) and exclude extreme high and low enrichment values
significant_results <- alz_results %>%
  filter(Enrichment_p < 0.05)
print(significant_results)
alz_filtered <- significant_results %>%
  filter(
    Enrichment < 100 & 
      Enrichment > -100  # Keep the reasonable values
  )

# Plot with enrichment and standard error bars
ggplot(alz_filtered, aes(x = reorder(Category, Enrichment), y = Enrichment)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error),
                width = 0.2, color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Enrichment Scores by Significant Annotation Category (Filtered)",
    x = "Category",
    y = "Enrichment"
  )


############## Enrichment on all anotations ##############
alz_filtered_all <- alz_results %>%
  filter(
    Enrichment < 100 & 
      Enrichment > -100  # Keep the reasonable values
  )
ggplot(alz_filtered_all, aes(x = reorder(Category, Enrichment), y = Enrichment)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error),
                width = 0.2, color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Enrichment Scores by Annotation Category",
       x = "Category",
       y = "Enrichment")
