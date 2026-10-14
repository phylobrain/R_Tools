
# Species proportion barplot
percent <- integration_object@meta.data %>%
  count(Cell_Identity, Species) %>%
  group_by(Cell_Identity) %>%
  mutate(Percent = n / sum(n) * 100) %>%
  ungroup()

ggplot(percent, aes(
  x = Cell_Identity,
  y = Percent,
  fill = factor(Species, levels = c("Species1", "Species2"))
)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = paste0(round(Percent), "%")),
    position = position_stack(vjust = 0.5)
  ) +
  coord_flip() +
  labs(x = "Cluster", y = "Percentage (%)", fill = "Species") +
  scale_fill_manual(values = c("steelblue3", "firebrick2")) +
  theme_classic()