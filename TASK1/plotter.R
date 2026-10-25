library(ggplot2)

files = list(
  "output/SOBI_S_10_40_500.csv",
  "output/SOBI_F_10_40_500.csv",
  "output/AMUSE_S_10_40_500.csv",
  "output/AMUSE_F_10_40_500.csv",
  "output/SOBI_S_10_40_1000.csv",
  "output/SOBI_F_10_40_1000.csv",
  "output/AMUSE_S_10_40_1000.csv",
  "output/AMUSE_F_10_40_1000.csv",
  "output/SOBI_S_10_40_2000.csv",
  "output/SOBI_F_10_40_2000.csv",
  "output/AMUSE_S_10_40_2000.csv",
  "output/AMUSE_F_10_40_2000.csv"
)

names = list(
  "SOBI_S_10_40_500",
  "SOBI_F_10_40_500",
  "AMUSE_S_10_40_500",
  "AMUSE_F_10_40_500",
  "SOBI_S_10_40_1000",
  "SOBI_F_10_40_1000",
  "AMUSE_S_10_40_1000",
  "AMUSE_F_10_40_1000",
  "SOBI_S_10_40_2000",
  "SOBI_F_10_40_2000",
  "AMUSE_S_10_40_2000",
  "AMUSE_F_10_40_2000"
)

results = lapply(1:12, function(i) {
  
  # Extract parameters from filename
  name = names[[i]]
  parts = strsplit(name, "_")[[1]]

  method = parts[1]
  bootstrap = parts[2]
  lag = as.integer(parts[3])
  T = as.integer(parts[5])
  
  # Read the three columns (L = 10, 20, 40)
  data = read.csv(files[[i]])
  
  data.frame(
    q_hat = unlist(data-1, use.names = FALSE),
    method = paste(method, bootstrap, sep = "_"),
    T = T,
    L = rep(c(10, 20, 40), each = nrow(data))
  )
})

results = do.call(rbind, results)
results$method_label = paste(results$method, results$L, sep="_")

# Method levels the same as Nordhausen 2018 to replicate results.
method_levels = c(
  "AMUSE_F_10", "AMUSE_F_20", "AMUSE_F_40",
  "AMUSE_S_10", "AMUSE_S_20", "AMUSE_S_40",
  "SOBI_F_10", "SOBI_F_20", "SOBI_F_40",
  "SOBI_S_10", "SOBI_S_20", "SOBI_S_40"
)

# Match the level to the method label in the data frame.
results$method_label = factor(
  results$method_label,
  levels = method_levels
)

p = ggplot(results, aes(x = method_label, fill = factor(q_hat, levels = 0:9))) +
  geom_bar() +
  coord_flip() +
  facet_grid(. ~ T) +
  labs(x = "Method", y = "Count", fill = expression(q)) +
  theme_bw() +
  theme(legend.position = "right", panel.grid.minor = element_blank())


ggsave("q_bar_chart.pdf", plot = p, width = 11.38, height = 4)
