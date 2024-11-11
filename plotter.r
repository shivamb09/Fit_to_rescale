library(ggplot2)
headers <- c("bp", "normalized")


file1 <- read.table("normalized_total.tsv", header = TRUE) 
file2 <- read.table("reference.hist", header = FALSE, col.names = headers)
file3 <- read.table("sub_sampled_1_plot.tsv", header = TRUE)


x_ticks <- seq(0, 700, by = 100) 
y_ticks <- seq(0, 0.025, by = 0.005) 


plot <- ggplot() +
 
  geom_line(data = file1, aes(x = bp, y = normalized, color = "Query"), size = 1, na.rm = TRUE) +
 
  geom_line(data = file2, aes(x = bp, y = normalized, color = "Reference"), size = 1, na.rm = TRUE) +
  
  geom_point(data = file3, aes(x = bp, y = normalized, color = "Rep 1"), size = 1.5, shape = 4, na.rm = TRUE) +
  
  
  scale_color_manual(values = c("Query" = "purple", "Reference" = "#1E90FF", "Rep 0" = "#00BFFF")) +
  

  scale_x_continuous(breaks = x_ticks, limits = c(0, 700)) +
  scale_y_continuous(breaks = y_ticks, limits = c(0, 0.025)) +
  
 
  labs(x = "bp", y = "normalized") +
  
 
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.5),
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),
    legend.title = element_blank(),
    legend.position = "top", 
    legend.text = element_text(size = 8) 
  )

ggsave("rescaled_plot.pdf", plot = plot, width = 8, height = 6)
