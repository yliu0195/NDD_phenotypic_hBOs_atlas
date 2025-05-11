##beginning Lu's data

setwd("E:/Users/thinkpad/工作/简历/Resume/找 Gleeson/Projects/Research In Progress/2023-11-14-Lu-PCA-plot")
raw<-read.csv(file="Lu-PCA-plot.csv",header=T)
library(cowplot)
library(dplyr)
library(ggplot2)

x <- raw[6:8]
pc <- prcomp(x)
df <- cbind(pc$x[,1:2], raw[,1]) %>% as.data.frame()
df$PC1 <- as.numeric(df$PC1) / (pc$sdev[1] * sqrt(nrow(raw)))
df$PC2 <- as.numeric(df$PC2) / (pc$sdev[2] * sqrt(nrow(raw)))
df$V3 <- as.factor(df$V3)
# plot_density
p1 <- ggplot(df, aes(PC1, PC2, colour = V3)) +
  geom_point(size = 3, aes(shape = V3)) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
               data = df[df$V3 == "1" | df$V3 == "2",], size = 1)

# Add density curves to y and x axis
xdens <- 
  axis_canvas(p1, axis = "x") + 
  geom_density(data = df, aes(x = PC1, fill = V3, colour = V3), alpha = 0.3)
ydens <-
  axis_canvas(p1, axis = "y", coord_flip = TRUE) + 
  geom_density(data = df, aes(x = PC2, fill = V3, colour = V3), alpha = 0.3) +
  coord_flip()
p1 %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()


# plot_box
p1 <- ggplot(df, aes(PC1, PC2, colour = V3)) +
  geom_point(size = 3, aes(shape = V3)) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
               data = df[df$V3 == "1" | df$V3 == "2",], size = 1)

# Add density curves to y and x axis
xdens <- 
  axis_canvas(p1, axis = "x") + 
  geom_boxplot(data = df, aes(x = PC1, fill = V3, colour = V3), alpha = 0.3,notch=TRUE)
ydens <-
  axis_canvas(p1, axis = "y", coord_flip = TRUE) + 
  geom_boxplot(data = df, aes(x = PC2, fill = V3, colour = V3), alpha = 0.3,notch=TRUE) +
  coord_flip()
p1 %>%
  insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  ggdraw()
    