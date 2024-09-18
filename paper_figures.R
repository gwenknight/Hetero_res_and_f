##### plots for paper 
library(tidyverse)
library(patchwork)

theme_set(theme_bw(base_size = 14))

## Example distribution
mres = 50
nfit = 50
x <- seq(1/mres,1,1/mres) # seq(from = 0, to = 1, length.out = mres)
y <- seq(1/nfit,1,1/nfit) #seq(from = 0, to = 1, length.out = nfit)
f <- function(x, y) dmvnorm(cbind(x, y), mean = c(0.5, 0.3),sigma = diag(c(1,1)/10))
z <- outer(x, y, FUN = f); z <- z/sum(z) # Generate and normalise
data <- expand.grid(x = x,y = y)
data$z = as.vector(z)
g1 <- ggplot(data, aes(x, y, fill= z)) + geom_tile() + 
  scale_fill_gradient(low="blue", high="red","Proportion of\nnew mutants") +
  scale_x_continuous("Relative resistance level",expand = c(0, 0)) +
  scale_y_continuous("Relative fitness level",expand = c(0, 0)) 

g2 <- ggplot(data %>% filter(y %in% seq(1/nfit, 1, 1/nfit)[seq(1,50,length = 9)]), aes(x=x, y = z)) + geom_bar(stat = "identity", aes(fill = z)) + 
  scale_fill_gradient(low="blue", high="red","Proportion of\nnew mutants") +
   facet_wrap(~y) + 
  scale_x_continuous("Relative resistance level",expand = c(0, 0)) +
  scale_y_continuous("Proportion of new mutants",expand = c(0, 0)) 

(g1 + theme(legend.position = "none")) + g2 + 
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')
ggsave("figures/initial_prop_eg.pdf", width = 12, height = 7)
