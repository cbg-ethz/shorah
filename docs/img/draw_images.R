library(tidyverse)
library(reshape2)
library(xkcd)
set.seed(246)
### Already have read in fonts (see previous answer on how to do this)
loadfonts()

plot_type <- 'shotgun'
#plot_type <- 'amplicon'

### Set up the data
line_len <- 30
gen_start <- 1
gen_stop <- 20

if(plot_type == 'amplicon'){
  reg_start <- 7
  reg_stop <- 13
  read_len <- 6
  read_start <- rnorm(n=10, mean=reg_start, sd=0.07)
  read_stop <- rnorm(n=10, mean=reg_stop, sd=0.07)
}


if(plot_type == 'shotgun'){
  reg_start <- 3
  reg_stop <- 18
  read_start <- runif(n=10, min=reg_start, max=reg_stop-read_len)
  read_stop <- read_start + rnorm(n=10, mean=read_len, sd=0.1)
}

### Everything in a list that will be melted
datalist <- list()

### start with the genome and the region
datalist[[1]] <- data.frame(x=seq(gen_start, gen_stop, length.out=line_len),
                            y=rep(13, line_len),
                            obj = rep('genome', line_len))

datalist[[2]] <- data.frame(x=seq(reg_start, reg_stop, length.out=line_len),
                            y=rep(12, line_len),
                            obj = rep('region', line_len))


### then add the reads
for (i in 1:10){
  datalist[[i + 2]] <- data.frame(x=seq(read_start[i], read_stop[i],
                                    length.out=line_len),
                              y=rep(i, line_len),
                              obj=rep('read', line_len))
}

melted <- melt(datalist, id.vars=c("x", "obj"))

cbPalette <- c("#E69F00", "#F35B04", "#00A6A6", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
### Plot the chart
p <- ggplot(melted, aes(x, value, group=value, colour=obj)) +
  geom_smooth(position=position_jitter(h = 0.02), size=1.5) +
  labs(x="", y="") + 
  scale_colour_manual(name="", values=cbPalette) +
  theme_xkcd() + theme(axis.ticks = element_blank(), axis.text = element_blank()) + 
  ggtitle(plot_type)
p

#ggsave("xkcd_ggplot.pdf", plot=p, width=8, height=5)
#embed_fonts("xkcd_ggplot.pdf", outfile="xkcd_ggplot_embed.pdf")
