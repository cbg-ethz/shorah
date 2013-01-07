library(ggplot2)
library(reshape)
library(extrafont)
set.seed(2)
### Already have read in fonts (see previous answer on how to do this)
loadfonts()

### XKCD theme
theme_xkcd <- theme(
  panel.background = element_rect(fill=NA), 
  axis.ticks = theme_blank(), # element_line(colour="white"),
  panel.grid = theme_blank(), # element_line(colour=NA),
  axis.text.y = element_text(colour=NA), 
  axis.text.x = element_text(colour=NA),
  text = element_text(size=16, family="Humor Sans")
)

### Set up the data
line_len <- 30
gen_start <- 1
gen_stop <- 20
reg_start <- 7
reg_stop <- 13
read_len <- 6
read_start <- rnorm(n=10, mean=reg_start, sd=0.07)
read_stop <- rnorm(n=10, mean=reg_stop, sd=0.07)

### uncomment these lines for shotgun sequencing
#read_start <- runif(n=10, min=reg_start, max=reg_stop-read_len)
#read_stop <- read_start + rnorm(n=10, mean=read_len, sd=0.1)

### Everything in a list that will be melted
datalist <- list()

### start with the genome and the region
datalist[[1]] <- data.frame(x=seq(gen_start, gen_stop, length.out=line_len),
                             y=rep(13, line_len), colour=rep("blue", line_len))

datalist[[2]] <- data.frame(x=seq(reg_start, reg_stop, length.out=line_len),
                            y=rep(12, line_len), colour=rep("green", line_len))


### then add the reads
for (i in 1:10){
  datalist[[i + 2]] <- data.frame(x=seq(read_start[i], read_stop[i],
                                    length.out=line_len),
                              y=rep(i, line_len), colour=rep("orange", line_len))
}

melted <- melt(datalist, id.vars=c("x", "colour"))
### Plot the chart
p <- ggplot(melted, aes(x, value, group=value, colour=colour)) +
 geom_line(position=position_jitter(h = 0.02), size=1.5) + theme_xkcd +
 labs(x="", y="") + 
 scale_colour_discrete(name="", breaks=c("blue", "green", "orange"),
                       labels=c("genome", "region", "reads")) +
 theme(legend.text = element_text(colour="black", size = 16, face = "bold"))

p

ggsave("xkcd_ggplot.pdf", plot=p, width=8, height=5)
embed_fonts("xkcd_ggplot.pdf", outfile="xkcd_ggplot_embed.pdf")