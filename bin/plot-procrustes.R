#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)

# parse arguments -----------------------------------------------------------------------
parser = ArgumentParser(description='Plot procrustes analysis outcome')
parser$add_argument('-i', '--input', type='character', help="Input data frame with columns: j,id,tot.var,procrust.cor,var.log")
parser$add_argument('-o1', '--output1', type='character', help="Output figure")
parser$add_argument('-o2', '--output2', type='character', help="Output figure 2")
parser$add_argument('--label', action='store_true')
parser$add_argument('--threshold', type='double', default=1, help="Procrustes cor >= threshold to be labelled")
parser = parser$parse_args()

# read input
input = fread(parser$input, header=T)

# plot figures --------------------------------------------------------------------------
png(filename=parser$output1, width=1000, height=380, res=100)
par(mfrow=c(1,3))
pos = which( input$procrust.cor > as.numeric(quantile(input$procrust.cor, parser$threshold)) )
col = ifelse(input$procrust.cor > as.numeric(quantile(input$procrust.cor, parser$threshold)),"red", "black")
labels = rep('', nrow(input))
labels[pos] = input$id[pos]
plot(
    x=input$j, 
    y=input$procrust.cor, 
    col=col, 
    xlab="Ref gene index", 
    ylab="Procrustes correlation"
    )
if (parser$label) {
    text(
        x=input$j,
        y=input$procrust.cor,
        labels=labels
    )
}
plot(
    x=input$j, 
    y=input$var.log, 
    col=col,
    xlab="Ref gene index", 
    ylab="ALR ref variance"
    )
if (parser$label) {
    text(
        x=input$j,
        y=input$var.log,
        labels=labels
    )
}
plot(
    x=input$var.log, 
    y=input$procrust.cor, 
    col=col,
    xlab="ALR ref variance", 
    ylab="Procrustes correlation")
if (parser$label) {
    text(
        x=input$var.log,
        y=input$procrust.cor,
        labels=labels
    )
}

data = data.frame(x=input$id, y=input$procrust.cor)
data = data[order(data$y),]
data$x = factor(data$x, levels=data$x)
g = ggplot(data, aes(x=x, y=y)) +
      geom_point() + 
    #   geom_segment( aes(x=x, xend=x, y=0, yend=y)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=5)) +
      xlab("Ref gene") + 
      ylab("Procrustes correlation") + 
      ylim(0.8,1)
ggsave(parser$output2, plot=g, width=18, height=3)

# png(filename=parser$output2, width=1000, height=380, res=100)
# barplot(
#     input$procrust.cor,
#     xlab="Ref gene",
#     ylab="Procrustes correlation",
#     names.arg = input$id
#     )
