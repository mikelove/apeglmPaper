library(iCOBRA)

# example: load the Bottomly simulation data, replicate 1
cd <- readRDS("cobradata_b1.rds")

# you can use the iCOBRA app to explore the NBR vs FSR curves
# (there are no null genes in this simulation, so no meaningful FDR)

# 1) launch the iCOBRA app
# 2) Click "Truth"
# 3) Select column containing continuous truth: 'lfc'
#    optionally, Select variable to stratify by: 'mean.strat'
# 4) Click "Results"
# 5) Click "Start calculation!"
# 6) Click "NBR vs FSR"
# 7) Click "Plot settings" to customize, e.g. add curve, change threshold
COBRAapp(cd)

# to calculate iCOBRA plots in R:
cp <- calculate_performance(cd, cont_truth = "lfc", svalthrs=c(.01,.05),
                            aspects=c("scatter","fsrnbr","fsrnbrcurve"))

# use the same color palette as the paper plots
cplot <- prepare_data_for_plot(cp, colorscheme="hue_pal")

# number of genes over FSR (apeglm, ashr, DESeq2)
plot_fsrnbrcurve(cplot, xaxisrange=c(0,.1))

# the table giving nominal and observed FSR
fsrnbr(cp)

# scatter plot of LFCs
plot_scatter(cplot, pointsize=0.5)
