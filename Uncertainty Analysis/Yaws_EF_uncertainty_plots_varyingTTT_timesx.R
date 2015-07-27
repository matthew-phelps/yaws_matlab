## intro
rm(list = ls())
mac <- "/Users/Matthew/Google Drive/Yaws project/MATLAB/Uncertainty Analysis"
pc <- 'C:/Users/wrz741/Google Drev/Yaws project/MATLAB/Uncertainty Analysis'

setwd(mac)
library (SDMTools) # used for calculating weighted standard deviations
library (MASS) # used for fiting distributions
library (vcd) # used for goodness-of-fit tests
library (ggplot2)
library (stats)
library (reshape) # for renaming variables\
library (R.matlab) # for MATLAB files






# 1 treatment round -------------------------------------------------------


# load point estimates
Point.estimates.1x <-readMat("pe_1x.mat")$PE.1x

# load results of uncertainty analysis:
results.1x <- readMat("results1x_n.mat")$results
ef.1x <- readMat("results1x_Ef.mat")$results

# find min and max for each coverage level
infect.range = matrix(0, ncol(results.1x), 4)
for (i in 1:ncol(results.1x)){
    infect.range[i,1] = min(results.1x[,i])
    infect.range[i,2] = max(results.1x[,i])
    infect.range[i,3] = min(ef.1x[,i])
    infect.range[i,4] = max(ef.1x[,i])
}
infect.range = as.data.frame(infect.range)


Point.estimates.1x = as.data.frame(Point.estimates.1x)

x.axis = c("50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90", "95%")

plot.title = c("1 treatment round")

plot.1x <- ggplot(infect.range, aes( x = as.numeric(rownames(infect.range)))) +
    geom_ribbon( aes( ymin = V1, ymax = V2), alpha = 0.4 ) +
    geom_ribbon( aes( ymin = V3, ymax = V4), alpha = 0.4 ) +
    geom_line( data = Point.estimates.1x, aes( y = V1), size = 0.8) +
    xlab("Coverage level") +
    ylab("Percentage infected") +
    ggtitle (bquote(atop(.(plot.title[1])))) +
    geom_hline( aes( yintercept = 0.0016), linetype="longdash", alpha = 0.5, color ="darkred") +
    scale_x_discrete( labels=c("50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90%", "95%")) 



plot.1x








# 2 treatment rounds 3 months ------------------------------------------------------


# load point estimates
Point.estimates.2x <-readMat("pe_2x.mat")$PE.2x

# load results of uncertainty analysis:
results.2x <- readMat("results2x_n_3mo.mat")$results
ef.2x <- readMat("results2x_Ef.mat")$results

# find min and max for each coverage level
infect.range = matrix(0, ncol(results.1x), 4)
for (i in 1:ncol(results.2x)){
    infect.range[i,1] = min(results.2x[,i])
    infect.range[i,2] = max(results.2x[,i])
    infect.range[i,3] = min(ef.2x[,i])
    infect.range[i,4] = max(ef.2x[,i])
}
infect.range = as.data.frame(infect.range)


Point.estimates.2x = as.data.frame(Point.estimates.2x)

x.axis = c("50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90", "95%")

plot.title = c("2 treatment round")

plot.2x <- ggplot(infect.range, aes( x = as.numeric(rownames(infect.range)))) +
    geom_ribbon( aes( ymin = V1, ymax = V2), alpha = 0.4 ) +
    geom_ribbon( aes( ymin = V3, ymax = V4), alpha = 0.4 ) +
    geom_line( data = Point.estimates.2x, aes( y = V1), size = 0.8) +
    xlab("Coverage level") +
    ylab("Percentage infected") +
    ggtitle (bquote(atop(.(plot.title[1])))) +
    geom_hline( aes( yintercept = 0.0016), linetype="longdash", alpha = 0.5, color ="darkred") +
    scale_x_discrete( labels=c("50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90%", "95%"))+
    scale_y_continuous(limits=c(0, .04))


plot.2x




# 2 treatment rounds 3 months ------------------------------------------------------


# load point estimates
Point.estimates.2x <-readMat("pe_2x.mat")$PE.2x

# load results of uncertainty analysis:
results.2x <- readMat("results2x_n_9mo.mat")$results
ef.2x <- readMat("results2x_Ef.mat")$results

# find min and max for each coverage level
infect.range = matrix(0, ncol(results.1x), 4)
for (i in 1:ncol(results.2x)){
  infect.range[i,1] = min(results.2x[,i])
  infect.range[i,2] = max(results.2x[,i])
  infect.range[i,3] = min(ef.2x[,i])
  infect.range[i,4] = max(ef.2x[,i])
}
infect.range = as.data.frame(infect.range)


Point.estimates.2x = as.data.frame(Point.estimates.2x)

x.axis = c("50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90", "95%")

plot.title = c("2 treatment round")

plot.2x <- ggplot(infect.range, aes( x = as.numeric(rownames(infect.range)))) +
  geom_ribbon( aes( ymin = V1, ymax = V2), alpha = 0.4 ) +
  geom_ribbon( aes( ymin = V3, ymax = V4), alpha = 0.4 ) +
  geom_line( data = Point.estimates.2x, aes( y = V1), size = 0.8) +
  xlab("Coverage level") +
  ylab("Percentage infected") +
  ggtitle (bquote(atop(.(plot.title[1])))) +
  geom_hline( aes( yintercept = 0.0016), linetype="longdash", alpha = 0.5, color ="darkred") +
  scale_x_discrete( labels=c("50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90%", "95%"))+
  scale_y_continuous(limits=c(0, .04))


plot.2x



# 3 treatment rounds ------------------------------------------------------


# load point estimates
Point.estimates.3x <-readMat("pe_3x.mat")$PE.3x

# load results of uncertainty analysis:
results.3x <- readMat("results3x_n_3mo.mat")$results
ef.3x <- readMat("results3x_Ef.mat")$results

# find min and max for each coverage level
infect.range = matrix(0, ncol(results.1x), 4)
for (i in 1:ncol(results.2x)){
  infect.range[i,1] = min(results.3x[,i])
  infect.range[i,2] = max(results.3x[,i])
  infect.range[i,3] = min(ef.3x[,i])
  infect.range[i,4] = max(ef.3x[,i])
}

infect.range = as.data.frame(infect.range)
Point.estimates.3x = as.data.frame(Point.estimates.3x)

x.axis = c("50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90", "95%")

plot.title = c("3 treatment rounds")

plot.3x <- ggplot(infect.range, aes( x = as.numeric(rownames(infect.range)))) +
    geom_ribbon( aes( ymin = V1, ymax = V2), alpha = 0.4) +
    geom_ribbon( aes( ymin = V3, ymax = V4), alpha = 0.4) +
    geom_line( data = Point.estimates.3x, aes( y = V1), size = 0.8) +
    xlab("Coverage level") +
    ylab("Percentage infected") +
    ggtitle (bquote(atop(.(plot.title[1])))) +
    geom_hline( aes( yintercept = 0.0016), linetype="longdash", alpha = 0.5, color ="darkred") +
    scale_x_discrete( labels=c("50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90%", "95%")) +
    scale_y_continuous(limits=c(0, .04))


plot.3x


# # 3 treatment rounds ------------------------------------------------------
# 
# 
# # load point estimates
# Point.estimates.3x <-readMat("pe_3x.mat")$PE.3x
# 
# # load results of uncertainty analysis:
# results.3x <- readMat("results3x_n_9mo.mat")$results
# ef.3x <- readMat("results3x_Ef.mat")$results
# 
# # find min and max for each coverage level
# infect.range = matrix(0, ncol(results.1x), 4)
# for (i in 1:ncol(results.2x)){
#   infect.range[i,1] = min(results.3x[,i])
#   infect.range[i,2] = max(results.3x[,i])
#   infect.range[i,3] = min(ef.3x[,i])
#   infect.range[i,4] = max(ef.3x[,i])
# }
# 
# infect.range = as.data.frame(infect.range)
# # Point.estimates.3x = as.data.frame(Point.estimates.3x)
# 
# x.axis = c("50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90", "95%")
# 
# plot.title = c("3 treatment rounds")
# 
# plot.3x <- ggplot(infect.range, aes( x = as.numeric(rownames(infect.range)))) +
#   geom_ribbon( aes( ymin = V1, ymax = V2), alpha = 0.4) +
#   geom_ribbon( aes( ymin = V3, ymax = V4), alpha = 0.4) +
#   geom_line( data = Point.estimates.3x, aes( y = V1), size = 0.8) +
#   xlab("Coverage level") +
#   ylab("Percentage infected") +
#   ggtitle (bquote(atop(.(plot.title[1])))) +
#   geom_hline( aes( yintercept = 0.0016), linetype="longdash", alpha = 0.5, color ="darkred") +
#   scale_x_discrete( labels=c("50%", "55%", "60%", "65%", "70%", "75%", "80%", "85%", "90%", "95%")) +
#   scale_y_continuous(limits=c(0, .04))
# 
# 
# plot.3x






# Final multiplot ---------------------------------------------------------
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

multiplot (plot.1x, plot.2x, plot.3x, cols=1)
