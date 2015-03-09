#  http://hnm.stat.cmu.edu/2014-07%20IMPS%20WORKSHOP/01%20Intro%20to%20Networks%20in%20R/intro-to-networks-demos.r
#
#
# intro-to-networks-demos.r
#
# Wed Mar 05 11:15:53 2014
#
# This file contains the R code to run (most of) the examples in the
# lecture notes "Introduction to Social Network Models".
#

###################################################################

# SECTION 1

# begin by loading in libraries, helper functions and data...

# 1. The igraph library, for exploratory data analysis (EDA) of
#    social network data...
#

install.packages('igraph', repos='http://lib.stat.cmu.edu/R/CRAN/')
library(igraph)

# 2. Some helper functions that make using igraph and the Pitts &
# Spillane data easier...
#
# This also creates the small network "small.example" which is used
# for many of the illustrations in the lecture notes...
#

source("intro-to-networks-helper-functions.r")

# 3. The Pitts & Spillane (2009) data...
#

load("schools-advice-data.Rdata")

#
# at this point you should see the following from ls():
#
# > ls()
# [1] "advice.mat"       "all.vars.mat"     "edge.centrality"  "edge.vars.mat"   
# [5] "extract.X"        "extract.Y"        "node.centrality"  "school.vars.mat" 
# [9] "small.example"    "stack.X"          "stack.Y"          "teacher.vars.mat"
#
# These objects are as follows:
#
# DATA SETS:
#
# advice.mat -- a 3 dimensional array containing 15 adjacency matrices
#               for the 15 teacher advice networks from Pitts &
#               Spillane (2009).
#
#               In each adjacency matrix, a 1 in position (i,j) means
#               that teacher i goes to teacher j for advice.
#
# edge.vars.mat -- a 3 dimensional array containing the edge/dyad
#               covariates for each of the 15 advice networks from
#               P&S.  The variables are:
#
#               same.yrs.in.schl 
#               same.innov.attitude 
#               teach.same.grade
#
# teacher.vars.mat -- a 3 dimensional array containing the teacher/node
#               covariates for each of the 15 advice networks from
#               P&S.  The variables are:
#
#               yrs.tchg.sender 
#               tchr.trust.sender 
#               yrs.tchg.recvr 
#               tchr.trus.recvr
#
# school.vars.mat --  a 3 dimensional array containing the school/network
#               covariates for each of the 15 advice networks from
#               P&S.  The variables are:
#
#               catholic 
#               school.size
#
# all.vars.mat -- a 3 dimensional array containing all of the
#               variables above.  Unfortunately it is not labelled so
#               it is difficult to tell which columns are which
#               variables.
#
#               Instead of using this, we will just assemble the
#               covariates we need with the helper functions.
#
# small.example -- a small advice network labelled like the
#               Pitts-Spillane networks, for making managable examples
#               for the slides.
#
#               Nb., small.example is an adjacency matrix, that will
#               have to be converted to an igraph network object for
#               most of our work.
#
# HELPER FUNCTIONS:
#
# node.centrality -- a helper function that prints out several
#               measures of node centrality for each of the nodes in
#               the graph G.
#
# edge.centrality -- a helper function that prints out edge-betweenness
#               for each edge in the graph G, and adds labels to make
#               it easier to see which betweenness number goes with
#               which edge.
#
# extract.X -- a helper function to extract the predictor variables
#              from any one network from the Pitts and Spillane data
#
# extract.Y -- a helper function to extract the response variable
#              (presence or absence of ties) from any one network from
#              the Pitts and Spillane data
#
# stack.X --  a helper function to extract the predictor variables
#              from any subset of networks from the Pitts and Spillane data
#
# stack.Y -- a helper function to extract the response variable
#              (presence or absence of ties) from any subset of
#              networks from the Pitts and Spillane data

########################################################################

# SECTION 2

# Let's convert small.example to a graph and graph it:
#

G <- graph.adjacency(small.example)

plot(G)

# if you are in Rstudio, you should type "x11()" first, so all the plots
# in this tutorial will go to a larger and more legible plotting window.)

# there are endless variations on this basic plot.
#
# a couple of good "stock" ones are
#
#       plot(G,layout=layout.fruchterman.reingold)
#
# and
#
#       plot(G,layout=layout.kamada.kawai)
#
# in this particular case because the graph is so simple, they produce
# almost the same result.  In general they do not.
#
# If you want to explore others, and try moving the points around
# yourself, try this:
#
# tkplot(G)
#
# try clicking and dragging nodes in the widnow that pops up when you
# do this.   ...and just close the window when you are done...


# small.example itself looks like this:

small.example

# (you can also view small.example dynamically in a spreadsheet-like
# window, by doing this:)

View(small.example)

# (just close the window when you are done)

# If we want to make a sociomatrix plot of it, we do this:

n <- dim(small.example)[1] # how many nodes
image(1:n,1:n,small.example,col=c("white","black"))
box()

# unfortunately that is not in the same orientation as the matrix in R
# is, so we will flip it around, and relabel it...

n <- dim(small.example)[1] # how many nodes
initials <- dimnames(small.example)[[1]]
image(1:n,1:n,t(small.example)[,n+1-(1:n)],col=c("white","black"),
      xlab="",ylab="",axes=F)
axis(1,at=1:n,labels=initials)
axis(2,at=1:n,labels=NA)
text(0,n+1-(1:n),labels=initials,crt=90)
box()

#####################################################################

# SECTION 3

# now let's calculate some measures of centrality -- which nodes (or
# edges) are in the "center" of the graph?
#
# the helper function node.centrality() collects together a few common
# node centrality statistics:

node.centrality(G)

# and the gelper function edge.centrality() prints a list of all the
# edges and their edge betweenness centralities.
#
# Because there are usually a lot of edges, just typing
#
#          edge.centrality(G)
#
# will just scroll through a lot of centrality statistics (and
# sometimes with no way to go back to the first few)...
#
# instead you can redirect the output to a text file, which you can
# look at with your favorite text editor later:

sink("edge-centralities.txt") # open the file
edge.centrality(G)            # write the edge centrality data onto it
sink()                        # close and save the file

# the following will work on my computer (because I have installed the
# text editor "emacs") but it probably won't work on yours:

system("emacs edge-centralities.txt")

# instead, just navigate to the folder where this text file was
# created, and double click on it to open it with whatever the default
# editor on your computer is...

###########################################################################

# SECTION 4

# OK, time for some communitiy detection!

# Here is how to make the two plots in the slides...

kk.layout <- layout.kamada.kawai(G) # make a consistent layout for the graphs

par(mfrow=c(1,2)) # make two plots side by side

# make the left plot
com <- edge.betweenness.community(G)
V(G)$color <- com$membership+1
plot(G, vertex.label.dist=1.5,layout=kk.layout,
     main="Edge-Betweenness Communities",
     edge.arrow.size=0.25)

# make the right plot
com <- walktrap.community(G)
V(G)$color <- com$membership+1
plot(G, vertex.label.dist=1.5,layout=kk.layout,
     main="Walktrap Communities",
     edge.arrow.size=0.25)

# ta da!!

# There are lots of different heuristic community-detection algorithms
# in the igraph package.  Here's some code that lets you compare
# several at once:

fn <- list(edge.betweenness=edge.betweenness.community,
           walktrap=walktrap.community,
           spinglass=spinglass.community,
           label.propagation=label.propagation.community,
           optimal=optimal.community)


g <- G
kk.layout <- layout.kamada.kawai(g)

par(mfrow=c(2,3))
for (i in 1:length(fn)) {
  com <- fn[[i]](g)
  V(g)$color <- com$membership+1
  oldmar <- par(mar=c(0,0,2,0))
  plot(g, vertex.label.dist=1.5,layout=kk.layout,main=names(fn)[i],
       edge.arrow.size=0.25)
  par(oldmar)
}

# by asking help for each blah.blah.community function, you can learn
# how each set of heuristics works...

########################################################################

# SECTION 5

# Let's make plots of all 15 teacher advice networks in the
# Pitts-Spillane data...

# first we'll make sociograms -- plots of the networks as graphs

# the basic idea is simple...

par(mfrow=c(3,5))
for (nwk in 1:15) {
  plot(graph.adjacency(advice.mat[[nwk]]))
}

# ...but ugly 'cause the networks are too squashed...
#
# ...so we give them some breathing room, and label the plots while
# we're at it...

par(mfrow=c(3,5))
for (nwk in 1:15) {
  oldmar <- par(mar=c(0,0,2,0))
  plot(graph.adjacency(advice.mat[[nwk]]),
       vertex.label=NA,
       edge.arrow.size=0.25,
       main=paste("School",nwk))
#  box()
  par(oldmar)
}

# we could also make sociomatrix plots... same idea...

# simple but ugly plots:

par(mfrow=c(3,5))
for (nwk in 1:15) {
    k <- dim(advice.mat[[nwk]])[1]
  image(1:k,1:k,advice.mat[[nwk]],col=c("white","black"))
}

# and plots that have been allowed to breath a little...

par(mfrow=c(3,5))
for (nwk in 1:15) {
  oldmar <- par(mar=c(2,2,2,1))
  k <- dim(advice.mat[[nwk]])[1]
  image(1:k,1:k,advice.mat[[nwk]],col=c("white","black"),
        main=paste("School",nwk))
  box()
  par(oldmar)
}

# (well, ok, still a bit busy, but now at least you can see the
# effects of different network sizes, different tie densities, etc.,
# across the 15 networks....)

############################################################################

# SECTION 6 --- EXERCISES!!

# Now you should explore some of the Pitts & Spillane data a little bit...
#

# try something like this:

example <- advice.mat[[11]]  # you don't have to choose school 11 -- choose
                             # any school from 1 to 15, or try it out on
                             # a couple of different schools...

example
example[is.na(example)] <- 0
example

G <- graph.adjacency(example)  # to make an igraph network object

# now you should:
#
# (1) make a plot of "example" as a sociomatrix (black and white squares)
#
# (2) make a plot of G as a network (sociogram)
#
# (3) compute the node.centrality() and edge.centrality() statistics
# and compare the results to your plot in (2) -- does it all make
# sense?  [Careful -- edge.centrality() makes a lot of output and you
# may wish to dump it to a text file with sink()!!]
#
# (4) try out some community detection.  I suggest trying the
# walktrap.community() function first, and then maybe try others.
# Don't forget to plot your results by imitating the code above!
#


############################################################################

# SECTION 7

# Fitting the Erdos-Renyi-Gilbert model to "small.example"...

y <- small.example   # adjacency matrix

diag(y) <- NA        # we don't want to include self-ties or self-non-ties
                     # in our analyses

y <- c(y)            # convert y from a matrix to a vector, for regression

e.r.g <- glm(y ~ 1, family=binomial) # logistic regression, intercept only...

summary(e.r.g)      # check out the fitted model... 

# Call:
# glm(formula = y ~ 1, family = binomial)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -0.9374  -0.9374  -0.9374   1.4381   1.4381  
# 
# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.5947     0.2202  -2.701  0.00692 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 117.15  on 89  degrees of freedom
# Residual deviance: 117.15  on 89  degrees of freedom
#   (10 observations deleted due to missingness)
# AIC: 119.15
# 
# Number of Fisher Scoring iterations: 4

exp(-0.5947) / (1 + exp(-0.5947))   # see what the estimated tie density is...
#[1] 0.3555572                      # and it matches the "direct" estimate!

#####################################################################

# SECTION 8

# Fitting the "sender-receiver" model to "small.example"...

# setting up y is the same as before:


y <- small.example   # adjacency matrix

diag(y) <- NA        # we don't want to include self-ties or self-non-ties
                     # in our analyses

y <- c(y)            # convert y from a matrix to a vector, for regression

# We need to write this model as y ~ X\beta, where X contains a dummy
# (0/1) variable for each teacher considered as a *sender*, and for
# each teacher considered as a *receiver*.  The code to create X
# follows.  Please feel free to explore it, ask questions about it,
# etc. after the workshop...

n <- dim(small.example)[1]
base <- matrix(0,nrow=n,ncol=n)
X <- NULL
for (i in 1:n) {              # dummies for each teacher as a "sender"
  tmp <- base
  tmp[i,] <- rep(1,n)
  X <- cbind(X,c(tmp))
}
for(i in 1:n) {               # dummies for each teacher as a "receiver"
  tmp <- base
  tmp[,i] <- rep(1,n)
  X <- cbind(X,c(tmp))
}
dimnames(X)[[2]] <- c(paste("a.",dimnames(small.example)[[2]],sep=""),
                      paste("b.",dimnames(small.example)[[2]],sep=""))

# you could View(X) now to see what you got, but it's a pretty big
# matrix and it might not be so obvious at first...)

# but now we can fit the model...

ab.model <- glm(y ~ X,family=binomial)

summary(ab.model)

# Call:
# glm(formula = y ~ X, family = binomial)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -1.7865  -0.8488  -0.6210   0.9162   2.4920  
# 
# Coefficients: (2 not defined because of singularities)
#               Estimate Std. Error z value Pr(>|z|)  
# (Intercept) -1.548e+00  1.170e+00  -1.323   0.1859  
# Xa.KPJ       6.592e-02  1.111e+00   0.059   0.9527  
# Xa.WQM      -9.162e-17  1.105e+00   0.000   1.0000  
# Xa.FOM       7.760e-01  1.070e+00   0.725   0.4685  
# Xa.SAE      -1.580e+00  1.410e+00  -1.120   0.2627  
# Xa.NYZ       5.705e-01  1.075e+00   0.531   0.5956  
# Xa.YAW       6.465e-01  1.079e+00   0.599   0.5492  
# Xa.EVN       1.103e+00  1.072e+00   1.028   0.3037  
# Xa.BWV       3.503e-01  1.077e+00   0.325   0.7450  
# Xa.WAP      -6.827e-01  1.182e+00  -0.578   0.5635  
# Xa.REK              NA         NA      NA       NA  
# Xb.KPJ       6.117e-01  1.117e+00   0.547   0.5840  
# Xb.WQM       1.904e-15  1.176e+00   0.000   1.0000  
# Xb.FOM       1.721e+00  1.108e+00   1.554   0.1202  
# Xb.SAE       9.936e-01  1.083e+00   0.917   0.3589  
# Xb.NYZ       6.855e-02  1.179e+00   0.058   0.9536  
# Xb.YAW       6.852e-01  1.120e+00   0.612   0.5405  
# Xb.EVN       1.513e-01  1.176e+00   0.129   0.8977  
# Xb.BWV       2.851e+00  1.221e+00   2.335   0.0195 *
# Xb.WAP      -5.232e-02  1.170e+00  -0.045   0.9643  
# Xb.REK              NA         NA      NA       NA  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 117.147  on 89  degrees of freedom
# Residual deviance:  97.687  on 71  degrees of freedom
#   (10 observations deleted due to missingness)
# AIC: 135.69
# 
# Number of Fisher Scoring iterations: 5

####################################################################

# SECTION 9

# Let's try some dyadic independence models to see if any of the
# covariates matter for predicting ties in the Pitts & Spillane
# networks...

# The helper functions "extract.Y()" and "extract.X()" will set up a Y
# vector and an X matrix for doing regression on any one of the 15
# Pitts-Spillane networks.  You just have to specify a number from 1
# to 15 in the parentheses.  It's better than constructing these
# things by hand as we did above!

# Here's the example from the slides:

Y <- extract.Y(1)  # extract Y and X for network #1
X <- extract.X(1)

dimnames(X)[[2]]   # remind me what the predictor variables are... ?
# [1] "same.yrs.in.schl"    "same.innov.attitude" "teach.same.grade"   
# [4] "yrs.tchg.sender"     "tchr.trust.sender"   "yrs.tchg.recvr"     
# [7] "tchr.trus.recvr"     "catholic"            "school.size"

test01.glm <- glm(Y ~ teach.same.grade,
                  family=binomial,data=as.data.frame(X))

summary(test01.glm)

# Call:
# glm(formula = Y ~ teach.same.grade, family = binomial, data = as.data.frame(X))
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -0.5792  -0.3348  -0.3348  -0.3348   2.4123  
# 
# Coefficients:
#                  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       -2.8536     0.1626 -17.549  < 2e-16 ***
# teach.same.grade   1.1532     0.2877   4.009 6.11e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 441.85  on 869  degrees of freedom
# Residual deviance: 427.59  on 868  degrees of freedom
#   (30 observations deleted due to missingness)
# AIC: 431.59
# 
# Number of Fisher Scoring iterations: 5


Y <- extract.Y(11)  # extract Y and X for network #11
X <- extract.X(11)

dimnames(X)[[2]]   # remind me what the predictor variables are... ?

test11.glm <- glm(Y ~ teach.same.grade,
                  family=binomial,data=as.data.frame(X))

summary(test11.glm)

# Call:
# glm(formula = Y ~ teach.same.grade, family = binomial, data = as.data.frame(X))
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -0.5078  -0.5078  -0.5078  -0.5078   2.0553  
# 
# Coefficients:
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -1.9833     0.2754  -7.202 5.95e-13 ***
# teach.same.grade  -15.5828  1398.7210  -0.011    0.991    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 93.470  on 131  degrees of freedom
# Residual deviance: 91.474  on 130  degrees of freedom
#   (12 observations deleted due to missingness)
# AIC: 95.474
# 
# Number of Fisher Scoring iterations: 16


# And if you wanted to see how all the variables work together, you
# could try something like this:

Y <- extract.Y(3)  # extract Y and X for network #3
X <- extract.X(3)

test03.glm <- glm(Y ~ X, family=binomial)

summary(test03.glm)

# Call:
# glm(formula = Y ~ X, family = binomial)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -1.0319  -0.5713  -0.4224  -0.3009   2.6887  
# 
# Coefficients: (2 not defined because of singularities)
#                      Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -1.89344    1.38890  -1.363  0.17280    
# Xsame.yrs.in.schl    -0.37940    0.26486  -1.432  0.15201    
# Xsame.innov.attitude  0.45379    0.26017   1.744  0.08113 .  
# Xteach.same.grade     0.59449    0.27777   2.140  0.03234 *  
# Xyrs.tchg.sender     -0.12730    0.04361  -2.919  0.00351 ** 
# Xtchr.trust.sender    0.01710    0.31501   0.054  0.95672    
# Xyrs.tchg.recvr      -0.17853    0.04562  -3.913  9.1e-05 ***
# Xtchr.trus.recvr      0.32637    0.31278   1.043  0.29673    
# Xcatholic                  NA         NA      NA       NA    
# Xschool.size               NA         NA      NA       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 452.12  on 599  degrees of freedom
# Residual deviance: 414.72  on 592  degrees of freedom
#   (25 observations deleted due to missingness)
# AIC: 430.72
# 
# Number of Fisher Scoring iterations: 5

######################################################################

# SECTION 10

# Finally, although it is not illustrated in the slides, it is also
# possible to stack the data from several networks together to get a
# "joint analysis" of the predictors of tie formation, across several
# networks.
#
# This is a kind of "poor man's" version of a Hierarchical Network
# Model (HNM).  It is deficient from an HNM in two ways:
#
# (1) the coefficients are not allowed to vary at all from one network
#     to the next; this doesn't seem right in light of the fact that
#     e.g. networks 1 and 11 don't both show effects for teaching in
#     the same grade;
#
# (2) It does not contain any latent structure to account for
#     higher-order structure in the data, that is not simply emergent
#     from a dyadic independence model.
#
# But it is still interesting to fit, and maybe to compare to an HNM
# fit later...

# we use the "stack.Y" and "stack.X" helper functions to stack the
# data from several networks together for one big analysis.
#
# We will just illustrate it for analysis of all X variables on 15
# networks at once:

Y <- stack.Y(1:15)
X <- stack.X(1:15)

testAll.glm <- glm(Y ~ X, family=binomial)

summary(testAll.glm)

# Call:
# glm(formula = Y ~ X, family = binomial)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -1.4532  -0.4647  -0.3689  -0.2788   2.7243  
# 
# Coefficients:
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -0.647215   0.304436  -2.126  0.03351 *  
# Xsame.yrs.in.schl     0.290673   0.064318   4.519 6.20e-06 ***
# Xsame.innov.attitude  0.335707   0.057354   5.853 4.82e-09 ***
# Xteach.same.grade     1.041810   0.058991  17.660  < 2e-16 ***
# Xyrs.tchg.sender     -0.013314   0.004208  -3.164  0.00156 ** 
# Xtchr.trust.sender    0.060697   0.058571   1.036  0.30007    
# Xyrs.tchg.recvr      -0.012133   0.004192  -2.894  0.00380 ** 
# Xtchr.trus.recvr     -0.026933   0.058290  -0.462  0.64405    
# Xcatholic            -1.182524   0.111148 -10.639  < 2e-16 ***
# Xschool.size         -0.019818   0.001511 -13.119  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 10134  on 15985  degrees of freedom
# Residual deviance:  9330  on 15976  degrees of freedom
#   (433 observations deleted due to missingness)
# AIC: 9350
# 
# Number of Fisher Scoring iterations: 5

########################################################################

# SECTION 11

# Now it is your turn!
#
# Use the helper functions
#
# extract.Y()   # for single networks
# extract.X()
#
# and
#
# stack.Y()     # for groups of networks
# stack.X() 
#
# to explore regression analysis of the ties in the teacher networks.
#
# try it on different networks, or different groups of networks, and
# compare your results.
#

########################################################################

#
# This is the end of the demos!
#

########################################################################
