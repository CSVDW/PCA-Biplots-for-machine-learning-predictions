#
#
#                 Script for Paper
#    Biplots for understanding model predictions in digital soil mapping
#
#
#    Author: Stephan van der Westhuizen
#
#    Example data sets can be downloaded from github
#
#    Example data set: (1) LUCAS_aggregated, (2) biplot_data
#
#########################################################################################
#define dat - the aggregated LUCAS data (aggregated by median for each country)

dat = read.table("LUCAS_aggregated.txt", sep=",", header=T)
dat = dat[,-7] #exclude coarse fragments (CF)

#Define p, number of variables, and r, the lower dimension approximation
p = ncol(dat) 
r = 2

#Define J matrix
J = diag(c(1,1,rep(0,p-r)))

#standardise data to have zero-means and unit variances
X = scale(dat,scale=T,center=T)

#Perform principle component analysis (PCA), that used singular value decomposition (SVD)
#Set scale and center parameters to FALSE as data was already standardised 
my_pca = prcomp(X, scale. = F, center = F)

#call U, D and V
U = my_pca$x
D = my_pca$sdev^2
V = my_pca$rotation

#Calculate Xhat
Xhat <- X%*%V%*%J%*%t(V)

#Construct Gabriel biplot (arrows from the origin)
XV = X%*%V
XVJ = X%*%V%*%J #scores
VJ = V%*%J      #loadings

#variances explained by first two PCs
sigma = cov(XV)
sum(diag(sigma%*%J))/sum(diag(sigma))
#OR
pc1 = D[1]/sum(D)
pc2 = D[2]/sum(D)
pc1 + pc2

#predictivity of axes
axes_predictivity=diag(t(Xhat)%*%Xhat)/diag(t(X)%*%X)
(axes_predictivity_sorted = sort(axes_predictivity,decreasing=T))


#plot scores
plot(XVJ, xlab=paste0("PC1 (", round(pc1*100,1), "%)"), 
     ylab=paste0("PC2 (", round(pc2*100,1), "%)"), 
     main="Gabriel biplot")
#add country names
text(XVJ, row.names(XVJ))

# scaling factor for axes
f = .55
#add arrows (axes / scores)
arrows(0,0,VJ[,1]*D[1]/f,
       VJ[,2]*D[2]/f,
       length = 0.1,
       lwd=  2,
       angle = 20,
       col = "darkgreen")

# Plot annotations
text(VJ[,1]*D[1]/(f-.04),
     VJ[,2]*D[2]/(f-.04),
     rownames(V),
     col = "darkgreen",
     cex = 1.2)


###############################################
### Gower biplot (predictive axes)

#based on 
#https://corybrunson.github.io/ordr/articles/ordr.html
install.packages("remotes")
library(remotes)
remotes::install_github("corybrunson/ordr")
library(ordr)
library(viridis)
library(ggplot2)

d = read.csv("biplot_data.csv") #same data used for Shiny app
Xnames = colnames(d)[-1] #names of covariates
p = length(Xnames)
prc_dat = d[,c(Xnames,"SOC")] #rearrange data so that the first p columns are the covariates
breaks = quantile(prc_dat[,"SOC"],probs=seq(0,1,length.out=7)) #could categorise soc predictions with quantile function, used 7 breaks here
labels = paste0(round(breaks,1))[-1]
prc_dat$response_cat <- cut(prc_dat[,"SOC"], breaks = breaks, labels = labels, include.lowest=T)
prc = ordinate(prc_dat,model=~prcomp(.,center=T,scale.=T),cols=1:p) #ignore potential warning notification

#predictive axes with points
ggbiplot(prc, axis.type = "predictive", axis.percents = T) +
  theme_biplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank()) +
  geom_rows_point() +
  geom_rows_point(aes(col = response_cat),alpha=.9) + scale_color_brewer(palette = "RdYlBu") +
  geom_cols_axis(aes(label = name, center = center,scale=scale, label_dodge = 0.001)) + 
  guides(col=guide_legend(title="SOC", reverse=T))


#   For alpha bags the user can refer to 
#
#   https://stackoverflow.com/questions/31893559/r-adding-alpha-bags-to-a-2d-or-3d-scatterplot
#
#   these functions are based on functions in the aplpack package (see below for functions and more examples)


StatBag <- ggplot2::ggproto("Statbag",
                            ggplot2::Stat,compute_group = 
                              function(data, scales, prop = 0.5) {
                     plothulls_ <-function(x, y, fraction, n.hull = 1, col.hull, lty.hull,
                                lwd.hull, ensity=0, ...){
                         if(ncol(x) == 2){ y <- x[, 2]; x <- x[, 1] }
                         n <- length(x)
                         if(!missing(fraction)) {
                           n.hull <- 1
                           if(missing(col.hull)) col.hull <- 1
                           if(missing(lty.hull)) lty.hull <- 1
                           if(missing(lwd.hull)) lwd.hull <- 1
                           x.old <- x; y.old <- y
                           idx <- chull(x, y); x.hull <- x[idx]; y.hull <- y[idx]
                           for(i in 1:(length(x)/3)){
                             x <- x[-idx]; y <- y[-idx]
                             if((length(x)/n) < fraction){
                               return(cbind(x.hull, y.hull))
                             }
                             idx <- chull(x, y); x.hull <- x[idx]; y.hull <- y[idx];
                           }
                         }
                         if(missing(col.hull)) col.hull <- 1:n.hull
                         if(length(col.hull)) col.hull <- rep(col.hull, n.hull)
                         if(missing(lty.hull)) lty.hull <- 1:n.hull
                         if(length(lty.hull)) lty.hull <- rep(lty.hull, n.hull)
                         if(missing(lwd.hull)) lwd.hull <- 1
                         if(length(lwd.hull)) lwd.hull <- rep(lwd.hull, n.hull)
                         result <- NULL
                         for(i in 1:n.hull){
                           idx <- chull(x, y); x.hull <- x[idx]; y.hull <- y[idx]
                           result <- c(result, list(cbind(x.hull, y.hull)))
                           x <- x[-idx]; y <- y[-idx]
                           if(0 == length(x)) return(result)
                         }
                         result
                       } 
                     the_matrix <- matrix(data = c(data$x, data$y), ncol = 2)
                     setNames(data.frame(plothulls_(the_matrix, fraction = prop)),
                              nm = c("x", "y"))},
                   required_aes = c("x", "y"))


stat_bag <- function(mapping = NULL, data = NULL, geom = "polygon",
                     position = "identity", na.rm = FALSE, show.legend = NA,
                     inherit.aes = TRUE, prop = 0.5, alpha = 0.3, ...) {
  ggplot2::layer(
    stat = StatBag, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, prop = prop, alpha = alpha, ...)
  )}


# biplot with alpha-bags


ggbiplot(prc, axis.type = "predictive", axis.percents = T) +
  theme_biplot() +  theme_dark() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  stat_bag(aes(fill=response_cat),prop=.9,alpha=.5) +  #prop is the percentage for the alpha-bag
  scale_fill_brewer(palette = "RdYlBu", na.translate=F) +
  geom_cols_axis(aes(label = name, center = center,scale=scale, label_dodge = -0.04)) +
  scale_colour_brewer(palette = "RdYlBu") + 
  guides(fill=guide_legend(title = "SOC", reverse=T)) 




######################################
#
#     END of Script
#
######################################
