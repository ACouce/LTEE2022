# .. This script builds and plots a linear model for mutation prevalence

# .. build linear model
Rblocks[Rblocks==-Inf]<-NA
model50 <- lm(nobs50 ~ Rblocks$R + Rblocks$M + Rblocks$K + Rblocks$size, data = Rblocks)

# .. print summary statistics
print(summary(model50))

# .. load useful 3D graphics library
library("plot3D")

# .. create custom palette
mycol <- rgb(190, 190, 190, max = 255, alpha = 60, names = "blue50")

# .. set the x, y, and z variables for 3D plot
x <- as.vector(Rblocks$R)
y <- -as.vector(Rblocks$size)
z <- as.vector(Rblocks$nobs50)

# .. compute the linear regression 
fit <- lm(z ~ x + y)

# .. create grid from the x and y values (min to max) and predict values for every point, which will produce the regression plane
grid.lines = 10
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), nrow = grid.lines, ncol = grid.lines)

# .. create the fitted points for droplines to the surface
forcol <- x
fitpoints <- predict(fit)

# .. make list of point sizes scaled to represent target size
sizes<-((Rblocks$size)/1600)^0.7

# .. initialize graphical output
png(file="3d.png", width = (15.63*0.9)*(2/1.75), height = (8.54)/1.75, units = 'in', res = 300)
par(mfrow=c(1,2), family="sans", cex=1)

# .. make the 3D scatter plot
scatter3D(x, y, z, pch = 19,  colvar=forcol, col= ramp.col(col = c("dodgerblue3","#7F648F", "#7F648F", "#ff3f0f", "#ff3f0f"), n = 48), cex = sizes, ticktype = "detailed", colkey = FALSE, xlab = ' ', ylab = ' ', zlab = ' ', theta = -138, phi = 30, bty="b", cex.lab=2, cex.axis=0.5, lwd=1.5, type='h', surf = list(x = x.pred, y = y.pred, z = z.pred, facets = TRUE,  col=mycol, border='darkgrey', lwd=0.3), yaxt="n")

# .. close graphical output
dev.off()
