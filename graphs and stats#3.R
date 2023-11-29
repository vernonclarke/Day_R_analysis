##############################################################################################################
# This code will check if the relevant packages are already installed in R. 
# If not, a prompt will appear requiring the user to select a download site for installing these repositories. 
# This installation is only required once for new versions of R. 
# Once/if installed, the packages are then loaded into the current R environment.

# The output from running each analysis is included (as commented out code after the routines that generated it.

rm( list=ls(all=TRUE ) )
load_required_packages <- function(packages){
	new.packages <- packages[!(packages %in% installed.packages()[,'Package'])]
	if (length(new.packages)) install.packages(new.packages)
	invisible(lapply(packages, library, character.only=TRUE))
}	  
required.packages <- c('svglite', 'lme4')
load_required_packages(required.packages) 
##############################################################################################################
# set the working directory
mypath <- '/documents/Repositories/Day_R_analysis/figs and csv files' 
wd <- paste0(getwd(), mypath)
setwd(wd)
##############################################################################################################
# Custom functions
custom_boxplot <- function(data, wid=1, cap=0.5, xlab = 'membrane potential (mV)', 
ylab = 'PSP amplitude (mV)', xrange=c(-70,-50), yrange=c(-10,15), 
lwd=0.8, type=6) {
	x <- data$x
	y <- data$y
	unique_x <- unique(data$x)
	xrange <- xrange + c(-wid, wid)
	plot(1, type='n', ylim=yrange, xlim=xrange, xlab=xlab, ylab=ylab, xaxt='n', yaxt='n', bty='n', lwd=lwd)
	
	for (i in 1:length(unique_x)) {
	current_x <- unique_x[i]
	d <- data$y[data$x == current_x]
	
	q1 <- quantile(d, probs=0.25, type=type)
	q3 <- quantile(d, probs=0.75, type=type)
	iqr <- q3 - q1  # Calculate IQR
	
	lower_bound <- q1 - 1.5 * iqr  # Lower bound for outliers
	upper_bound <- q3 + 1.5 * iqr  # Upper bound for outliers
	
	# Exclude outliers
	d_filtered <- d[d >= lower_bound & d <= upper_bound]
	
	median_val <- median(d_filtered)
	min_val <- min(d_filtered)
	max_val <- max(d_filtered)
	
	rect(current_x - wid, q1, current_x + wid, q3, col='white', lwd=lwd)
	segments(current_x, q1, current_x, min_val, lwd=lwd)
	segments(current_x, q3, current_x, max_val, lwd=lwd)
	segments(current_x - cap, min_val, current_x + cap, min_val, lwd=lwd)
	segments(current_x - cap, max_val, current_x + cap, max_val, lwd=lwd)
	segments(current_x - wid*1.1, median_val, current_x + wid*1.1, median_val, col='black', lwd=3*lwd)
	}
	axis(1, at=unique_x, labels=unique_x)
	axis(2)
}
	    
# if random effects model is singular
isSingular.fun <- function(formula, data){
	mod <- suppressMessages(lmer(formula=formula, data=data))
	isSingular(mod)
	} 
	
fun.plot = function(data, wid=1, cap=0.5, xlab = 'membrane potential (mV)', ylab = 'PSP amplitude (mV)', xrange=c(-70,-50), yrange=c(-10,15), lwd=0.8, amount=0.5, p.cex=0.25, type=6,  regression=TRUE, silent=FALSE){
	
	# Fit the model using lmer
	# model_lmer <- lmer(y ~ x + (1|s))
	x <- data$x
	y <- data$y
	s <- data$s
	formula  <- y ~ x + (1|s)
	mixed <- !isSingular.fun(formula, data)
	mod <-  if (mixed) lmer(y ~ x + (1|s)) else lm(y ~ x)
	formula <- if (!mixed) y ~ x else formula
	if (mixed){
		fixed_effects <- fixef(mod)
		m <- fixed_effects[[2]]
		c <- fixed_effects[[1]]
	}else{
		coeffs <- coef(mod)
		m <- coeffs[[2]]
		c <- coeffs[[1]]
	}
	r2_values <- R2.calculator(formula, data)
	
	if (!silent){
		cat('model is ', format(formula), '\n')

		print(summary(mod))

		# x intercept
		cat('x intercept is ', format(-c/m), 'mV', '\n')
	    
		if (mixed){
			cat('rsqr (marginal) ', format(r2_values[1]), ' rsqr (conditional) ', format(r2_values[2]))
		}else{
			cat('rsqr (multiple) ', format(r2_values[1]), ' rsqr (adjusted) ', format(r2_values[2]))
		}
	}
	# Add a column to data for jittered x-values
	set.seed(42)
	data$x_jitter <- jitter(data$x, amount=amount)
	
	# Boxplot
	custom_boxplot(data, wid=wid, cap=cap, xlab=xlab, ylab=ylab, xrange=xrange, yrange=yrange, lwd=lwd, type=type)
	# Plot individual data points with reduced jitter, reduced size, and unfilled circles without x and y axes
	points(data$x_jitter, data$y, pch=19, bg='transparent', col='darkgray', lwd=lwd/2, cex=p.cex)
	# Connect data points within subjects with gray dotted lines
	line=TRUE
	if (line){
		subjects <- unique(data$s)
		for(subj in subjects){
  			subset_data <- data[data$s == subj, ]
  			lines(subset_data$x_jitter, subset_data$y, col='darkgray', lwd=lwd, lty=3)  # lty=2 for dotted line
		}
	}
	
	if (regression){
		# Predict y values using the model for new data with matching grouping factor levels
		y_pred <- m * unique(x) + c
		# Add the line of best fit
		lines(unique(x), y_pred, col='black', lwd=lwd, lty=1)
	}
	
	# list(reversal=-c/m, r2_values =r2_values)	
}

R2.calculator <- function(formula, data) {

	# Convert the model formula to a string
	formula_str <- deparse(formula)
	
	# Fit the model
	if (grepl('\\|', formula_str)) {
		# Model has random effects (intercepts or slopes)
		model <- lmer(formula, data = data)
		includeRandomEffect <- TRUE
	} else {
		# Model is a simple linear model
		model <- lm(formula, data = data)
		includeRandomEffect <- FALSE
	}
	
	# Variance explained by fixed effects (sigma^2_f)
	y_pred_fixed <- predict(model, re.form = NA)
	varFixed <- var(y_pred_fixed)
	
	if (includeRandomEffect) {
		# Variance for random effect (sigma^2_a)
		varRandom <- sum(sapply(VarCorr(model), function(v) v[1]))
		
		# Residual variance (sigma^2_epsilon)
		varResid <- sigma(model)^2
		
		# Calculate R-squared
		R2_marginal <- varFixed / (varFixed + varRandom + varResid)
		R2_conditional <- (varFixed + varRandom) / (varFixed + varRandom + varResid)
		
		return(list(marginal = R2_marginal, conditional = R2_conditional))

	} else {
		# For lm model, R-squared is simply the variance of the predicted values 
		# over the total variance
		totalVar <- var(data$y)
		R2 <- varFixed / totalVar
	
		# Calculate adjusted R-squared manually
  			n <- nrow(data)
  			p <- length(coef(model)) - 1  # number of predictors, excluding intercept
		adjR2 <- 1 - (1 - R2) * (n - 1) / (n - p - 1)

  			return(list(R2 = R2, adjustedR2 = adjR2))
	}
}


# simple function to import data from a ''csv' file
#  if NA is zero imports all exlcude excludes those subjects s
import.fun <- function(name, exclude=NA){
  	df <- read.csv(paste0(name, '.csv'))
  
  	# Exclude rows based on 's' values
  	if (!is.na(exclude[1])) {  # Check if the first element of 'exclude' is not NA
    		df <- df[!(df$s %in% exclude), ]
  	}
  
  	# Assuming the first two columns are always 's' and 'x'
  	fixed_colnames <- c('s', 'x')
  
  	# Check the number of remaining columns after 's' and 'x'
  	num_y_cols <- ncol(df) - length(fixed_colnames)
  
  	if (num_y_cols == 1) {
    		y_colnames <- 'y'
  	} else {
    		y_colnames <- paste0('y', 1:num_y_cols)
  	}
  
  	colnames(df) <- c(fixed_colnames, y_colnames)
  	return(df)
}

# Wrapper to perform the Wilcoxon Signed-Rank Test etc
wilcox.f <- function(data, group1, group2, paired=TRUE, alternative='two.sided', exact=NULL){
    	# Extract data for the first group
    	x <- data[data$x == group1,][order(data[data$x == group1,]$s),]$y
    	# Extract data for the second group
    	y <- data[data$x == group2,][order(data[data$x == group2,]$s),]$y
    	# Perform the Wilcoxon Signed-Rank Test
    	wilcox.test(x, y, paired = paired, alternative = alternative, exact = exact)
}

# Functions for FigS1	
fun.plot.S1 = function(){
	plot(dataS1$'A+B', dataS1$'C', xlab = 'linear prediction (mV)', ylab = 'actual combination', bty='n', pch=20, col='black', xlim=c(0,40), ylim=c(0,40))

	# Define the points for the y=x line
	xline = yline = seq(0, 40, 0.1)

	# Add shading below the line using polygon
	polygon(c(0, xline, 40), c(0, yline, 0), col='lightgray', border=NA)

	# Overlay the y=x line on top of the shaded region
	lines(xline, yline, col='black', lwd=1)

	# Replot your data points on top to ensure they're not covered by the polygon
	points(dataS1$'A+B', dataS1$'C', pch=20, col='black')
}

output.fun <- function(data, type=6, MAD = FALSE){
	unique_x <- unique(data$x)
	out <- sapply(1:length(unique_x), function(ii){
		current_x <- unique_x[ii]
        d <- data$y[data$x == current_x]
        
        q1 <- quantile(d, probs=0.25, type=type)
        q3 <- quantile(d, probs=0.75, type=type)
        iqr <- q3 - q1  # Calculate IQR
        
        lower_bound <- q1 - 1.5 * iqr  # Lower bound for outliers
        upper_bound <- q3 + 1.5 * iqr  # Upper bound for outliers
        
        # Exclude outliers
        d_filtered <- if (MAD) d else d[d >= lower_bound & d <= upper_bound] # do NOT remove outliers for MAD

        median_val <- median(d_filtered)
        min_val <- min(d_filtered)
        max_val <- max(d_filtered)
        
        
        if (MAD){
        	mad_value <- mad(d_filtered, constant = 1)
        	c(min_val[[1]], median_val-mad_value, median_val, median_val+mad_value, max_val[[1]])	
        }else{
        	c(min_val[[1]], q1[[1]], median_val, q3[[1]], max_val[[1]])	
        }
                     
    })
    if (MAD){
    	rownames(out) <- c('min','lm','median','um','max')
    }else{
    	rownames(out) <- c('min','q1','median','q3','max')
    }
    colnames(out) = unique_x
    return(out)
}

plot_error_bars <- function(X, Y, color, lwd, xrange, yrange) {
	x_q1 <- X[1]
	x_median <- X[2]
	x_q3 <- X[3]
	
	y_q1 <- Y[1]
	y_median <- Y[2]
	y_q3 <- Y[3]
	
	wid.x <- diff(xrange) / 50
	wid.y <- diff(yrange) / 50
	
	# Plot median points
	points(x_median, y_median, pch=19, col=color)
	
	# Error bars for x
	segments(x_q1, y_median, x_q3, y_median, col=color, lwd=lwd)
	segments(x_q1, y_median - 0.5*wid.y, x_q1, y_median + 0.5*wid.y, col=color, lwd=lwd)
	segments(x_q3, y_median - 0.5*wid.y, x_q3, y_median + 0.5*wid.y, col=color, lwd=lwd)
	
	# Error bars for y
	segments(x_median, y_q1, x_median, y_q3, col=color, lwd=lwd)
	segments(x_median - 0.5*wid.x, y_q1, x_median + 0.5*wid.x, y_q1, col=color, lwd=lwd)
	segments(x_median - 0.5*wid.x, y_q3, x_median + 0.5*wid.x, y_q3, col=color, lwd=lwd)	  
}

fun.plot2 = function(data12, data13){
	plot(NULL, xlim=xrange, ylim=yrange, xlab='x', ylab='y', type='n', bty='n')

	box12 <- output.fun(data12, MAD=TRUE)
	box13 <- output.fun(data13, MAD=TRUE)

	# Box13 Data
	X <- box13[2:4,2]
	Y <- box12[2:4,2]
	plot_error_bars(X, Y, color='indianred', lwd=lwd, xrange=xrange, yrange=yrange)

	# Box12 Data
	X <- box13[2:4,3]
	Y <- box12[2:4,3]
	plot_error_bars(X, Y, color='black', lwd=lwd, xrange=xrange, yrange=yrange)
    
	points(data13$y[data13$x == 2], data12$y[data12$x == 2],pch=20, col='indianred')
	points(data13$y[data13$x == 3], data12$y[data12$x == 3], pch=20, col='black')
}
##############################################################################################################
# Fig 2
# import data sets from stored csv

data2F <- import.fun('data2F')
data1 <- data2F[, c('s', 'x', 'y1')]; colnames(data1)[3] <- 'y'
data2 <- data2F[, c('s', 'x', 'y2')]; colnames(data2)[3] <- 'y'
data3 <- data2F[, c('s', 'x', 'y3')]; colnames(data3)[3] <- 'y'

data2C <- import.fun('data2C')
data4 <- data2C[, c('s', 'x', 'y1')]; colnames(data4)[3] <- 'y'
data5 <- data2C[, c('s', 'x', 'y2')]; colnames(data5)[3] <- 'y'

data6 <- import.fun('data2G')
data7 <- import.fun('data2J')

# FIG2C
dev.new(width=6 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data4, ylab='PSC amplitude (pA)', yrange=c(-20,25), p.cex=0.6)

# model is  y ~ x + (1 | s) 
# Linear mixed model fit by REML ['lmerMod']
# Formula: y ~ x + (1 | s)

# REML criterion at convergence: 133.7

# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -1.60033 -0.65121  0.03389  0.45512  2.00762 

# Random effects:
#  Groups   Name        Variance Std.Dev.
#  s        (Intercept)  4.723   2.173   
#  Residual             10.125   3.182   
# Number of obs: 25, groups:  s, 5

# Fixed effects:
#             Estimate Std. Error t value
# (Intercept)   53.880      5.524   9.754
# x              0.964      0.090  10.711

# Correlation of Fixed Effects:
#   (Intr)
# x 0.978 
# x intercept is  -55.89212 mV 
# rsqr (marginal)  0.7652421  rsqr (conditional)  0.8399144

fun.plot(data5, yrange=c(-20,25), p.cex=0.6)

# model is  y ~ x 

# Call:
# lm(formula = y ~ x)

# Residuals:
#    Min     1Q Median     3Q    Max 
#  -2.27  -1.74  -0.27   1.29   5.20 

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 42.10000    3.58113   11.76 3.33e-11 ***
# x            0.70600    0.05928   11.91 2.57e-11 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 2.096 on 23 degrees of freedom
# Multiple R-squared:  0.8605,	Adjusted R-squared:  0.8544 
# F-statistic: 141.9 on 1 and 23 DF,  p-value: 2.566e-11

# x intercept is  -59.63173 mV 
# rsqr (multiple)  0.8604881  rsqr (adjusted)  0.8544223

# FIG2F
dev.new(width=9 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,3), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data1, p.cex=0.6)

# model is  y ~ x + (1 | s) 
# Linear mixed model fit by REML ['lmerMod']
# Formula: y ~ x + (1 | s)

# REML criterion at convergence: 95.1

# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -1.3959 -0.7732 -0.1343  0.8783  1.9779 

# Random effects:
#  Groups   Name        Variance Std.Dev.
#  s        (Intercept) 0.1413   0.3758  
#  Residual             2.2187   1.4895  
# Number of obs: 25, groups:  s, 5

# Fixed effects:
#             Estimate Std. Error t value
# (Intercept) 46.44000    2.55085   18.21
# x            0.75400    0.04213   17.90

# Correlation of Fixed Effects:
#   (Intr)
# x 0.991 
# x intercept is  -61.59151 mV 
# rsqr (marginal)  0.9261828  rsqr (conditional)  0.9306014

fun.plot(data2, p.cex=0.6)

# model is  y ~ x 

# Call:
# lm(formula = y ~ x)

# Residuals:
#    Min     1Q Median     3Q    Max 
#  -3.58  -0.60  -0.08   0.96   2.96 

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 43.32000    2.73269   15.85 7.15e-14 ***
# x            0.70400    0.04523   15.56 1.05e-13 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.599 on 23 degrees of freedom
# Multiple R-squared:  0.9133,	Adjusted R-squared:  0.9095 
# F-statistic: 242.2 on 1 and 23 DF,  p-value: 1.053e-13

# x intercept is  -61.53409 mV 
# rsqr (multiple)  0.9132883  rsqr (adjusted)  0.9095182

fun.plot(data3, p.cex=0.6)

# model is  y ~ x 

# Call:
# lm(formula = y ~ x)

# Residuals:
#    Min     1Q Median     3Q    Max 
#  -3.14  -0.90  -0.14   0.88   2.84 

# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 46.26000    2.38701   19.38 9.60e-16 ***
# x            0.75200    0.03951   19.03 1.42e-15 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.397 on 23 degrees of freedom
# Multiple R-squared:  0.9403,	Adjusted R-squared:  0.9377 
# F-statistic: 362.3 on 1 and 23 DF,  p-value: 1.42e-15

# x intercept is  -61.51596 mV 
# rsqr (multiple)  0.9403001  rsqr (adjusted)  0.9377044

# FIG2GJ
dev.new(width=6 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data6, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(-70, -55), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)

fun.plot(data7, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(-70, -55), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)

# FIG2GJ: statistcal tests
wilcox.f(data=data6, group1=1, group2=2)

# 	Wilcoxon signed rank test with continuity correction

# data:  x and y
# V = 11, p-value = 0.4076
# alternative hypothesis: true location shift is not equal to 0

# Warning message:
# In wilcox.test.default(x, y, paired = paired, alternative = alternative,  :
#   cannot compute exact p-value with ties

wilcox.f(data=data7, group1=1, group2=2)

# 	Wilcoxon signed rank test with continuity correction

# data:  x and y
# V = 21, p-value = 0.03552
# alternative hypothesis: true location shift is not equal to 0

# Warning message:
# In wilcox.test.default(x, y, paired = paired, alternative = alternative,  :
#   cannot compute exact p-value with ties
##############################################################################################################
# Fig3C
data8 <- import.fun('data3CA')
data9 <- import.fun('data3CB')

dev.new(width=6 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data8, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)

fun.plot(data9, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)

# Fig3C: statistcal tests
wilcox.f(data=data8,group1=1, group2=2)

# 	Wilcoxon signed rank exact test

# data:  x and y
# V = 36, p-value = 0.007813
# alternative hypothesis: true location shift is not equal to 0

wilcox.f(data=data9,group1=1, group2=2)

# 	Wilcoxon signed rank exact test

# data:  x and y
# V = 21, p-value = 0.03125
# alternative hypothesis: true location shift is not equal to 0

# Fig3F
data3F <- import.fun('data3F')
data10 <- data3F[, c('s', 'x', 'y1')]; colnames(data10)[3] <- 'y'
data11 <- data3F[, c('s', 'x', 'y2')]; colnames(data11)[3] <- 'y'

dev.new(width=6 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data10, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)
fun.plot(data11, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)


# Fig3F: statistcal tests; performing Mann U as not enough values for paired (n = 5 pairs is required)
wilcox.f(data=data10,group1=1, group2=2, paired=FALSE)

# 	Wilcoxon rank sum exact test

# data:  x and y
# W = 16, p-value = 0.02857
# alternative hypothesis: true location shift is not equal to 0

wilcox.f(data=data11,group1=1, group2=2, paired=FALSE)

# 	Wilcoxon rank sum exact test

# data:  x and y
# W = 16, p-value = 0.02857
# alternative hypothesis: true location shift is not equal to 0

# nb result IS, in fact, identical to the previous one (simply because pairs in data10 and data11 go in identical directions with no ties)
##############################################################################################################
# Fig4EF 
data4E <- import.fun('data12')
data12 <- data4E[, c('s', 'x', 'y1')]; colnames(data12)[3] <- 'y'
data13 <- data4E[, c('s', 'x', 'y2')]; colnames(data13)[3] <- 'y'

dev.new(width=9 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data12, yrange=c(0,35), xrange=c(0.5,3.5), xlab='', ylab='', wid=0.2, cap=0.1, amount=0, p.cex=0.6, regression=FALSE, silent=TRUE)
fun.plot(subset(data13, x != 1), yrange=c(0,0.25), xrange=c(0.5,3.5), xlab='', ylab='', wid=0.2, cap=0.1, amount=0, p.cex=0.6, regression=FALSE, silent=TRUE)

# stats for Fig4E (dataset combined)
wilcox.f(data=data12,group1=2, group2=3)

# 	Wilcoxon signed rank exact test

# data:  x and y
# V = 0, p-value = 0.0009766
# alternative hypothesis: true location shift is not equal to 0

wilcox.f(data=data13,group1=2, group2=3)

# 	Wilcoxon signed rank exact test

# data:  x and y
# V = 0, p-value = 0.0009766
# alternative hypothesis: true location shift is not equal to 0

# nb result IS, in fact, identical to the previous one (simply because pairs in data12 and data13 go in identical directions with no ties)
##############################################################################################################
# data for figS1 
dataS1 <- read.csv('data14.csv')
colnames(dataS1) = c('A+B', 'C')

# FigS1
dev.new(width=4.5 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,1), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot.S1()

	
# Initial settings
lwd = 0.8; xrange = c(0,0.25); yrange = c(0,35)

# FigS1
dev.new(width=4.5 ,height=4,noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,1), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot2(data12, data13)
##############################################################################################################
# saving graphs: if plotsave <- TRUE then any generated graphs are saved as svgs into working folder.
# The final graphs were created from these svg files using Adobe Illustrator.
plotsave <- TRUE  
if (plotsave) {	
	svglite(paste0('Fig2F1 ', gsub(':', '-', Sys.time()), '.svg'), width=2,height=3.75, pointsize=10)
	fun.plot(data1, silent=TRUE)
	dev.off()

	svglite(paste0('Fig2F2 ', gsub(':', '-', Sys.time()), '.svg'), width=2,height=3.75, pointsize=10)
	fun.plot(data2, silent=TRUE)
	dev.off()

	svglite(paste0('Fig2F3 ', gsub(':', '-', Sys.time()), '.svg'), width=2,height=3.75, pointsize=10)
	fun.plot(data3, silent=TRUE)
	dev.off()

	svglite(paste0('Fig2C_1 ', gsub(':', '-', Sys.time()), '.svg'), width=2.5,height=2.75, pointsize=10)
	fun.plot(data4, ylab='PSC amplitude (pA)', yrange=c(-20,25), silent=TRUE)
	dev.off()
	
	svglite(paste0('Fig2C_2 ', gsub(':', '-', Sys.time()), '.svg'), width=2.5,height=2.75, pointsize=10)
	fun.plot(data5, ylab='PSC amplitude (pA)', yrange=c(-20,25), silent=TRUE)
	dev.off()
	
	svglite(paste0('Fig2G ', gsub(':', '-', Sys.time()), '.svg'), width=2.5,height=2.75, pointsize=10)
	fun.plot(data6, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(-70, -55), amount=0.05, regression=FALSE, silent=TRUE)
	dev.off()

	svglite(paste0('Fig2J ', gsub(':', '-', Sys.time()), '.svg'), width=2.5,height=2.75, pointsize=10)
	fun.plot(data7, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(-70, -55), amount=0.05, regression=FALSE, silent=TRUE)
	dev.off()

	svglite(paste0('Fig3C1 ', gsub(':', '-', Sys.time()), '.svg'), width=2.2,height=3.70, pointsize=10)
	fun.plot(data8, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)
	dev.off()

	svglite(paste0('Fig3C2 ', gsub(':', '-', Sys.time()), '.svg'), width=2.2,height=3.70, pointsize=10)
	fun.plot(data9, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)
	dev.off()

	svglite(paste0('Fig3F1 ', gsub(':', '-', Sys.time()), '.svg'), width=2.2,height=3.70, pointsize=10)
	fun.plot(data10, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)
	dev.off()

	svglite(paste0('Fig3F2 ', gsub(':', '-', Sys.time()), '.svg'), width=2.2,height=3.70, pointsize=10)
	fun.plot(data11, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)
	dev.off()

	svglite(paste0('Fig4E1 ', gsub(':', '-', Sys.time()), '.svg'), width=2.2,height=3.70, pointsize=10)
	fun.plot(data12, yrange=c(0,35), xrange=c(0.5,3.5), xlab='', ylab='', wid=0.3, cap=0.15, amount=0, p.cex=0.6, regression=FALSE, silent=TRUE)
	dev.off()
	
	svglite(paste0('Fig4E2 ', gsub(':', '-', Sys.time()), '.svg'), width=2.2,height=3.70, pointsize=10)
	fun.plot(subset(data13, x != 1), yrange=c(0,0.25), xrange=c(0.5,3.5), xlab='', ylab='', wid=0.3, cap=0.15, amount=0, p.cex=0.6, regression=FALSE, silent=TRUE)
	dev.off()
	
	svglite(paste0('Fig4E3 ', gsub(':', '-', Sys.time()), '.svg'), width=2.2,height=3.50, pointsize=10)
	fun.plot2(data12, data13)
	dev.off()
	
	svglite(paste0('FigS1 ', gsub(':', '-', Sys.time()), '.svg'), width=3.25,height=3.70, pointsize=10)
	fun.plot.S1()
	dev.off()
}
##############################################################################################################
