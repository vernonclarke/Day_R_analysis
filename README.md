## R analysis for Day et al., 2023

Recreates graphical outputs and statistical analyses in the manuscript:

Day M., Belal M., Surmeier W. C., Melendez A., Wokosin, D., Tkatch T., Clarke V. R. J. and Surmeier D. J.

State-dependent GABAergic regulation of striatal spiny projection neuron excitability (2023)

## Table of Contents
- [Initial Set Up](#initial-set-up)
- [Performing Analysis](#performing-analysis)
- [Functions](#functions)
  - [WBplot](#wbplot)
  - [R2calc](#r2calc)
- [References](#references)
 
## Initial Set Up

  The analyses were conducted in the R graphical user interface (GUI): R version 4.3.1 – 'Beagle Scouts'. 

  R can be downloaded [here](https://www.R-project.org/). If you prefer to work with `RStudio`, it can be downloaded [here](https://posit.co/products/open-source/rstudio/). 

  Only the R console was used for analysis. This code should work in `RStudio` although this has not been explicitly tested. 

  All the required code is included in the R file: `graphs and statistics.R`
  
## Performing Analysis

  In order for the R code to work, it is necessary to load various packages within the R environment.

  The following steps 1-3 should be executed prior to any analysis. 
  
  1. **Load Packages**

     This code will check if the relevant packages are already installed in R. If not, a prompt will appear requiring the user to select a download site for installing these repositories. This installation is only required once for new versions of R. Once/if installed, the packages are then loaded into the current R environment.

     The required packages are [`svglite`](https://cran.r-project.org/web/packages/svglite/index.html) and [`lme4`](https://www.rdocumentation.org/packages/nlme/versions/3.1-163/topics/lme) .

     **Nb** as the function `load_required_packages` checks if the packages are in the library, installs them if they are not AND subsequently loads them into the current R environment, it must be run everytime analysis is performed.


```R
rm( list=ls(all=TRUE ) )
load_required_packages <- function(packages){
	new.packages <- packages[!(packages %in% installed.packages()[,'Package'])]
	if (length(new.packages)) install.packages(new.packages)
	invisible(lapply(packages, library, character.only=TRUE))
}	  
required.packages <- c('svglite', 'lme4')
load_required_packages(required.packages) 
```
    
2. **Set Initial Settings**

   Set the working directory:

```R
mypath <- '/yourpath/Day_R_analysis/figs and csv files' 
wd <- paste0(getwd(), mypath)
setwd(wd)
``` 
3. **Required Custom Functions**
   
   These custom-written functions are required to make the graphs and perform the statistical analyses.

```R
WBplot <- function(data, wid=1, cap=0.5, xlab = 'membrane potential (mV)', ylab = 'PSP amplitude (mV)', xrange=c(-70,-50), yrange=c(-10,15), lwd=0.8, type=6) {
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
	    
# is random effects model is singular?
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
	r2_values <- R2calc(formula, data)
	
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
	WBplot(data, wid=wid, cap=cap, xlab=xlab, ylab=ylab, xrange=xrange, yrange=yrange, lwd=lwd, type=type)
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

R2calc <- function(formula, data) {
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

plot.error.bars <- function(X, Y, color, lwd, xrange, yrange) {
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
	plot.error.bars(X, Y, color='indianred', lwd=lwd, xrange=xrange, yrange=yrange)

	# Box12 Data
	X <- box13[2:4,3]
	Y <- box12[2:4,3]
	plot.error.bars(X, Y, color='black', lwd=lwd, xrange=xrange, yrange=yrange)
    
	points(data13$y[data13$x == 2], data12$y[data12$x == 2],pch=20, col='indianred')
	points(data13$y[data13$x == 3], data12$y[data12$x == 3], pch=20, col='black')
}
```

4. **Data Analysis**

   Having run all the previous code, the following code performs the analysis and makes the graphs for all the electrophysiological experiments found in the manuscript.

```R
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
dev.new(width=6 ,height=4, noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data4, ylab='PSC amplitude (pA)', yrange=c(-20,25), p.cex=0.6)

fun.plot(data5, yrange=c(-20,25), p.cex=0.6)

# FIG2F
dev.new(width=9 ,height=4, noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,3), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data1, p.cex=0.6)

fun.plot(data2, p.cex=0.6)

fun.plot(data3, p.cex=0.6)

# FIG2GJ
dev.new(width=6 ,height=4, noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data6, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(-70, -55), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)

fun.plot(data7, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(-70, -55), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)

# FIG2GJ: statistcal tests
wilcox.f(data=data6, group1=1, group2=2)
wilcox.f(data=data7, group1=1, group2=2)
```
	
```R
# Fig3C
data8 <- import.fun('data3CA')
data9 <- import.fun('data3CB')

dev.new(width=6 ,height=4, noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data8, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)

fun.plot(data9, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)

# Fig3C: statistcal tests
wilcox.f(data=data8,group1=1, group2=2)
wilcox.f(data=data9,group1=1, group2=2)

# Fig3F
data3F <- import.fun('data3F')
data10 <- data3F[, c('s', 'x', 'y1')]; colnames(data10)[3] <- 'y'
data11 <- data3F[, c('s', 'x', 'y2')]; colnames(data11)[3] <- 'y'

dev.new(width=6 ,height=4, noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data10, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)
fun.plot(data11, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6, regression=FALSE, silent=TRUE)

# Fig3F: statistcal tests; performing Mann U as not enough values for paired (n = 5 pairs is required)
wilcox.f(data=data10,group1=1, group2=2, paired=FALSE)

# nb result IS, in fact, identical to the previous one (simply because pairs in data10 and data11 go in identical directions)
wilcox.f(data=data11,group1=1, group2=2, paired=FALSE)

```
```R
# simple eg to illustrate why p values in Fig3F are identical:
test1 <- data11
# Set a seed for reproducibility
set.seed(42) 

# Replace y with normally distributed numbers
test1$y <- rnorm(nrow(test1))

# Add 10 if x = 1, and 20 if x = 2
test1$y <- test1$y + ifelse(test1$x == 1, 10, ifelse(test1$x == 2, 20, 0)) # all values in level 2 are larger than any value in level 1

wilcox.f(data=test1,group1=1, group2=2, paired=FALSE)
```
 
```R
# Fig4EF 
data4E <- import.fun('data12')
data12 <- data4E[, c('s', 'x', 'y1')]; colnames(data12)[3] <- 'y'
data13 <- data4E[, c('s', 'x', 'y2')]; colnames(data13)[3] <- 'y'

dev.new(width=9 ,height=4, noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
fun.plot(data12, yrange=c(0,35), xrange=c(0.5,3.5), xlab='', ylab='', wid=0.2, cap=0.1, amount=0, p.cex=0.6, regression=FALSE, silent=TRUE)
fun.plot(subset(data13, x != 1), yrange=c(0,0.25), xrange=c(0.5,3.5), xlab='', ylab='', wid=0.2, cap=0.1, amount=0, p.cex=0.6, regression=FALSE, silent=TRUE)

# stats for Fig4E (dataset combined)
wilcox.f(data=data12,group1=2, group2=3)
wilcox.f(data=data13,group1=2, group2=3)
```

```R
# simple eg to illustrate why p values in Fig4EF are identical:
test2 <- data12
# Replace y with numbers generated from a normal distribution
set.seed(42) # for reproducibility
test2$y <- rnorm(nrow(test2))

# Add 10, 20, or 30 to y based on the value of x
test2$y <- test2$y + ifelse(test2$x == 1, 10, ifelse(test2$x == 2, 20, 30)) # direction of change for each pair is the same

wilcox.f(data=test2,group1=2, group2=3)
```
```R	
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
```

5. **Save all the figures in the same directory as the raw data**

   Saving graphs: if `plotsave <- TRUE` then any generated graphs are saved as svgs into working folder 

```R
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
```	

## Functions

### WBplot

The function, `WBplot` is a customised function to create whisker-and-box plots. 

This is a standardized way of displaying the distribution of data based on a five-number summary:

- Minimum: The smallest data point, excluding any outliers.
- First quartile (q1): The data point below which 25% of the data fall.
- Median (q2 or second quartile): The data point that divides the data into two halves. 50% of the data fall below the median, and 50% of the data fall above it.
- Third quartile (q3): The data point below which 75% of the data fall.
- Maximum: The largest data point, excluding any outliers.
- Outliers are defined as values more extreme than q1 - 1.5 * iqr  for lower and q3 + 1.5 * iqr for upper bound limits where the inter-quartile range is defined as iqr = q3 - q1.

Any outliers are removed and the default setting for calculating the quartiles is `type = 6`.

In R's `quantile` function, there are 9 types of quantile algorithms, named type 1 to type 9. These algorithms use different methods to calculate specified quantiles (for more information [see](https://doi.org/10.2307/2684934)). 

- Type 1: Inverse of the empirical distribution function.
- Type 2: Similar to type 1 but with averaging at discontinuities.
- Type 3: SAS definition: nearest even order statistic.
- Type 4: Linear interpolation of the empirical cdf.
- Type 5: Piecewise linear function where the knots are the values of order statistics.
- Type 6: Linear interpolation of the expectations based on the order statistics (default).
- Type 7: Linear interpolation of the modes based on the order statistics.
- Type 8: Linear interpolation between the points that capture the α percent and 1-α percent of the data.
- Type 9: Linear interpolation of the approximate medians for order statistics.   

**Nb** The R function, `boxplot` is not used to make whisker-and-box plots because this is not prefered method used by most graphics software. The native R function calls `boxplot.stats` which, in turn, calls `stats::fivenum` to calculate the medium iqr, minimum and maximum based on Tukey's five-number summary definition.

John Tukey's 'hinges', which are used in his five-number summary and for drawing boxplots, are similar to quartiles but can be calculated in a way that's slightly different from any of the standard quantile methods in R. Tukey's original definition involved using the median to split the data set and then finding the median of the lower and upper halves. If the data set or data half contains an odd number of points, the median is included in both halves.

How Tukey's hinges are usually computed:

- The lower hinge is the median of the lower half of the data set (not including the overall median if the number of data points is odd).
- The upper hinge is the median of the upper half of the data set (again not including the overall median if the number of data points is odd).

This method is somewhat akin to R's 'Type 1' method for calculating quantiles, also known as the 'inverted empirical distribution function'. 

The `quantile` function in R with the option `type=7` in R's default `boxplot`, the applied method is close to, but not exactly the same as, Tukey's original 'hinge' method. 

`GraphPad Prism` default seems to calculate quartiles using the method which corresponds to 'Type 6' in R's `quantile` function. For further information [see](https://www.graphpad.com/support/faq/how-prism-computes-percentiles/)

A comparison of `WBplot` with 'Type 6' , 'Type 7' and R's native `boxplot` function is illustrated by running the following code:

```R
data4E <- import.fun('data12')
data12 <- data4E[, c('s', 'x', 'y1')]; colnames(data12)[3] <- 'y'

dev.new(width=9 ,height=4, noRStudioGD=TRUE)
par(mar=c(1, 1, 1, 1), mfrow=c(1,3), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)

WBplot(data=data12, wid=0.2, cap=0.1, xlab = '', ylab = '', xrange=c(0.5,3.5), yrange=c(0,35), lwd=0.8)
# Add text to the plot
text(x = 2.5, y = 0, labels = 'WBplot; type=6')

WBplot(data=data12, wid=0.2, cap=0.1, xlab = '', ylab = '', xrange=c(0.5,3.5), yrange=c(0,35), lwd=0.8, type=7)
text(x = 2.5, y = 0, labels = 'WBplot; type=7')

boxplot(y ~ x, data = data12, main = '', xlab = '', ylab = '', ylim=c(0,35), xlim=c(0.5,3.5), boxwex=0.35, lwd=0.8, col='white', border='black', frame=FALSE)
text(x = 2.5, y = 0, labels = 'boxplot')
```


### Summary

- **The native R function, `boxplot` calculates whisker-and-box plots based on Tukey's original 'hinges' method by calling `stats::fivenum`**.
- **The function `WBplot` calculates quartiles using R function `quantile`. This can be set to type = 1 to 9**.
- **The default in `WBplot` is type = 6 which should produce similar results to `GraphPad Prism`**.
- **For results closer (but not identical) to Tukey's 'hinges' method / R's native `boxplot`, set type = 7**.

### Random Mixed-Effects Model

**Linear Regression** is performed using the package `lmer`. The function determines whether the fit of the model is singular. If not then it fits by random mixed-effects model. In `lmer` terminology, the formula for this is **y ~ x + (1|s)**. This formula specifies how the dependent variable 'y' is modeled in relation to the independent (or fixed-effect) predictor variable 'x' and the random effect of the subject 's'. 

- y: This is the dependent variable you are trying to model or predict.
-  ~: The tilde separates the dependent variable from the independent variables and random effects.
- x: This is the independent (or fixed-effect) variable. The model will estimate how 'y' varies with 'x'.
- +: The plus sign indicates that you are including more terms in the model.
- (1|s): This is a random intercept for subject s. In other words, each subject is allowed to have its own intercept value of y that is randomly distributed around the overall mean of y.

### Linear Model
**y ~ x**: the formula specifies how the dependent variable y is modeled in relation to the predictor variable x. 

Fits to the random mixed-effects models may be singular if:

- There is some redundancy in random effects: the model is too complex for the data with random effects that do not contribute much variance. As an example, the random intercepts for each level of subject, s may be very similar. This can happen if there are too few levels in s or if the data within each level of s does not vary much.

- The variance of the random effects or the residual variance is estimated to be near zero. This suggests that the random effect might not be necessary, as it doesn't explain a significant amount of variability in the data.

If fits of **y ~ x + (1|s)** are singular, the function simplifies the model to linear regression **y ~ x**.

### R2calc

The function `R2calc` calculates the appropriate value of $R^2$.

In the context of linear mixed-effect models, $R^2$ can be a bit more complex to define and interpret than in standard linear regression. 

There are actually two commonly reported $R^2$ values for linear mixed models:
	
- **Marginal $R^2$**: Represents the variance explained by the fixed effects alone.
- **Conditional $R^2$**: Represents the variance explained by both the fixed and random effects.

The conditional $R^2$ is always equal to or larger than the marginal $R^2$ since it also includes the variance explained by the random effects.

It is suitable for datasets where the first level corresponds to an independent (or fixed-effect) variable and the second level to some grouping/clustering factor (such as subjects with repeated measurements). 

The method for calculating $R^2$ in linear mixed-effect models is nicely summarised in [Nakagawa et al., 2017](http://dx.doi.org/10.1098/rsif.2017.0213) and is reproduced here.

For Normal/Gaussian error distributions, the model is specified as follows:

$$ y_{ij} = b_0 + \sum_{h=1}^{p} b_h x_{hij} + a_i + \epsilon_{ij} $$

where:
- $y_{ij}$ is the $j^{th}$ observation of the $i^{th}$ subject.
- $b_0$ is the intercept (or grand mean).
- $b_h$ is the fixed effect coefficient for the $h^{th}$ predictor.
- $x_{hij}$ is the $j^{th}$ value for the $i^{th}$ subject for the $h^{th}$ predictor.
- $a_i$ is an subject-specific effect, assumed to be normally distributed in the population with mean 0 and variance $\sigma^2_a$.
- $\epsilon_{ij}$ is an observation-specific residual, assumed to be normally distributed in the population with mean 0 and variance $\sigma^2_\epsilon$.

For this model, two types of $R^2$ can be defined:

**Marginal $R^2$ or $R^2_{LMM(m)}$**:

$$ R^2_{LMM(m)} = \frac{\sigma^2_f}{\sigma^2_f + \sigma^2_a + \sigma^2_\epsilon} $$

This represents the proportion of the total variance explained by the fixed effects.

**Conditional $R^2$ or $R^2_{LMM(c)}$**:

$$ R^2_{LMM(c)} = \frac{\sigma^2_f + \sigma^2_a}{\sigma^2_f + \sigma^2_a + \sigma^2_\epsilon} $$

This represents the proportion of variance explained by both fixed and random effects.

Where:
- $\sigma^2_f$ is the variance explained by fixed effects, calculated as $\text{var} \left( \sum_{h} b_h x_{hij} \right)$.

In cases where the fit reverts to linear regression (i.e. linear mixed-effect fits of **y ~ x + (1|s)** are singular), this function returns the standard value of $R^2$ and the adjusted value $R^2_{adj}$.

**$R^2$ for Linear Model or $R^2_{LM}$**:

$$ R^2_{LM} = \frac{\sigma^2_f}{\sigma^2_t} $$

where $\sigma^2_f$ is the variance of the fitted (predicted) values and $\sigma^2_t$ is the total variance of the response variable ($y$). 

This represents the proportion of the total variance in the response variable that is explained by the linear model.

**Adjusted $R^2$ for Linear Model or $R^2_{adj, LM}$**:

$$ R^2_{adj, LM} = 1 - \left( \frac{(1 - R^2_{LM}) \times (n - 1)}{n - p - 1} \right) $$

where $R^2_{LM}$ is the $R^2$ for the linear model, $n$ is the number of observations, and $p$ is the number of predictors excluding the intercept. 

Adjusted R-squared accounts for the number of predictors in the model and is especially useful when comparing models with different numbers of predictors.


## References

[Bates D, Maechler M, Bolker B, Walker S (2015). Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.](https://doi.org/10.18637/jss.v067.i01)

[Hyndman R. J. and Fan Y (1996), Sample Quantiles in Statistical Packages. The American Statistician, Vol. 50(4), 361-365.](https://doi.org/10.2307/2684934)

[Nakagawa S, Johnson P. C. D. and Schielzeth H (2017). The coefficient of determination Rsqr and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. J. R. Soc. 14 14(134):20170213.](http://dx.doi.org/10.1098/rsif.2017.0213)

[R Core Team (2023). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.](https://www.R-project.org/)

[Wickham H, Henry L, Pedersen T, Luciani T, Decorde M, Lise V (2023). _svglite: An 'SVG' Graphics Device_. R package version 2.1.2.](https://CRAN.R-project.org/package=svglite)



This code was written by Vernon Clarke.

The provided code was executed on a `MacBook M2 pro 32GB`. 

For queries related to this repository, please open an [issue](https://github.com/vernonclarke/Day_R_analysis/issues) or [email](WOPR2@proton.me) me directly





