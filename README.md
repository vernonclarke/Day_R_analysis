### R analysis for Day et al., 2023

Recreates graphical outputs and statistical analyses in the manuscript 

## Table of Contents
- [Initial Set Up](#initial-set-up)
- [Performing Analysis](#performing-analysis)

- [Functions](#functions)
  - [custom_boxplot](#custom_boxplot)


 
## Initial Set Up

The final analysis and figures presented in the manuscript were generated using R. 

The analyses were conducted in the R graphical user interface (GUI):
  - R version 4.3.1 – "Beagle Scouts"
  - [R Statistical Software](https://www.R-project.org/)

  ### Setting up
  
  Only the R console was used for analysis. 
  
  If you prefer to work with `RStudio`, it can be downloaded [here](https://posit.co/products/open-source/rstudio/). The provided code should work although this has not been explicitly tested.
  
## Performing Analysis

  In order for the R code to work, it is necessary to load various packages within the R environment.

  The following steps 1-3 should be executed prior to any analysis. 
  
  1. **Load Packages**

     This code will check if the relevant packages are already installed in R. If not, a prompt will appear requiring the user to select a download site for installing these repositories. This installation is only required once for new versions of R. Once/if installeed, the packages are then loaded into the current R environment.

     The required packages are [`MuMIn`](https://cran.r-project.org/web/packages/MuMIn/index.html), [`svglite`](https://cran.r-project.org/web/packages/svglite/index.html) and [`lme4`](https://www.rdocumentation.org/packages/nlme/versions/3.1-163/topics/lme) .

     


```R
	rm( list=ls(all=TRUE ) )
	load_required_packages <- function(packages){
		new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
		if (length(new.packages)) install.packages(new.packages)
		invisible(lapply(packages, library, character.only=TRUE))
	}	  
	required.packages <- c('MuMIn', 'svglite', 'lme4')
	load_required_packages(required.packages) 
```
    
2. **Initial Settings**

     - set the working directory 

     - saving graphs: if `plotsave <- TRUE` then any generated graphs are saved as svgs into working folder 

```R
     mypath <- '/yourpath/Day_R_analysis/figs and csv files' 
     wd <- paste0(getwd(), mypath)
     setwd(wd)

     plotsave <- TRUE  
``` 
3. **Required Custom-written Functions**

   The custom functions are required to make the graphs etc.

```R
	custom_boxplot <- function(data, wid=1, cap=0.5, xlab = 'membrane potential (mV)', 
	                           ylab = 'PSP amplitude (mV)', xrange=c(-70,-50), yrange=c(-10,15), 
	                           lwd=0.8, type=6) {
		x <- data$x
		y <- data$y
		unique_x <- unique(data$x)
		xrange <- xrange + c(-wid, wid)
		plot(1, type="n", ylim=yrange, xlim=xrange, xlab=xlab, ylab=ylab, xaxt="n", yaxt="n", bty='n', lwd=lwd)
		
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
		
		rect(current_x - wid, q1, current_x + wid, q3, col="white", lwd=lwd)
		segments(current_x, q1, current_x, min_val, lwd=lwd)
		segments(current_x, q3, current_x, max_val, lwd=lwd)
		segments(current_x - cap, min_val, current_x + cap, min_val, lwd=lwd)
		segments(current_x - cap, max_val, current_x + cap, max_val, lwd=lwd)
		segments(current_x - wid*1.1, median_val, current_x + wid*1.1, median_val, col="black", lwd=3*lwd)
		}
		axis(1, at=unique_x, labels=unique_x)
		axis(2)
	}
		    
	# if random effects model is singular
	isSingular.fun <- function(formula, data){
		mod <- suppressMessages(lmer(formula=formula, data=data))
		isSingular(mod)
		} 
		
	fun.plot = function(data, wid=1, cap=0.5, xlab = 'membrane potential (mV)', ylab = 'PSP amplitude (mV)', xrange=c(-70,-50), yrange=c(-10,15), lwd=0.8, amount=0.5, p.cex=0.25, type=6){
		
		# Fit the model using lmer
		# model_lmer <- lmer(y ~ x + (1|s))
		x <- data$x
		y <- data$y
		s <- data$s
		formula  <- y ~ x + (1|s)
		mixed <- !isSingular.fun(formula, data)
		mod <-  if (mixed) lmer(y ~ x + (1|s)) else lm(y ~ x)
		formula <- if (!mixed) y ~ x else formula
		cat("model is ", format(formula), "\n")
		print(summary(mod))
	
		if (mixed){
			fixed_effects <- fixef(mod)
			m <- fixed_effects[[2]]
			c <- fixed_effects[[1]]
		}else{
			coeffs <- coef(mod)
			m <- coeffs[[2]]
			c <- coeffs[[1]]
		}
		
		# x intercept
		cat("x intercept is ", format(-c/m), "mV", "\n")
	    
		# use MuMIn to evaluate r2 
		r2_values <- r.squaredGLMM(mod)
		# print(r2_values)
	        cat("rsqr (marginal) ", format(r2_values[1]), " rsqr (conditional) ", format(r2_values[2]))
		# cat("rsqr (conditional) ", format(r2_values[2]))
	
		# Add a column to data for jittered x-values
		set.seed(42)
		data$x_jitter <- jitter(data$x, amount=amount)
		
		# Boxplot
		custom_boxplot(data, wid=wid, cap=cap, xlab=xlab, ylab=ylab, xrange=xrange, yrange=yrange, lwd=lwd, type=type)
		# Plot individual data points with reduced jitter, reduced size, and unfilled circles without x and y axes
		points(data$x_jitter, data$y, pch=19, bg="transparent", col="darkgray", lwd=lwd/2, cex=p.cex)
		# Connect data points within subjects with gray dotted lines
		line=TRUE
		if (line){
			subjects <- unique(data$s)
			for(subj in subjects){
	  			subset_data <- data[data$s == subj, ]
	  			lines(subset_data$x_jitter, subset_data$y, col="darkgray", lwd=lwd, lty=3)  # lty=2 for dotted line
			}
		}
		
		# Predict y values using the model for new data with matching grouping factor levels
		y_pred <- m * unique(x) + c
		# Add the line of best fit
		lines(unique(x), y_pred, col="black", lwd=lwd, lty=1)
		
		# list(reversal=-c/m, r2_values =r2_values)	
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

	# Perform the Wilcoxon Signed-Rank Test
	fun.wilcox <- function(data, paired=TRUE, alternative='two.sided', exact=NULL){
		x <- data[data$x == 2,][order(data[data$x == 2,]$s),]$y
		y <- data[data$x == 3,][order(data[data$x == 3,]$s),]$y
		wilcox.test(x, y, paired = paired, alternative=alternative, exact=exact)
	}

	# Perform the Wilcoxon Signed-Rank Test
	fun.wilcox2 <- function(data, paired=TRUE, alternative='two.sided', exact=NULL){
		x <- data[data$x == 1,][order(data[data$x == 1,]$s),]$y
		y <- data[data$x == 2,][order(data[data$x == 2,]$s),]$y
		wilcox.test(x, y, paired = paired, alternative=alternative, exact=exact)
	}
```

4. **Data Analysis**

   Having run all the previous code, this code performs the analysis and makes the graphs in the ms

```R
	# data for fig 2
	# import data sets from stored csv
	
	data2F <- import.fun('data2F')
	data1 <- data2F[, c("s", "x", "y1")]; colnames(data1)[3] <- "y"
	data2 <- data2F[, c("s", "x", "y2")]; colnames(data2)[3] <- "y"
	data3 <- data2F[, c("s", "x", "y3")]; colnames(data3)[3] <- "y"
	
	data2C <- import.fun('data2C')
	data4 <- data2C[, c("s", "x", "y1")]; colnames(data4)[3] <- "y"
	data5 <- data2C[, c("s", "x", "y2")]; colnames(data5)[3] <- "y"
	
	data6 <- import.fun('data2G')
	data7 <- import.fun('data2J')

	# FIG2C
	dev.new(width=6 ,height=4,noRStudioGD=TRUE)
	par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
	fun.plot(data4, ylab='PSC amplitude (pA)', yrange=c(-20,25), p.cex=0.6)

	fun.plot(data5, yrange=c(-20,25), p.cex=0.6)

	# FIG2F
	dev.new(width=9 ,height=4,noRStudioGD=TRUE)
	par(mar=c(1, 1, 1, 1), mfrow=c(1,3), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
	fun.plot(data1, p.cex=0.6)

	fun.plot(data2, p.cex=0.6)

	fun.plot(data3, p.cex=0.6)

	# FIG2GJ
	dev.new(width=6 ,height=4,noRStudioGD=TRUE)
	par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
	fun.plot(data6, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(-70, -55), amount=0.05, p.cex=0.6)

	fun.plot(data7, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(-70, -55), amount=0.05, p.cex=0.6)

	# FIG2GJ: statistcal tests
	fun.wilcox2(data6)
	fun.wilcox2(data7)

```
	
```R
	# Fig3C
	data8 <- import.fun('data3CA')
	data9 <- import.fun('data3CB')

	dev.new(width=6 ,height=4,noRStudioGD=TRUE)
	par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
	fun.plot(data8, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6)

	fun.plot(data9, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6)

	# Fig3C: statistcal tests
	fun.wilcox2(data8)
	fun.wilcox2(data9)

	
	# Fig3F
	data3F <- import.fun('data3F')
	data10 <- data3F[, c("s", "x", "y1")]; colnames(data10)[3] <- "y"
	data11 <- data3F[, c("s", "x", "y2")]; colnames(data11)[3] <- "y"

	dev.new(width=6 ,height=4,noRStudioGD=TRUE)
	par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
	fun.plot(data10, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6)
	fun.plot(data11, wid=0.25, cap=0.125, xrange=c(0.5, 2.5), yrange=c(0, 5), amount=0.05, p.cex=0.6)
	
	
	# Fig3F: statistcal tests
	fun.wilcox2(data10, paired = FALSE)	
	# nb result is in fact identical (pairs go in identical directions; performing Mann U as not enough values for paired (n = 5 pairs is required)
	fun.wilcox2(data11, paired = FALSE)
	
```

 
```R
 	# Fig4EF 
	data4E <- import.fun('data12')
	data12 <- data4E[, c("s", "x", "y1")]; colnames(data12)[3] <- "y"
	data13 <- data4E[, c("s", "x", "y2")]; colnames(data13)[3] <- "y"

	dev.new(width=9 ,height=4,noRStudioGD=TRUE)
	par(mar=c(1, 1, 1, 1), mfrow=c(1,2), oma = c(2, 2, 2, 0), ps=10, cex = 0.9, cex.main = 0.9)
	fun.plot(data12, yrange=c(0,35), xrange=c(0.5,3.5), xlab='', ylab='', wid=0.2, cap=0.1, amount=0, p.cex=0.6)
	fun.plot(subset(data13, x != 1), yrange=c(0,0.25), xrange=c(0.5,3.5), xlab='', ylab='', wid=0.2, cap=0.1, amount=0, p.cex=0.6)
	
	
	# stats for Fig4E (dataset combined)
	fun.wilcox(data12)
	fun.wilcox(data13)
```

```R	
	# data for figS1 
	dataS1 <- read.csv('data14.csv')
	colnames(dataS1) = c('A+B', 'C')
	
	# plot figures
	



```







## Functions

### custom_boxplot

The function, `custom_boxplot` is a customised function to create boxplots. 

The boxplot function creates a box-and-whisker plot, which is a standardized way of displaying the distribution of data based on a five-number summary:

1. Minimum: The smallest data point, excluding any outliers.
2. First quartile (Q1): The data point below which 25% of the data fall.
3. Median (Q2 or second quartile): The data point that divides the data into two halves. 50% of the data fall below the median, and 50% of the data fall above it.
4. Third quartile (Q3): The data point below which 75% of the data fall.
5. Maximum: The largest data point, excluding any outliers.
6. Outliers are defined as values more extreme than Q1 - 1.5 * iqr  for lower and Q3 + 1.5 * iqr for upper bound limits.

Any outliers are removed and the default setting for calculating the quartiles is type = 6.

In R's `quantile` function, there are 9 types of quantile algorithms, named type 1 to type 9. These methods are defined to give different treatments for the lower and upper tails and whether they should be exclusive or inclusive. 

Type 1: Inverse of empirical distribution function.

Type 2: Similar to type 1 but with averaging at discontinuities.

Type 3: SAS definition: nearest even order statistic.

Type 4: Linear interpolation of the empirical cdf.

Type 5: Piecewise linear function where the knots are the values of order statistics.

Type 6: Linear interpolation of the expectations based on the order statistics (default).

Type 7: Linear interpolation of the modes based on the order statistics.

Type 8: Linear interpolation between the points that capture the α percent and 1-α percent of the data.

Type 9: Linear interpolation of the approximate medians for order statistics.   


NOTE The R function, `boxplot` is not used to make whisker-and-box plots because this is not the method used by most graphics software. The native R function calls `boxplot.stats` which, in turn, calls `stats::fivenum` to calculate the medium iqr and min and max based on Tukey's five-number summary definition.

John Tukey's 'hinges' which are used in his five-number summary and for drawing boxplots, are similar to quartiles but can be calculated in a way that's slightly different from any of the standard quantile methods in R. Tukey's original definition involved using the median to split the data set and then finding the median of the lower and upper halves. If the data set or data half contains an odd number of points, the median is included in both halves.

How Tukey's hinges are usually computed:

1. The lower hinge is the median of the lower half of the data set (not including the overall median if the number of data points is odd).

2. The upper hinge is the median of the upper half of the data set (again not including the overall median if the number of data points is odd).

This method is somewhat akin to R's 'Type 1' method for calculating quantiles, also known as the "inverted empirical distribution function." 

The `quantile` function in R with the option `type=1` in R's default `boxplot`, the applied method is close to, but not exactly the same as, Tukey's original 'hinge' method. 

`GraphPad Prism` default seems to calculate quartiles using the method that is commonly taught, which corresponds to 'Type 7' in R's `quantile` function

### Summary

**The native R function, `boxplot` calculates whisker-and-box plots based on Tukey's original 'hinges' method by calling `stats::fivenum`**

**The function `custom_boxplot` calculates quartiles using R function `quantile`. This can be set to type = 1 to 9**

**The default in `custom_boxplot` is type = 6 which should produce similar results to `GraphPad Prism`; for results closer (but not identical) to Tukey's 'hinges' method / R's native `boxplot`, set type = 1**

**Linear Regression** is performed using the package `lmer`. The function determines whether the fit of the model is singular. If not then it fits by random mixed effects model. In `lmer` terminology, the formula for this is **y ~ x + (1|s)**. This formula specifies how the dependent variable 'y' is modeled in relation to the independent (or fixed-effect) predictor variable 'x' and the random effect of the subject 's'. 

### Random Mixed Effects Model
**y ~ x + (1|s)**: the formula specifies how the dependent variable y is modeled in relation to the predictor variable x and the random effect of the subject s. 

1. y: This is the dependent variable you are trying to model or predict.
2. ~: The tilde separates the dependent variable from the independent variables and random effects.
3. x: This is the independent (or fixed-effect) variable. The model will estimate how y varies with x.
4. +: The plus sign indicates that you are including more terms in the model.
5. (1|s): This is a random intercept for subject s. In other words, each subject is allowed to have its own baseline value of y that is randomly distributed around the overall mean of y.

In the context of linear mixed models, $R^2$ can be a bit more complex to define and interpret than in standard linear regression. There are actually two commonly reported $R^2$ values for linear mixed models:
	
1. Marginal $R^2$: Represents the variance explained by the fixed effects alone.
2. Conditional $R^2$: Represents the variance explained by both the fixed and random effects.

The conditional $R^2$ is always equal to or larger than the marginal $R^2$ since it also includes the variance explained by the random effects.

### Linear Model
**y ~ x**: the formula specifies how the dependent variable y is modeled in relation to the predictor variable x. 

Fits to the random mixed effects model may be singular if:

1. There is some redundancy in random effects: the model is too complex for the data with random effects that do not contribute much variance. For example, the random intercepts for each level of subject, s may not vary enough from each other. This can happen if there are too few levels in s or if the data within each level of s does not vary much.

2. The variance of the random effects or the residual variance is estimated to be near zero. This suggests that the random effect might not be necessary, as it doesn't explain a significant amount of variability in the data.

If fits of **y ~ x + (1|s)** are singular, the function simplifies the model to linear regression **y ~ x**.

The returned values for conditional $R^2$ and marginal $R^2$ are the same in a linear model since the two definitions effectively converge to a common value. 

The linear model only includes fixed effects and there are no random effects to consider. Therefore, the variance explained by the fixed effects (marginal $R^2$ ) is the total variance explained by the model (and is usually just referred to as $R^2$).

     






 
Note: the use of X11 (including tcltk) requires XQuartz (version 2.8.5 or later). 

always re-install XQuartz when upgrading your macOS to a new major version. 




n the boxplot.default function in R, the calculation of the hinges (which correspond to the first and third quartiles of the data) is not done directly within this function. Instead, the boxplot.default function calls another function, boxplot.stats, for each group of data it is plotting. This boxplot.stats function is where the actual computation of the hinges and other statistics necessary for the boxplot takes place.

The boxplot.stats function calculates the hinges based on the quartiles of the data. The quartiles are typically calculated using a method that is similar to the Type 7 quantile algorithm in R, which is a part of the quantile function. This method involves interpolating between data points to compute quartiles, which is a common approach in statistical analysis.








