### R analysis for Day et al., 2023

Recreates graphical outputs and statistical analyses



Note: the use of X11 (including tcltk) requires XQuartz (version 2.8.5 or later). 

always re-install XQuartz when upgrading your macOS to a new major version. 


## Data Analysis

this code should pull data from csv files and analyse it/produce boxplots etc and recreate relevant parts of the figures
returns svg files that are then used to produce figures in Adobe Illustrator

The final analysis and figures presented in the manuscript were generated using R. 

The analyses were conducted in the R graphical user interface (GUI):
  - R version 4.3.1 – "Beagle Scouts"
  - [R Statistical Software](https://www.R-project.org/)

  1. **Load Packages**

     The required packages are [`MuMIn`](https://cran.r-project.org/web/packages/MuMIn/index.html), [`svglite`](https://cran.r-project.org/web/packages/svglite/index.html) and [`lme4`](https://www.rdocumentation.org/packages/nlme/versions/3.1-163/topics/lme) 

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
3. **Required Custom Functions**

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
     ```
     

customised function to create boxplots

In R's quantile() function, there are 9 types of quantile algorithms, named type 1 to type 9. 

These methods are defined to give different treatments for the lower and upper tails and whether they should be exclusive or inclusive. 

Type 1: Inverse of empirical distribution function.

Type 2: Similar to type 1 but with averaging at discontinuities.

Type 3: SAS definition: nearest even order statistic.

Type 4: Linear interpolation of the empirical cdf.

Type 5: Piecewise linear function where the knots are the values of order statistics.

Type 6: Linear interpolation of the expectations based on the order statistics (default).

Type 7: Linear interpolation of the modes based on the order statistics.

Type 8: Linear interpolation between the points that capture the α percent and 1-α percent of the data.

Type 9: Linear interpolation of the approximate medians for order statistics.   

quantile.inc in excel is 7 

NOTE not using R function boxplot to make whisker and box because this is NOT the method used by most graphics software 

the r function boxplot calls stats::fivenum to calculate the medium iqr and min and max; based on Tukey's five-number summary definition

John Tukey's "hinges," which are used in his five-number summary and for drawing boxplots, are similar to quartiles but can be calculated in a way that's slightly different from any of the standard quantile methods in R. Tukey's original definition involved using the median to split the data set and then finding the median of the lower and upper halves. If the data set or data half contains an odd number of points, the median is included in both halves.

How Tukey's hinges are usually computed:

1. The lower hinge is the median of the lower half of the data set (not including the overall median if the number of data points is odd).

2. The upper hinge is the median of the upper half of the data set (again not including the overall median if the number of data points is odd).

This method is somewhat akin to R's Type 1 method for calculating quantiles, also known as the "inverted empirical distribution function." 

The `quantile` function in R with the option `type=1` in R's default boxplot(), the applied method is close to, but not exactly the same as, Tukey's original method for hinges. 

GraphPad Prism generally calculates quartiles using the method that is commonly taught, which corresponds to "Type 7" in R's quantile() function

R default and also here is type=6; should produce similar results to GraphPad




#     notes on random mixed effects model
#     y ~ x + (1|s):
#     formula specifies how the dependent variable y is modeled in relation to the predictor variable x and the random effect of the subject s 
#     y: This is the dependent variable you are trying to model or predict.
#     ~: The tilde separates the dependent variable from the independent variables and random effects.
#     x: This is the independent (or fixed-effect) variable. The model will estimate how y varies with x.
#     +: The plus sign indicates that you are including more terms in the model.
#     (1|s): This is a random intercept for subject s. In other words, each subject is allowed to have its own baseline value of y that is randomly distributed around the overall mean of y.

    
    
# if random effects model is singular
isSingular.fun <- function(formula, data){
	mod <- suppressMessages(lmer(formula=formula, data=data))
	isSingular(mod)
	} 
	
fun.plot = function(data, wid=1, cap=0.5, xlab = 'membrane potential (mV)', ylab = 'PSP amplitude (mV)', xrange=c(-70,-50), yrange=c(-10,15), lwd=0.8, amount=0.5, p.cex=0.25, type=6){
	
	# Fit the model using lme4
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

	
	# In the context of linear mixed models, Rsqr can be a bit more complex to define and interpret than in standard linear regression. 
	# There are actually two commonly reported Rsqr values for linear mixed models:

	#     Marginal R2: Represents the variance explained by the fixed effects alone.	
	#     Conditional R2: Represents the variance explained by both the fixed and random effects.
    
	# use MuMIn to evaluate r2 if necessary
	r2_values <- r.squaredGLMM(mod)
	# print(r2_values)
    cat("rsqr (marginal) ", format(r2_values[1]), " rsqr (conditional) ", format(r2_values[2]))
	# cat("rsqr (conditional) ", format(r2_values[2]))
	# 	 The boxplot function creates a box-and-whisker plot, which is a standardized way of displaying the distribution of data based on a five-number summary:
	#    Minimum: The smallest data point, including any outliers.
	#    First quartile (Q1): The data point below which 25% of the data fall.
	#    Median (Q2 or second quartile): The data point that divides the data into two halves. 50% of the data fall below the median, and 50% of the data fall above it.
	#    Third quartile (Q3): The data point below which 75% of the data fall.
	#    Maximum: The largest data point, including any outliers.

	# Add a column to data for jittered x-values
	set.seed(42)
	data$x_jitter <- jitter(data$x, amount=amount)
	
	# Boxplots for each x value with outliers shown
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

# simple import function if NA is zero imports all exlcude excludes those subjects s
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
