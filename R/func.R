
#' @title 	weighted_x_sd
#' @descrip Calculates weighted mean and standard deviation across two groups 
#' @param x1 Mean of group 1
#' @param sd1 Standard deviation of group 1
#' @param n1 Sample size of group 1
#' @param m2 Mean of group 2
#' @param sd2 Standard deviation of group 2
#' @param n2 Sample size of group 2

weighted_x_sd <- function(x1, sd1, n1, x2, sd2, n2){
	     mi <- c(x1, x2)
	   	 ni <- c(n1, n2)
	   	sdi <- c(sd1, sd2)

		# Calculate weighted mean
		m.weighted <- sum(ni*mi)/sum(ni)

		#Calculate weighted SD accounting for between group differences in mean. Formula assumes normality! See simulations
		sd.weighted <- sqrt((sum((ni-1) * sdi^2) + sum(ni*(mi - m.weighted)^2)) / (sum(ni) - 1))

		#tot.n <- sum(ni)

		weighted <- c(m.weighted, sd.weighted)

	return(weighted)
}

#' @title calc_weighted
#' @descrip Calculates weighted mean and standard deviation across two groups across a data set
#' @param data Dataset that contains the means (x), standard deviations (sd), and sample size (n) for the two groups need

calc_weighted <- function(data){
	tmp <- as.data.frame(t(apply(data, 1, function(x) 
									weighted_x_sd( x1  = x[1], 
										   		   sd1 = x[2], 
										   		   n1  = x[3], 
										   		   x2  = x[4], 
										   		   sd2 = x[5], 
										   		   n2  = x[6]))
									)
							)
	return(tmp)
}


#' @title 	calc_weighted_factors
#' @descrip Calculates weighted mean across four groups / factors
#' @param m1_1 Character identify the column for the mean of group 1,1
#' @param m1_2 Character identify the column mean of group 1,2
#' @param m2_1 Character identify the column mean of group 2,1
#' @param m2_2 Character identify the column mean of group 2,2
#' @param sd1_1 Character identify the column identifying the standard deviation of group 1,1
#' @param sd1_2 Character identify the column identifying the standard deviation of group 1,2
#' @param sd2_1 Character identify the column identifying the standard deviation of group 2,1
#' @param sd2_2 Character identify the column identifying the standard deviation of group 2,2
#' @param n1_1 Character identify the column identifying the sample size of group 1,1
#' @param n1_2 Character identify the column identifying the sample size of group 1,2
#' @param n2_1 Character identify the column identifying the sample size of group 2,1
#' @param n2_2 Character identify the column identifying the sample size of group 2,2
#' @example
#' \donttest{
#' data <- data.frame(m1_1 = c(15, 20), sd1_1 = c(2, 3), n1_1 = c(20, 30), m1_2 = c(15, 15), sd1_2 = c(2, 6), n1_2 = c(20, 30), m2_1 = c(50, 70), sd2_1 = c(10, 8), n2_1 = c(40, 30), m2_2 = c(25, 30) , sd2_2 = c(8, 6), n2_2 = c(80, 30))
#' 
#' calc_weighted_factors(m1_1 = "m1_1", sd1_1 = "sd1_1", n1_1 = "n1_1", m1_2 = "m1_2", sd1_2 = "sd1_2", n1_2 = "n1_2", m2_1 = "m2_1", sd2_1 = "sd2_1", n2_1 = "n2_1", m2_2 = "m2_2", sd2_2 = "sd2_2", n2_2 = "n2_2", data = data)
#' 
#' calc_weighted_factors(m1_1 = "m1_1", sd1_1 = "sd1_1", n1_1 = "n1_1", m1_2 = "m1_2", sd1_2 = "sd1_2", n1_2 = "n1_2", m2_1 = "m2_1", sd2_1 = "sd2_1", n2_1 = "n2_1", m2_2 = "m2_2", sd2_2 = "sd2_2", n2_2 = "n2_2", data = data, append = FALSE)
#' }
#' @export

calc_weighted_factors <- function(m1_1, sd1_1, n1_1, m1_2, sd1_2, n1_2, m2_1, sd2_1, n2_1, m2_2, sd2_2, n2_2, data, append = TRUE) {

		data_grp1 <- data[,c(m1_1, sd1_1, n1_1, m1_2, sd1_2, n1_2)]
		data_grp2 <- data[,c(m2_1, sd2_1, n2_1, m2_2, sd2_2, n2_2)]
		data_grp3 <- data[,c(m1_1, sd1_1, n1_1, m2_1, sd2_1, n2_1)]
		data_grp4 <- data[,c(m1_2, sd1_2, n1_2, m2_2, sd2_2, n2_2)]

		weighted_grp1 <- calc_weighted(data_grp1)
		colnames(weighted_grp1) <- c("m1*", "sd1*")

		weighted_grp2 <- calc_weighted(data_grp2)
		colnames(weighted_grp2) <- c("m2*", "sd2*")

		weighted_grp3 <- calc_weighted(data_grp3)
		colnames(weighted_grp3) <- c("m*1", "sd*1")

		weighted_grp4 <- calc_weighted(data_grp4)
		colnames(weighted_grp4) <- c("m*2", "sd*2")

		weighted_avg_sd <- cbind(weighted_grp1, weighted_grp2, weighted_grp3, weighted_grp4)

		if(append){
			data <- cbind(data, weighted_avg_sd)
			return(data)
		} else{
			weighted_avg_sd
		}	
}

#' @title 
#' @descrip Hedges' SMD interaction effect size
#' @param m1_1 Character identify the column for the mean of group 1,1
#' @param m1_2 Character identify the column mean of group 1,2
#' @param m2_1 Character identify the column mean of group 2,1
#' @param m2_2 Character identify the column mean of group 2,2
#' @param sd1_1 Character identify the column identifying the standard deviation of group 1,1
#' @param sd1_2 Character identify the column identifying the standard deviation of group 1,2
#' @param sd2_1 Character identify the column identifying the standard deviation of group 2,1
#' @param sd2_2 Character identify the column identifying the standard deviation of group 2,2
#' @param n1_1 Character identify the column identifying the sample size of group 1,1
#' @param n1_2 Character identify the column identifying the sample size of group 1,2
#' @param n2_1 Character identify the column identifying the sample size of group 2,1
#' @param n2_2 Character identify the column identifying the sample size of group 2,2


hedge_d_f1.f2 <- function(m1_1, sd1_1, n1_1, m1_2, sd1_2, n1_2, m2_1, sd2_1, n2_1, m2_2, sd2_2, n2_2, type = c("main", "interaction"), append = TRUE, data){
	library(magrittr)
	type = match.arg(type)

	if(type == "main"){
	
	hedge_SDpool <- sqrt(   ( ((data[,n1_1] - 1)*(data[,sd1_1]^2)) +  ((data[,n1_2] - 1)*(data[,sd1_2]^2)) + ((data[,n2_1] - 1)*(data[,sd2_1]^2)) + ((data[,n2_2] - 1)*(data[,sd2_2]^2))) / (data[,n1_1] + data[,n1_2] + data[,n2_1] + data[,n2_2] - 4) )
	  
	     hedge_d_f1 <- ((data[,m1_1] + data[,m1_2]) - (data[,m2_1] + data[,m2_2])) / (2*hedge_SDpool)
	         V_d_f1 <- (1/4) * ( (1/data[,n1_1]) + (1/data[,n1_2]) + (1/data[,n2_1]) + (1/data[,n2_2]) + (hedge_d_f1^2) / (2*(data[,n1_1] + data[,n1_2] + data[,n2_1] + data[,n2_2])))

	     hedge_d_f2 <- ((data[,m1_1] + data[,m2_1]) - (data[,m1_2] + data[,m2_2])) / (2*hedge_SDpool)
	         V_d_f2 <- (1/4) * ( (1/data[,n1_1]) + (1/data[,n1_2]) + (1/data[,n2_1]) + (1/data[,n2_2]) + (hedge_d_f2^2) / (2*(data[,n1_1] + data[,n1_2] + data[,n2_1] + data[,n2_2])))

	hedge_main <- round_df(data.frame(hedge_d_f1 = hedge_d_f1, V_d_f1 = V_d_f1, hedge_d_f2 = hedge_d_f2, V_d_f2 = V_d_f2), digits = 3)

	if(append){
		return(cbind(data, hedge_main))
	} else{
		return(hedge_main)
	}

	}


	if(type == "interaction"){

	}
}

#' @title 
#' @descrip lnRR main effects for each factor and interaction
#' @param m1_1 Character identify the column for the mean of group 1,1
#' @param m1_2 Character identify the column mean of group 1,2
#' @param m2_1 Character identify the column mean of group 2,1
#' @param m2_2 Character identify the column mean of group 2,2
#' @param sd1_1 Character identify the column identifying the standard deviation of group 1,1
#' @param sd1_2 Character identify the column identifying the standard deviation of group 1,2
#' @param sd2_1 Character identify the column identifying the standard deviation of group 2,1
#' @param sd2_2 Character identify the column identifying the standard deviation of group 2,2
#' @param n1_1 Character identify the column identifying the sample size of group 1,1
#' @param n1_2 Character identify the column identifying the sample size of group 1,2
#' @param n2_1 Character identify the column identifying the sample size of group 2,1
#' @param n2_2 Character identify the column identifying the sample size of group 2,2
#' @param m1*  Character identify the column identifying the weighted mean of group 1,*
#' @param m2*  Character identify the column identifying the weighted mean of group 2,*
#' @param sd1*  Character identify the column identifying the weighted standard deviation of group 1,*
#' @param sd2*  Character identify the column identifying the weighted standard deviation of group 2,*
#' @param n1*  Character identify the column identifying the total sample size of combined of group 1,*
#' @param n2*  Character identify the column identifying the total sample size of combined of group 2,*


lnRR_f1 <- function(m1_1, sd1_1, n1_1, m1_2, sd1_2, n1_2, m2_1, sd2_1, n2_1, m2_2, sd2_2, n2_2, m1, m2, sd1, sd2, data, type = "nakaMorris"){

	#if(type == "morris"){
	lnRRf1_morris <- log(( data[,m1_1]  + data[,m1_2] ) / ( data[,m2_1]  + data[,m2_2] ) )

	lnRRf1_naka <- 0.5 * log((data[,m1_1]*data[,m1_2]) / (data[,m2_1]*data[,m2_2]))

	lnRRf1_V_morris <- (( 1 / (data[,m1_1]  + data[,m1_2] ))^2 * (((data[,sd1_1])^2 / data[,n1_1]) + ((data[,sd1_2])^2 / data[,n1_2]))) + (( 1 / (data[,m2_1]  + data[,m2_2] ))^2 * (((data[,sd2_1])^2 / data[,n2_1]) + ((data[,sd2_2])^2 / data[,n2_2])))

	lnRRf1_naka_V <- (1/4) * (((data[,sd1_1])^2 / ((data[,m1_1])^2*data[,n1_1])) + ((data[,sd1_2])^2 / ((data[,m1_2])^2*data[,n1_2])) + ((data[,sd2_1])^2 / ((data[,m2_1])^2*data[,n2_2])) + ((data[,sd2_2])^2 / ((data[,m2_2])^2*data[,n2_2])))
 	#lnRRf2 <- log(( data[,m1_1]  + data[,m2_1] ) / 2 ) - log(( data[,m1_2]  + data[,m2_2] ) / 2 )
 	#lnRRf1f2 <- 

	#}

	#if(type == "naka"){
		  

		

	#}

		data$n1_w <- data[,n1_1] + data[,n1_2]
		data$n2_w <- data[,n2_1] + data[,n2_2]
		lnRR_weight <- metafor::escalc(measure = "ROM", m1i= data[,m1], m2i = data[,m2], sd1i= data[,sd1], sd2i= data[,sd2], n1i = data[,"n1_w"], n2i=data[,"n2_w"], append = FALSE)

	lnRR_data <- data.frame(lnRRf1_morris = lnRRf1_morris, lnRRf1_V_morris = lnRRf1_V_morris, lnRRf1_naka = lnRRf1_naka, lnRRf1_naka_V = lnRRf1_naka_V)

	if(type == "nakaMorris"){
	return(lnRR_data)
	}else{
		return(lnRR_weight)
	}
}

#' @title round_df
#' @descrip Rounding numeric columns in data frames
#' @param data Data frame wanting to be rounded
#' @param digits Number of digits numeric values are to be rounded to

round_df <- function(x, digits) {
    # round all numeric variables
    # x: data frame 
    # digits: number of digits to round
    numeric_columns <- sapply(x, mode) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x
}

## Here we want to test out and compare between three approaches, the first generating normal lnRR with weighted mean and SD, second, using lnRR from Morris with the arithmetic mean and third Nakagawa eta l's revision using the geometric mean. We should also test a geometric mean version as well. 
data <- data.frame(m1_1 = c(100, 20), sd1_1 = c(30, 30), n1_1 = c(500, 500), m1_2 = c(109, 15), sd1_2 = c(30, 30), n1_2 = c(500, 500), m2_1 = c(115, 70), sd2_1 = c(30, 30), n2_1 = c(500, 500), m2_2 = c(124, 30) , sd2_2 = c(30, 30), n2_2 = c(500, 500))

dat_weights <- calc_weighted_factors(m1_1 = "m1_1", sd1_1 = "sd1_1", n1_1 = "n1_1", m1_2 = "m1_2", sd1_2 = "sd1_2", n1_2 = "n1_2", m2_1 = "m2_1", sd2_1 = "sd2_1", n2_1 = "n2_1", m2_2 = "m2_2", sd2_2 = "sd2_2", n2_2 = "n2_2", data = data)

# Takes lnRR with the weighted mean values
lnRR_f1(m1_1 = "m1_1", sd1_1 = "sd1_1", n1_1 = "n1_1", m1_2 = "m1_2", sd1_2 = "sd1_2", n1_2 = "n1_2", m2_1 = "m2_1", sd2_1 = "sd2_1", n2_1 = "n2_1", m2_2 = "m2_2", sd2_2 = "sd2_2", n2_2 = "n2_2", m1= "m1*", m2= "m2*", sd1= "sd1*", sd2= "sd1*", data = dat_weights, type = "weighted")

# Takes lnRR with the weighted mean values
lnRR_f1(m1_1 = "m1_1", sd1_1 = "sd1_1", n1_1 = "n1_1", m1_2 = "m1_2", sd1_2 = "sd1_2", n1_2 = "n1_2", m2_1 = "m2_1", sd2_1 = "sd2_1", n2_1 = "n2_1", m2_2 = "m2_2", sd2_2 = "sd2_2", n2_2 = "n2_2", m1= "m1*", m2= "m2*", sd1= "sd1*", sd2= "sd1*", data = dat_weights, type = "nakaMorris")


### Here we are interested in testing whether, taking the marginalised weighted means and SDs impact the SMD estimates and how they compare with Hedges' formulation for factor 1 equation. Reason is that SD weighted contains both within and between SD, whereas Hedges just used within and this can affect things:

data <- data.frame(m1_1 = c(100, 20), sd1_1 = c(30, 30), n1_1 = c(500, 500), m1_2 = c(109, 15), sd1_2 = c(30, 30), n1_2 = c(500, 500), m2_1 = c(115, 70), sd2_1 = c(30, 30), n2_1 = c(500, 500), m2_2 = c(124, 30) , sd2_2 = c(30, 30), n2_2 = c(500, 500))

# calculate hedges main effect for f1 and f2, but we focus only on f1 here. 
hedge_d_f1.f2(m1_1 = "m1_1", sd1_1 = "sd1_1", n1_1 = "n1_1", m1_2 = "m1_2", sd1_2 = "sd1_2", n1_2 = "n1_2", m2_1 = "m2_1", sd2_1 = "sd2_1", n2_1 = "n2_1", m2_2 = "m2_2", sd2_2 = "sd2_2", n2_2 = "n2_2", data = data, type = "main")

# calculate the weighted mean and SD
dat_weights <- calc_weighted_factors(m1_1 = "m1_1", sd1_1 = "sd1_1", n1_1 = "n1_1", m1_2 = "m1_2", sd1_2 = "sd1_2", n1_2 = "n1_2", m2_1 = "m2_1", sd2_1 = "sd2_1", n2_1 = "n2_1", m2_2 = "m2_2", sd2_2 = "sd2_2", n2_2 = "n2_2", data = data)

        dat_weights$n1_w <- data[,n1_1] + data[,n1_2]
		dat_weights$n2_w <- data[,n2_1] + data[,n2_2]
		SMD_weight <- metafor::escalc(measure = "SMD", m1i= dat_weights[,"m1*"], m2i = dat_weights[,"m2*"], sd1i= dat_weights[,"sd1*"], sd2i= dat_weights[,"sd2*"], n1i = dat_weights[,"n1_w"], n2i=dat_weights[,"n2_w"], append = FALSE)

## Example code
weighted_x_sd(100, 30, 500, 109, 30, 500)
weighted_x_sd(115, 30, 500, 124, 30, 500)

weighted_x_sd(100, 30, 500, 115, 30, 500)
weighted_x_sd(109, 30, 500, 124, 30, 500)