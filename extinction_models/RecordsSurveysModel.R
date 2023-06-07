##########################################################################################
## CALCULATING THE PROBABILITY THAT A SPECIES IS EXTINCT, P(E)
## Modified from the R code provided in the Supplementary Data for:
## Inferring extinctions II: a practical, iterative model based on records and surveys.
## Biological Conservation 214:328-335 (2017) by
## Colin J. Thompson, Vira Koshkina, Mark A. Burgman, Stuart H.M. Butchart and Lewi Stone
## Updated by H.R. Akcakaya 19.June.2018 and 30.June.2019.  
##
##            To use this code, please see Red List Guidelines, section 11
##            http://www.iucnredlist.org/documents/RedListGuidelines.pdf
##
##########################################################################################

# Step 1: Fill in the Excel file template for your species. 
#         Rename the Excel file with your species name, such as Mammuthus primigenius.xlsx

# Step 2: Enter the same name below, in double-quotes, without the file extension
#         Enter the data folder (path), using forward slashes (/)

Species     <- "Mammuthus primigenius"
DataFolder  <- "C:/RedList/Extinction/test"
PX0         <- 1 # prob. that species is extant before the first record or survey year

# Step 3: Run this code

# Step 4: To interpret the results, see the Red List Guidelines


################### NO USER-SPECIFIED DATA BELOW THIS LINE ###################

if (!require(readxl)) install.packages('readxl')
library(readxl)

# Output file names
ImgFile <- paste(Species, "-Xt.jpg", sep = "")
Img2File <- paste(Species, "-PEx.jpg", sep = "")
XtFile  <- paste(Species, "-Xt.csv", sep = "")
EXFile  <- paste(Species, "-EXTINCT.csv", sep = "")
InputFile <- paste(Species, "-input.csv", sep = "") # all the input data including passive surveys

# load data files for recordings and surveys
setwd(DataFolder)
InputXL <- paste(Species, ".xlsx", sep = "") # 

recordings <- as.data.frame(read_excel(InputXL, sheet = "Records", range = anchored("A9", dim = c(1000, 4)), col_names = TRUE))
recordings <- recordings[complete.cases(recordings), ]

surveys <- as.data.frame(read_excel(InputXL, sheet = "Surveys", range = anchored("A14", dim = c(1000, 10)), col_names = TRUE))
surveys <- surveys[complete.cases(surveys), ]

passive.surveys <- as.data.frame(read_excel(InputXL, sheet = "Surveys", range = "B11:J12" , col_names = TRUE))
eps.passive <- unlist(c(passive.surveys[,1:3]))
names(eps.passive) <- NULL
pr.passive <- unlist(c(passive.surveys[,4:6]))
names(pr.passive) <- NULL
pi.passive <- unlist(c(passive.surveys[,7:9]))
names(pi.passive) <- NULL

Threats <- as.data.frame(read_excel(InputXL, sheet = "Threats", range = "B11:D13" , col_names = TRUE))
ThreatsX <- c(Threats[1,1]*Threats[2,1], Threats[1,2]*Threats[2,2], Threats[1,3]*Threats[2,3])

                         

CurrentYear <- read_excel(InputXL, sheet = "Threats", range = "B9" , col_names = FALSE)
CurrentYear <- unlist(c(CurrentYear))
names(CurrentYear)<-NULL

## Sampling of the values, counting st dev then setting

if (!require(stats)) install.packages('stats')
library(stats)

if (!require(plyr)) install.packages('plyr')
library(plyr)

if (!require(stringr)) install.packages('stringr')
library(stringr)

# function for calculating mid-range values -------------------------
pxt.recording = function(pci, pxt) {
	pci + (1 - pci) * pxt
}
pxt.survey = function(eps, pi, pr, pxt) {
	(1 - eps * pi * pr) * pxt
}


# function that calculates the Px
px.mid = function(){
	PXt = NULL
	PXt.min = NULL
	PXt.max = NULL
	PXt.best = NULL
	sd = NULL
	# PX0 = 1

	# first year t=1
	rec = recordings[recordings[, 'year'] == years[1],]
	pci.mid = (rec[, 'pci_lower'] + rec[, 'pci_upper']) / 2
	pci.min = rec[, 'pci_lower']
	pci.max = rec[, 'pci_upper']
	pci.best = rec[, 'pci_best']

	PXt[1] = pxt.recording(pci.mid, PX0)
	PXt.min[1] = pxt.recording(pci.min, PX0)
	PXt.max[1] = pxt.recording(pci.max, PX0)
	PXt.best[1] = pxt.recording(pci.best, PX0)

	n = 10000 #number of samples
	pxt.sam = rep(PX0, n)
	stdev = (rec[, 'pci_upper'] - rec[, 'pci_lower']) / 4
	pci.sam = rnorm(n, pci.mid, stdev)
	pxt.sam = pxt.recording(pci.sam, pxt.sam)

	sd[1] = sd(pxt.sam)


	for (t in 2:length(years)) {
		#calculating rj
		if (rec.year[t, 'type'])
			#if recording
		{
			rec = recordings[recordings[, 'year'] == years[t], ]

			#get mid point estimate
			pci.mid = (rec[, 'pci_lower'] + rec[, 'pci_upper']) / 2

			pci.min = rec[, 'pci_lower']
			pci.max = rec[, 'pci_upper']
			pci.best = rec[, 'pci_best']

			PXt[t] = pxt.recording(pci.mid, PXt[t - 1])

			PXt.min[t] = pxt.recording(pci.min, PXt.min[t - 1])
			PXt.max[t] = pxt.recording(pci.max, PXt.max[t - 1])
			PXt.best[t] = pxt.recording(pci.best, PXt.best[t - 1])

			#sample to get min and max bounds
			stdev = (rec[, 'pci_upper'] - rec[, 'pci_lower']) / 4
			pci.sam = rnorm(n, pci.mid, stdev)
			pxt.sam = pxt.recording(pci.sam, pxt.sam)
			sd[t] = sd(pxt.sam)


		}  else
			#if survey
		{
			sur = surveys[surveys[, 'year'] == years[t], ]
			eps.mid = (sur[, 'eps_lower'] + sur[, 'eps_upper']) / 2
			pi.mid = (sur[, 'pi_lower'] + sur[, 'pi_upper']) / 2
			pr.mid = (sur[, 'pr_lower'] + sur[, 'pr_upper']) / 2


			eps.min = sur[, 'eps_lower']
			pi.min = sur[, 'pi_lower']
			pr.min = sur[, 'pr_lower']

			eps.max = sur[, 'eps_upper']
			pi.max = sur[, 'pi_upper']
			pr.max = sur[, 'pr_upper']

			eps.best = sur[, 'eps_best']
			pi.best = sur[, 'pi_best']
			pr.best = sur[, 'pr_best']


			PXt[t] = pxt.survey(eps.mid, pi.mid, pr.mid, PXt[t - 1])
			PXt.min[t] = pxt.survey(eps.max, pi.max, pr.max, PXt.min[t - 1])
			PXt.max[t] = pxt.survey(eps.min, pi.min, pr.min, PXt.max[t - 1])
			PXt.best[t] = pxt.survey(eps.best, pi.best, pr.best, PXt.best[t - 1])


			eps.sam = rnorm(n, eps.mid, (sur[, 'eps_upper'] - sur[, 'eps_lower']) / 4)
			pi.sam = rnorm(n, pi.mid, (sur[, 'pi_upper'] - sur[, 'pi_lower']) / 4)
			pr.sam = rnorm(n, pi.mid, (sur[, 'pr_upper'] - sur[, 'pr_lower']) / 4)
			pxt.sam = pxt.survey(eps.sam, pi.sam, pr.sam, pxt.sam)
			sd[t] = sd(pxt.sam)
		}

	}
	# print(cbind(PXt-3*sd,PXt,PXt+3*sd))
	return (cbind(PXt_lower=(PXt - 3 * sd), PXt, PXt_upper=(PXt + 3 * sd), PXt.min, PXt.max, PXt.best))

}



###### Check column names in the input files and rename them if nesessary
# print the names of the columns in records and survey data files to check the column headers and make sure the correct information is being extracted

# coumn names used in the code
rec_names = c("year", "pci_lower", "pci_best", "pci_upper")
sur_names = c(
	"year",
	"eps_lower",
	"eps_best",
	"eps_upper",
	"pr_lower",
	"pr_best",
	"pr_upper",
	"pi_lower",
	"pi_best",
	"pi_upper"
)

## Check that data files have the right number of columns
if( length(rec_names) != dim(recordings)[2] ) stop(paste('Records data must have ',length(rec_names)," columns", sep="") )
if( length(sur_names) != dim(surveys)[2] ) stop(paste('Surveys data have ',length(sur_names)," columns", sep="") )


# Explanation of the data in the columns
required_collumns_records = c("Calendar year", "pci lower", "pci best guess","pci upper")
required_collumns_surveys = c("Calendar year", "epsilon lower","epsilon best guess", "epsilon upper", "pr lower","pr best guess",  "pr upper",  "pi lower","pi best guess",  "pi upper" )



#function that adds spaces to the right so all the names aare the same length
pretty = function(str_vec){ str_pad(str_vec, max(nchar(str_vec)), pad = " ",side='right') }



#if names don't match - rename them
if (any(sort(rec_names)!=sort(colnames(recordings)))) {

	cat(paste("The following Records data are extracted from \"",InputXL ,"\" and renamed :\n", sep="" ))
	cat(paste(pretty(required_collumns_records),"	:	",pretty(colnames(recordings)), " 	=>	",rec_names) ,sep="\n")
	colnames(recordings) = rec_names


} else {
	recordings=recordings[,rec_names]
	cat(paste("\nThe following Records data are extracted from \"",InputXL ,"\" :\n", sep="" ))
	cat(paste(pretty(required_collumns_records),"	:	",pretty(colnames(recordings))),sep="\n")

}


if (any(sort(sur_names)!=sort(colnames(surveys)))) {

	cat(paste("The following Survey data are extracted from \"",InputXL ,"\" and renamed: \n ", sep="" ))
	cat(paste(pretty(required_collumns_surveys),"	:	",pretty(colnames(surveys)), " 	=>	",sur_names) ,sep="\n")
	colnames(surveys) = sur_names

} else {

	surveys=surveys[,sur_names]
	cat(paste("\nThe following Records data are extracted from \"",InputXL ,"\" :\n", sep="" ))
	cat(paste(pretty(required_collumns_surveys),"	:	",pretty(colnames(surveys))),sep="\n")

}


# check if all input data are within 0 and 1 and lower bounds < upper bounds
if( any(recordings[, -which(names(recordings) %in% c('year'))] < 0 | recordings[,-which(names(recordings) %in% c('year'))] > 1) ) stop('Records data not between 0 and 1')
if( any(recordings[,'pci_lower'] > recordings[,'pci_upper'] )) stop('Pci: lower bound > upper bound ')

if( any(surveys[,-which(names(surveys) %in% c('year'))] < 0 | surveys[,-which(names(surveys) %in% c('year'))] > 1) ) stop('Survey data not between 0 and 1')
if( any(surveys[, 'eps_lower'] > surveys[,'eps_upper'] )) stop('Survey data Epsilon: lower bound > upper bound ')
if( any(surveys[,'pr_lower'] > surveys[,'pr_upper'] )) stop('Survey data Pr: lower bound > upper bound ')
if( any(surveys[,'pi_lower'] > surveys[,'pi_lower'] )) stop('Survey data Pi: lower bound > upper bound ')


if( eps.passive[1] > eps.passive[2] ) stop('Passive surveys Epsilon: lower bound > upper bound ')
if( pr.passive[1] > pr.passive[2]) stop('Passive surveys Pr: lower bound > upper bound ')
if( pi.passive[1] > pi.passive[2]) stop('Passive surveys Pi: lower bound > upper bound ')

# check that the earliest data point is a recording
if( min(recordings[, 'year']) > min(surveys[, 'year']) ) stop('The earliest data point must be a recording')



##### Add passive surveys

# create a sequence from the earliest record year to the latest record/survey/current year
years = seq(min(recordings[, 'year']),
						max(CurrentYear, recordings[, 'year'], surveys[, 'year']),
						by = 1)
Total.years = length(years)

#find years without either surveys or recordings and create a passive survey record for each of those years
pas.sur = c(eps.passive , pi.passive, pr.passive)

pas.sur.years = years[!years %in% surveys[, 'year'] &
												!years %in% recordings[, 'year']]
pas.surveys = cbind(pas.sur.years, t(replicate(length(pas.sur.years), pas.sur)))

colnames(pas.surveys) = sur_names
surveys = rbind(surveys, pas.surveys)


# making a dataset with both records and survey

# add a column type with 1 for records and 0 for surveys
rec_type = cbind(type = 1, recordings)
sur_type = cbind(type = 0, surveys)
data = rbind.fill(rec_type, sur_type)
data = data[order(data[, 'year']), ]

# save all input data
write.csv(data, InputFile)

# check if there is both a survey and a record for any of the years
duplicates = data[duplicated(data[, 'year']), ]


if (dim(duplicates)[1] != 0) {

	cat("\nFollowing data are ignored because a record for the same year already exists:\n")

	for (i in 1:dim(duplicates)[1]) {
		if (duplicates[i, 'type'] == 0)
			print(duplicates[i, c('year',"eps_lower", "eps_upper", "pr_lower",  "pr_upper",  "pi_lower",  "pi_upper" )])

		if (duplicates[i, 1] == 1)
			print(duplicates[i, c('year', "pci_lower",  "pci_upper" )])

	}
}

#remove duplicates from the data
data <-  data[!duplicated(data[, "year"]),]

# create a matrix of years on records and surveys. First column - calendar year; second column - 1 if record, 0 if survey
rec.year = data[, c("year", "type")]


# calculate probability of species being extant
PXt = px.mid()


# Prepare a table of results -----------------------------
PXt = cbind(years, PXt)

#check if
if (ncol(PXt) != 7) { stop('Number of columns in the result table is not equal 7')
} else {

	# save Pxt to a file
	write.csv(PXt, XtFile)


	# Calculate everything for plotting the result the results -----------------------------
	xx <- c(years, rev(years))
	yysd = c(PXt[, "PXt_lower"], rev(PXt[, "PXt_upper"]))
	yyint = c(PXt[, "PXt.min"], rev(PXt[, "PXt.max"]))
	par(family = "serif")

	# save the plot to a file
	jpeg(ImgFile, width = 900, height = 600)

	plot(
		years,
		PXt[, "PXt"],
		"l",
		ylim = c(0, 1),
		xaxt = 'n',
		xaxs = "i",
		yaxs = "i" ,
		xlab =	"Years",
		ylab = "P(X|t)"
	)

	polygon(xx, yyint, col = 'lightgrey', border = NA)
	polygon(xx, yysd, col = 'darkgrey', border = NA)

	lines(years, PXt[, 'PXt'], "l")
	lines(years, PXt[, 'PXt.best'], "l", col="red")

	axis(1, at = seq(min(years), max(years), by = 3), las = 2)
	axis(1, at = c(min(years), max(years)), las = 2)
	dev.off()

	## Output probability that the species is extinct, P(E), at the end of the record/survey period
	year <- PXt[nrow(PXt), 'years']
	PE_min <- 1 - PXt[nrow(PXt), 'PXt.max']
	PE_mid <- 1 - PXt[nrow(PXt), 'PXt']
	PE_max <- 1 - PXt[nrow(PXt), 'PXt.min']
	PE_best <- 1 - PXt[nrow(PXt), 'PXt.best']
	PE = cbind(year, PE_min, PE_best, PE_max)
	## PE = cbind(year, PE_min, PE_mid, PE_max, PE_best)
	row.names(PE) <- Species

	cat("\nProbability of extinction:\n")
	print(PE)
	
	# save PE to a file
	write.csv(PE, EXFile)
	
	# save the P(E) plot to a file
	jpeg(Img2File, width = 600, height = 600)
		windowsFonts(A = windowsFont("arial"))
	plot(
	  family="A",
	  font.lab=2,
	  pty="s",
	  PE_best,ThreatsX[2],
	  xlim = c(0, 1),
	  ylim = c(0, 1),
	  xlab = "P(E) from Records and Surveys Model",
	  ylab = "P(E) from Threats Model"

	)
	
	PEpoly.x <- c(0.5, 0.5, 1.0, 1.0)
	PEpoly.y <- c(0.5, 1.0, 1.0, 0.5)
	polygon(PEpoly.x, PEpoly.y, col="lightpink", border = NA, density=-1)	
	EXpoly.x <- c(0.9, 0.9, 1.0, 1.0)
	EXpoly.y <- c(0.9, 1.0, 1.0, 0.9)
	polygon(EXpoly.x, EXpoly.y, col="lightcoral", border = NA, density=-3)
	
	lines(c(PE_min,PE_max),rep(ThreatsX[2],2), col="black", lwd=5)
	lines(rep(PE_best,2),c(ThreatsX[1],ThreatsX[3]), col="black", lwd=5)
	points(PE_best,ThreatsX[2], col="red3", pch = 23, lwd=6)
	text(0.95,1.01,'EX', family="A",  font.lab=2, col="red4")
	text(0.73,1.01,'CR(PE)', family="A",  font.lab=2, col="red4")
	dev.off()
	

}


