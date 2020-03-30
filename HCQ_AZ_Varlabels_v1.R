###################################################################
# Project:    COVID-19, HCQ-AZ 
# This file:  Data in, with Stata variable labels & codebook
# When/Who:   25 Mar 2020, AAL
###################################################################

# Prelims
# rm(list=ls())
library(foreign)
library(here)

###################################################################
# Section I.			
###################################################################

a <- read.dta("Downloads/COVID19_HCQ_v7.dta") # need v11/12 file
var.labels <- attr(a,"var.labels")
data.key <- data.frame(var.name=names(a),var.labels)
data.key

###################################################################
# Section II.		
###################################################################


