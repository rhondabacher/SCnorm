## ----style-knitr, eval=FALSE, echo=FALSE, results='asis'-----------------
#    BiocStyle::latex()
#    render_sweave()

## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  install.packages('SCnorm_x.x.x.tar.gz', repos=NULL, type="source")
#  #OR
#  library(devtools)
#  devtools::install_github("rhondabacher/SCnorm")

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
  library(SCnorm)

## ---- eval=TRUE----------------------------------------------------------
data(ExampleData)
str(ExampleData)

## ---- eval=TRUE----------------------------------------------------------
Conditions = rep(c(1), each= 90)
str(Conditions)

## ---- eval=TRUE----------------------------------------------------------
checkCountDepth(Data = ExampleData, Conditions = Conditions, 
                OutputName = "check_exampleData", 
                FilterCellProportion = .1, NCores=3)

## ---- eval=TRUE----------------------------------------------------------

# Total Count normalization, Counts Per Million, CPM. 
ExampleData.Norm <- t((t(ExampleData) / colSums(ExampleData)) * 
                        mean(colSums(ExampleData))) 

checkCountDepth(Data = ExampleData, 
                NormalizedData = ExampleData.Norm, 
                Condition = Conditions, 
                OutputName = "check_exampleDataNorm",  
                FilterCellProportion = .1, 
                FilterExpression = 2, NCores=3)


## ---- eval=TRUE----------------------------------------------------------
Conditions = rep(c(1), each= 90)
DataNorm <- SCnorm(ExampleData, Conditions, 
                   OutputName = "MyNormalizedData",
                   SavePDF=TRUE, 
                   FilterCellNum = 10,
                   NCores=3)
str(DataNorm)

## ---- eval=FALSE---------------------------------------------------------
#  checkCountDepth(Data = umiData, Condition = Conditions,
#                  OutputName = "check_umi_scData",
#                  FilterCellProportion = .1, FilterExpression = 2)
#  
#  DataNorm <- SCnorm(umiData, Conditions,
#                     OutputName = "MyNormalizedUMIData",
#                     FilterCellNum = 10, PropToUse = .1,
#                     Thresh = .05, ditherCounts = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  DataNorm <- SCnorm(ExampleData, Conditions,
#                       OutputName = "MyNormalizedData",
#                       FilterCellNum = 10, useSpikes=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  DataNorm <- SCnorm(ExampleData, Conditions,
#                     OutputName = "MyNormalizedData",
#                     FilterCellNum = 10, withinSample = GC)
#  
#  DataNorm <- SCnorm(ExampleData, Conditions,
#                     OutputName = "MyNormalizedData",
#                     FilterCellNum = 10, withinSample = GeneLength)

## ----eval=FALSE----------------------------------------------------------
#    print(sessionInfo())

