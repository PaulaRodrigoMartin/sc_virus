####################################################################################################################
############# Run Cell type specififc log FC calculation with healthy control as median ###########
####################################################################################################################

# FC score is calculated in addMotifExpression.R
calc_FCscores 

print("Start running FC calc")

# call a script within the R code (need to call this script within addMotifExpression.R)

# OLD: 
cat("calculating expected (median) expression values for genes\n")  
sco$medianExp <- apply(sco$exp,1,median) # do for healthy cells within each cell type
cat("calculating fold change rank values for each sample\n")
var <- apply(sco$exp,2,function(x)order(x-sco$medianExp,decreasing = TRUE)) # do for each cell type

# Run script/function within miReact
source("./code/addFCscores.R") # this code should be within addMotifExpression.R
# Need to make FC calculation into an R function.
var <- calc_FCscores(exp)