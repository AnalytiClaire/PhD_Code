## Establishing statistical significance of John's 22 static module set ##

setwd('/Users/clairegreen/Documents/PhD/TDP-43') 
static_mod <- read.csv("staticFImodules_geneNames.csv")
J_stat_mod <- read.csv("J_Static_Mod.csv")
stat_mod_name <- static_mod[,1]

####check against random permutatuons

m=10000 #number of repetitions 
r <- c(1:m) #store repetition numbers in vector "r"

for (j in 1:m)
{
  random <- sample (stat_mod_name, size=22, replace=F)
  random <- length (which (Tracking.GM %in% random))
  r[j] <- random
}

test1 <- which (r > test) 
enrich [i,2] <- (length(test1)/m)

}

