
TDP_scores <- as.numeric(readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/FIXEDSTUFF/PrioritisationScores.txt"))
TDP_Act <- c(17, 11, 11, 5, 10)

PD_scores <- as.numeric(readLines("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood_Correct/PD_PrioritisationScores.txt"))
PD_Act <- c(17, 13, 6)

AA_scores <- as.numeric(readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/PARA/GeneExpression/Fixed_JENRA/scoreswogwas.txt"))
AA_Act <- as.numeric(readLines("/Users/clairegreen/Documents/PhD/Arthritis/Arthritis_Code/Results/PARA/GeneExpression/Fixed_JENRA/GWASscoreswogwas.txt"))


scores = TDP_scores
Act =TDP_Act
m = 100000
r <- vector()

for (i in 1:m){
  
  pick <- sample(x = scores, size = length(Act))
  r[i] <- mean(pick)
}

act <- mean(Act) 

testup <- which(r >= act) 
resultup <- sum((length(testup)+1))/(m+1) # calculate P value
resultup
mean <- mean(r)
mean
range <- range(r)
range


hist(r, 
     xlim = range(2:12),
     main = NULL, 
     xlab = "Average Score for 32 Random Genes")
abline(v = act, col = "red", lwd = 2)
