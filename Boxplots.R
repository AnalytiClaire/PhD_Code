####BOXPLOT ANALYSIS OF MICROARRAY DATA####

dev.off()
par(mfrow=c(3,2))

#C9orf72
boxplot(exp_C9.LCM, las=2, 
        col = c("red","red","red","red","red","red","red",
                "red","green","green","green"), main = "C9orf72")
legend(100, 100, legend=c("Disease", "Control"), col=c("red","green"))

#CHMP2B#

boxplot(exp_CHMP2B.LCM, las=2, 
        col = c("red","red","red","red","red","red","red"
                ,"green","green","green"), main = "CHMP2B")

#sALS#

boxplot(exp_SALS.LCM, las=2, 
        col = c("red","red","red","red","red","red","red"
                ,"green","green","green"), main = "sALS")

#FTLD#

boxplot(FTLD, las=2, 
        col = c("red","red","red","red","red","red","red",
                "red","red","red","red","red","red","red",
                "red", "red","green","green","green",
                "green","green","green","green","green"), main = "FTLD")
#VCP#

boxplot(VCP, las=2, 
        col = c("red","red","red","red","red","red","red"
                ,"green","green","green"), main = "VCP")

