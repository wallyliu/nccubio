# Yu-Wen Liu's bioInfo Homework1
# 2016/09/29

# read PAM1 from data
pam1 <- read.table("/Users/wally/Desktop/NCCU_master/1051/bioinformatric/nccubio/HW1/pam1.txt")

# construct PAM250 from PAM1
pam250 <- diag(20)
for( i in 1:250){
  pam250 <- pam250 %*% as.matrix(pam1/10000) 
}

# output PAM250 as a file
pam250 <- round( as.data.frame(pam250) * 10000, digit = 4)
pam250 <- rbind( colnames(pam250), pam250)
colnames(pam250) <- as.factor(colnames(pam1))
rownames(pam250) <- as.factor(c(" ", rownames(pam1)))
write.table( pam250, file = "pam250.txt", sep = "\t", col.names= FALSE)

