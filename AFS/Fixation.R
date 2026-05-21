rbinom(1, 100, 0.001)
colnames(freqsSNPs)[which(temp$Pop == i) + 4]

C_traj_fixed <- 1:10
for (i in 1:10) {
  tempC_Traj_Pop <- which(temp$Pop == i) + 4
  tempOne <- freqsSNPs[, tempC_Traj_Pop] == 1
  tempFixed <- rowSums(tempOne) == 4
  C_traj_fixed[i] <- sum(tempFixed)
}
#potentially pairwise prop.test on if the SNPs that looked "fixed" are also (edit from future gchat  prop.test isn't really the way to check if the distribution along the genome are random)
#random between pops
#I expect across pairwise comparisions and correcting for multiple comparision tests
#we either see no significance or minor signifiance because linkage exists

#secondarily if you wanted to be conservative you could match up to selection
#SNPs under selection that looked fixed in one or two pop(s) but move in the others
#are more likely to actually be fixed
#to compare difference in means. the greater dif in mean AFC from first to last
#of not-fixed pops, the greater the likelihood of it being really fixed

#three-ish categories possible: 1) Any time they are all "1", 2) Times they are "1" and also significant for selection
#3) times they are "1" and the mean AFC of non-fixed are "high" (more of a scale perse)
#point 3) is not gaurantee to be a perfect subset of point 2)

numFixed <- C_traj_fixed[i]
hist(which(tempFixed)[2:numFixed] - which(tempFixed)[1:(numFixed - 1)])$counts
hist(which(tempFixed)[2:numFixed] - which(tempFixed)[1:(numFixed - 1)])$breaks
#most "fixed alleles" do cluster together positionally along index
#would need to be refined to be along position
convPos <- freqsSNPs$convertedPosition
hist(convPos[which(tempFixed)[2:numFixed]] - convPos[which(tempFixed)[1:(numFixed - 1)]])$breaks
hist(convPos[which(tempFixed)[2:numFixed]] - convPos[which(tempFixed)[1:(numFixed - 1)]])$counts

hist(convPos[which(tempFixed)[2:numFixed]] - convPos[which(tempFixed)[1:(numFixed - 1)]], breaks = 1000)$counts[1:50]
hist(convPos[which(tempFixed)[2:numFixed]] - convPos[which(tempFixed)[1:(numFixed - 1)]], breaks = 1000)$breaks[1:50]