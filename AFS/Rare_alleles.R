# Standard sequencing: mean coverage ~90x, pool of 400 haploid alleles
x <- sapply(1:10000, function(x) {length(unique(sample(1:400, 90, replace = TRUE)))})
mean(x)

# Deep sequencing: mean coverage ~767x, pool of 400 haploid alleles
y <- sapply(1:10000, function(x) {length(unique(sample(1:400, 767, replace = TRUE)))})
mean(y)

# Minimum detectable allele frequency at standard depth
1 / floor(mean(x))

# Minimum detectable allele frequency at deep sequencing depth
1 / floor(mean(y))

# Minimum detectable frequency at count threshold of 3 in deep sequencing
3 / 767

cat("Mean unique alleles at standard depth (~90x):", mean(x), "\n")
cat("Mean unique alleles at deep sequencing depth (~767x):", mean(y), "\n")
cat("Min detectable frequency at standard depth:", 1/floor(mean(x)), "\n")
cat("Min detectable frequency at deep sequencing depth:", 1/floor(mean(y)), "\n")
cat("Min detectable frequency at count >= 3 threshold:", round(3/767, 6), "\n")