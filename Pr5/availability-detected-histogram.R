library(AHMbook)
# plot detections from point transect survey
tmp <- sim.pdata(N=150000, B=2, keep.all = TRUE)
truncation <- tmp$B
available <- tmp$d[tmp$d<=truncation]
hist(available, nclass = 100, xlab = "Distance (x)", 
     col = "grey", ylab="Number of animals",
     main = "Available (grey), detected (blue)",
     border="grey")
detected <- tmp$d[tmp$d<=truncation & tmp$y]
hist(detected, col = "blue", nclass=100,
     border="blue", add = TRUE)
