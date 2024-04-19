# 3 pairs uses the same number of treated and control as 
#  one 1-2 and one 2-1

# So, here are "fair" comparisons:
# 6L pairs
# 2L 1-2 plus 2L 2-1
# 3L pairs plus L 1-2 plus L 2-1 as full matching 
# Also, consider discarding the 3L pairs and doing L 1-2 plus L 2-1 in full matching 

# 2L 1-2 plus 2L 2-1 comes in first, 6L pairs comes in last, and discarding the pairs 
# in full matching is better than keeping them.


# For a simulated example, we don't need to
# distinguish 1-2 and 2-1, because one is a sign change of the other
# So the 2L 1-2 and 2L 2-1 is actually just 4L triples

library(sensitivitymv)


set.seed(1)

tau <- 0.5
delta <- tau * sqrt(2)
L <- 20000
y <- rnorm(L*6*3)
y <- matrix(y, L*6, 3)

# Add in the treatment effect
y[ ,1] <- y[ ,1] + delta

# 6L pairs
ypairs <- y[ ,1:2]

# 4L triples
ytriples <- y[1:(4*L), ]

# 3L pairs and 2L triples
yfull <- y[1:(5*L), ]
yfull[1:(3*L), 3] <- NA # 3L pairs

# Discard the 3L pairs from the full match
yfulldiscard <- yfull[(3*L+1):(5*L), ]


senmv(ypairs, gamma=3.4) 

senmv(ytriples, gamma=3.4)
senmv(ytriples, gamma=3.75)

senmv(yfull, gamma=3.4) 
senmv(yfull, gamma=3.6) 

senmv(yfulldiscard, gamma=3.4) 
senmv(yfulldiscard, gamma=3.75) 
