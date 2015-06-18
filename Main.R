source('EM.R')
source('bsplines.R')
source('simulation_4.R')
source('test_simulation.R')
library('ggplot2')

# 4 types , 7 basis dimension for each type  



beta <- matrix(100*rnorm(28), nrow = 7, ncol =4)

degree <- 4
knots <- nrow(beta) - degree -1
lb <- 0
ub <- 15
bb <- bspline_basis(lb,ub,knots,degree)

tbl_noise <- simulation_4(300,beta,bb)

beta_init <- matrix(100*rnorm(28),nrow = 7, ncol =4)
#beta_init <- beta + matrix(20*rnorm(28),nrow =7 ,ncol =4)

testmatrix <- EM(tbl_noise,beta_init,bb)

testdf <- test_simulation(tbl_noise,testmatrix,bb)
# 
 ggplot(testdf)+geom_line(aes(years_seen,pfvc_noise,color = 'blue')) + 
      geom_line(aes(years_seen, predict_pfvc,color ='red')) + 
       facet_wrap(~subtype) +
      scale_colour_manual(name = 'PFVC',
                          values = c('blue','red'),
                          labels = c('real','predict'))  
#                           
