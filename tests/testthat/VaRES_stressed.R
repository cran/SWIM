context("VaR and ES stressed")
library("SWIM")

################ stress ################
set.seed(0)
x <- as.data.frame(cbind(
  "normal" = rnorm(1000), 
  "gamma" = rgamma(1000, shape = 2)))

alpha <- c(0.8, 0.9)
q_ratio <- 1.05
s_ratio <- 1.1
k <- 1
res <- stress(type = "VaR ES", x = x, alpha = alpha, q_ratio = q_ratio, 
              s_ratio = s_ratio, k = k)

################ VaR ################
prob <- seq(0.1, 0.9, by = 0.1)

.VaR.stressed1 <- VaR_stressed(res, alpha = prob, xCol = k, wCol = 1, base = TRUE)
.VaR.stressed2 <- VaR_stressed(res, alpha = prob, xCol = k, wCol = 2, base = TRUE)
q_base <- as.numeric(quantile(x[, k], probs = prob, type = 1))
q_stressed1 <- quantile_stressed(res1, probs = prob, xCol = k, wCol = 1, type = "i/n")  
q_stressed2 <- quantile_stressed(res1, probs = prob, xCol = k, wCol = 2, type = "i/n")  

  
test_that("base", {
    # test that the values for base is equal to the R quantile function 
    expect_equal(as.numeric(.VaR.stressed1[, 2]), q_base)
    # test that the names are correct
    expect_equal(colnames(.VaR.stressed1), c(names(x)[k], "base"))
    # test that the values for base is equal to the R quantile function 
    expect_equal(as.numeric(.VaR.stressed2[, 2]), q_base)
    # test that the names are correct
    expect_equal(colnames(.VaR.stressed2), c(names(x)[k], "base"))
    
  })

test_that("stressed", {
  # test that the values for the stressed VaR is equal to 
  # the quantile_stressed function
  expect_equal(as.numeric(.VaR.stressed1[, 1]), as.numeric(q_stressed1))
  # test that the names are correct
  expect_equal(colnames(.VaR.stressed1)[1], colnames(q_stressed1))
  # test that the values for the stressed VaR is equal to 
  # the quantile_stressed function
  expect_equal(as.numeric(.VaR.stressed2[, 1]), as.numeric(q_stressed2))
  # test that the names are correct
  expect_equal(colnames(.VaR.stressed2)[1], colnames(q_stressed2))
})

################ ES ################
