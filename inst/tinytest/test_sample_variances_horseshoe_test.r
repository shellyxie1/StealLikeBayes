#This is a test file

# test_sample_horseshoe_variances.R
# Unit tests for the horseshoe prior variance sampler

# Test 1: Basic functionality - function runs without errors
  p <- 10
  coefs <- rnorm(p, 0, 0.5)
  theta <- rep(1, p)
  zeta <- 1
  nu <- rep(1, p)
  varpi <- 1
  
  expect_no_error(
    V_i <- sample_variances_horseshoe(coefs, theta, zeta, nu, varpi)
  )

# Test 2: Return structure
test_that("Function returns a vector of correct length", {
  p <- 10
  coefs <- rnorm(p, 0, 0.5)
  theta <- rep(1, p)
  zeta <- 1
  nu <- rep(1, p)
  varpi <- 1
  
  V_i <- sample_variances_horseshoe(coefs, theta, zeta, nu, varpi)
  
  expect_type(V_i, "double")
  expect_length(V_i, p)
  expect_true(all(V_i > 0))
})


# Test 3: Shrinkage behavior - large coefficients get less shrinkage
test_that("Large coefficients receive less shrinkage than small coefficients", {
  set.seed(123)
  p <- 50
  
  # Create coefficients with clear signal/noise separation
  coefs <- c(rep(3, 5), rep(0.01, 45))  # 5 large, 45 tiny
  
  theta <- rep(1, p)
  zeta <- 1
  nu <- rep(1, p)
  varpi <- 1
  
  # Run multiple iterations to get stable estimates
  for(i in 1:100) {
    V_i <- sample_variances_horseshoe(coefs, theta, zeta, nu, varpi)
    # theta, zeta, nu, varpi are updated by reference automatically
  }
  
  # Mean variance for large coefficients should be larger than for small ones
  mean_large <- mean(theta[1:5])
  mean_small <- mean(theta[6:50])
  
  expect_gt(mean_large, mean_small)
})


# Test 4: Reproducibility with set.seed
test_that("Results are reproducible with set.seed", {
  p <- 10
  coefs <- rnorm(p, 0, 0.5)
  
  # First run
  theta1 <- rep(1, p)
  zeta1 <- 1
  nu1 <- rep(1, p)
  varpi1 <- 1
  
  set.seed(42)
  V_i1 <- sample_variances_horseshoe(coefs, theta1, zeta1, nu1, varpi1)
  
  # Second run
  theta2 <- rep(1, p)
  zeta2 <- 1
  nu2 <- rep(1, p)
  varpi2 <- 1
  
  set.seed(42)
  V_i2 <- sample_variances_horseshoe(coefs, theta2, zeta2, nu2, varpi2)
  
  expect_equal(V_i1, V_i2)
  expect_equal(zeta1, zeta2)
  expect_equal(theta1, theta2)
})


# Run all tests
cat("\n=== Running Horseshoe Variance Sampler Tests ===\n\n")
test_dir(".", reporter = "summary")