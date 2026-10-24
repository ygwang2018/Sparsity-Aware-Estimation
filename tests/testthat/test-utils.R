test_that("beta vectorization round-trips", {
  B <- matrix(1:12, nrow = 3, byrow = TRUE)
  v <- SparsityAwareCATE:::vec_beta(B)
  B2 <- SparsityAwareCATE:::mat_beta(v, K = 3, p = 4)
  expect_equal(B2, B)
})

test_that("homogeneity structure is identified", {
  B <- rbind(c(1, 1, 0), c(1, 2, 0), c(1, 3, 0))
  s <- get_structure_from_beta(B, active_tol = 1e-8, homo_tol = 1e-8)
  expect_equal(s$homo, c(TRUE, FALSE, TRUE))
})

test_that("group-1 contrast has correct dimension", {
  a <- make_beta1_contrast(j = 2, K = 4, p = 5)
  expect_equal(dim(a), c(4L, 5L))
  expect_equal(sum(a), 1)
  expect_equal(a[1, 2], 1)
})
