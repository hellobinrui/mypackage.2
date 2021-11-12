test_that("square works", {
  expect_equal(square(3), 9)
  expect_equal(square(1:3), c(1, 4, 9))
})
