test_that("example values match closed-form reference", {
  expect_equal(deflvol(7.8, 6.4, 0.93), 9.335018, tolerance = 1e-5)
})

test_that("closed and integral methods agree", {
  CR <- c(7.8, 7.5, 8.0)
  IR <- c(6.4, 6.2, 6.5)
  defl_amp <- c(0.93, 0.90, 0.95)
  closed <- deflvol(CR, IR, defl_amp, method = "closed")
  integral <- deflvol(CR, IR, defl_amp, method = "integral")
  expect_equal(integral, closed, tolerance = 1e-4)
})

test_that("output length matches recycled input length", {
  expect_length(deflvol(c(7.8, 7.5), 6.4, 0.93), 2)
  expect_length(deflvol(numeric(0), numeric(0), numeric(0)), 0)
})

test_that("invalid inputs return NA", {
  expect_true(is.na(deflvol(NA, 6.4, 0.93)))
  expect_true(is.na(deflvol(7.8, -1, 0.93)))
  expect_true(is.na(deflvol(7.8, 6.4, Inf)))
})

test_that("separated spheres return 0", {
  expect_equal(deflvol(1, 1, -10), 0)
})
