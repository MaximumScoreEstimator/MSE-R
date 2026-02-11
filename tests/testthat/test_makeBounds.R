# makeBounds: box constraints for optimization

test_that("makeBounds symmetric default", {
    b <- makeBounds(5, 10)
    expect_equal(length(b$lower), 4)
    expect_equal(length(b$upper), 4)
    expect_equal(b$lower, rep(-10, 4))
    expect_equal(b$upper, rep(10, 4))
})

test_that("makeBounds asymmetric", {
    b <- makeBounds(5, 3, -4)
    expect_equal(length(b$lower), 4)
    expect_equal(length(b$upper), 4)
    expect_equal(b$lower, rep(-4, 4))
    expect_equal(b$upper, rep(3, 4))
})

test_that("makeBounds numAttrs 2", {
    b <- makeBounds(2, 1)
    expect_equal(length(b$lower), 1)
    expect_equal(length(b$upper), 1)
    expect_equal(b$lower, -1)
    expect_equal(b$upper, 1)
})
