context("lagr")

test_that("str_length is number of characters", {
    expect_that(str_length("a"), equals(1))
    expect_that(str_length("ab"), equals(2))
    expect_that(str_length("abc"), equals(3))
})