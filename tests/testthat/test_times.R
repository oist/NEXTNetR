test_that("exponential_time", {
  tt <- episimR_exponential_time(1.0);
  r <- episimR_time_sample(100, tt, 0, 1);
  print(r)
})