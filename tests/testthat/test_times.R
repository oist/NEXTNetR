test_that("exponential_time", {
  tt <- exponential_time(1.0);
  r <- time_sample(1000, tt, 0, 1);
  expect_gt(ks.test(r, "pexp", rate=1)$p.value, 1e-3)
})

test_that("generic_time", {
  tt <- generic_time(
    function(tau) df(tau, df1=1, df2=1),
    function(tau) pf(tau, df1=1, df2=1, lower.tail=FALSE), FALSE,
    function(q) qf(q, df1=1, df2=1, lower.tail=FALSE), FALSE,
    NULL, 0.0
  );
  r <- time_sample(1000, tt, 0, 1);
  expect_gt(ks.test(r, "pf", df1=1, df2=1)$p.value, 1e-3)
})
