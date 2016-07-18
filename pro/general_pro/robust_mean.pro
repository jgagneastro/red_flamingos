Function robust_mean, x_in, NAN=nan
  forward_function robust_bad
  return, mean(robust_bad(x_in),/nan)
End