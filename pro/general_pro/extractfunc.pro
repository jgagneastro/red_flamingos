;+
function extractfunc, xvals, datav, varv, profv, eval, coeffv, opvar

  gl = where(profv ne 0. and finite(profv))
  ;profv /= total(profv[gl], /NAN, /DOUBLE)
  t1 = profv[gl]/varv[gl]
  denom = total( t1 * profv[gl], /NAN, /DOUBLE) ; avoid recalc
  opt   = total( t1 * datav[gl], /NAN, /DOUBLE) / denom
  ;opvar = 1d0 / denom
  ;if abs(total(  profv[gl], /NAN, /DOUBLE) - 1d0) gt 1d-8 then stop
  opvar = total(  profv[gl], /NAN, /DOUBLE) / denom
  return, opt
  
  ;Using flux per pixel instead (f_lambda = weighted mean of "data - sky")
;  gl = where(profv ne 0. and finite(profv))
;  weights =  profv[gl]^2/varv[gl]
;  denom = total(weights, /nan, /double)
;  opt = total(weights * datav[gl], /nan, /double) / denom
;  opvar = total(weights^2 * varv[gl], /nan, /double) / denom^2
;  return, opt
  
end