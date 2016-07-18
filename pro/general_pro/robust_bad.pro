Function robust_bad, vec_in, BAD_VALUE=bad_value, INDICES=indices, NSIG=nsig, SILENT=silent
;  Enleve les positions des vecteurs qui sont trop deviantes, comme robust_sigma utilise
  
  compile_opt hidden
  forward_function avg
  
  if ~keyword_set(nsig) then nsig = 6.
  
  vec = vec_in
  Y = vec
  EPS = 1.0E-20
  Y0  = MEDIAN(Y,/EVEN)
  MAD = MEDIAN( ABS(Y-Y0), /EVEN )/0.6745
  
  ; If the MAD=0, try the MEAN absolute deviation:
  IF MAD LT EPS THEN MAD = AVG( ABS(Y-Y0) )/.80
  IF MAD LT EPS THEN BEGIN
    if ~keyword_set(silent) then message, 'ROBUST_BAD : Distribution trop bizarre !', /CONTINUE
    if keyword_set(indices) then return, -1
    return, vec_in
  ENDIF
  
  ; Now the biweighted value:
  U   = (Y-Y0)/(nsig*MAD)
  UU  = U*U
  BAD   = WHERE(UU GT 1.0, BADCOUNT)
  if keyword_set(indices) then return, bad
  
  if ~keyword_set(bad_value) then bad_value = !values.f_nan
  if badcount eq n_elements(vec) then message, 'ROBUST_BAD : Distribution trop bizarre !'
  if badcount ne 0 then vec[bad] = bad_value 
  
  return, vec
  
End