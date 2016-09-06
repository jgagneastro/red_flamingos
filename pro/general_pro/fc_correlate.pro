Function fc_correlate, vec1, vec2, lag, WEIGHTS=weights, INTEGERS=integers
  
  forward_function interpol2
  
  ;Shifts vec2 and compares it to vec1
  ;WEIGHTS are applied w/r to vec1 and vec2 is being shifted.
  if n_params() lt 2L then $
    message, ' Problem !'
  n1 = n_elements(vec1)
  n2 = n_elements(vec2)
  if n1 ne n2 then $
    message, ' Problem !'
  x2 = dindgen(n2)
  if ~keyword_set(weights) then weights = fltarr(n1)+1.
  
  addinterp = 3.
  if ~keyword_set(lag) then $
    lag = (dindgen(n2/5.*addinterp+1)-n2/5.*addinterp/2.)/(addinterp*2.)
  
  nlag = n_elements(lag)
  corrfun = dblarr(nlag) + !values.d_nan
  for i=0L, nlag-1 do begin
    if keyword_set(integers) then $
      veclag = shift(vec2,-lag[i]) else $
      veclag = interpol2(vec2,x2,x2-lag[i], BADVALUE=!values.f_nan, /REPAIRNANS)
    vec_cross = vec1 * veclag
    ;nfinite = total(finite(vec_cross))
    ;corrfun[i] = total(vec_cross,/nan)/double(nfinite)
    weights_i = weights
    bad = where(~finite(vec_cross), nbad)
    if nbad ne 0L then weights_i[bad] = 0.
    corrfun[i] = total(vec_cross*weights,/nan)/total(weights,/nan)
  endfor
  
  return, corrfun - min(corrfun,/nan)
  
End