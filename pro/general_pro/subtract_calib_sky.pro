Function subtract_calib_sky, calib, NSAMPLE=nsample
  nx = n_elements(calib)
  if ~keyword_set(nsample) then $
    nsample = 45L
  baseline = fltarr(nx)+!values.f_nan
  for si=0L, nx-1L do begin & $
    si_ind = lindgen(nsample) - nsample/2L + si & $
    if si_ind[0L] lt 0L then si_ind -= min(si_ind,/nan) & $
    if si_ind[-1L] gt (nx-1L) then si_ind -= (max(si_ind) - (nx-1L)) & $
    baseline[si] = min(calib[si_ind],/nan) & $
  endfor
  ;plot, calib
  ;oplot, baseline, color=255
  return, calib-baseline
End