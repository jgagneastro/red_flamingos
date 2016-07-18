Function interpol2, fp, xp, x, _REF_EXTRA=extra, BADVALUE=badvalue, REPAIRNANS=repairnans
  forward_function interpol
  ;on_error,2
  if badvalue eq !NULL then badvalue = !values.f_nan
  ;Identify initial NaN values
  if keyword_set(repairnans) then begin
    if repairnans eq 2L then $
      bad_init = where(~finite(fp) or ~finite(shift(fp,1)) or ~finite(shift(fp,-1)), nbad_init) else $
      bad_init = where(~finite(fp), nbad_init)
  endif
  f = interpol(fp,xp,x,_EXTRA=extra, /NAN)
  bad = where(x gt max(xp,/nan) or x lt min(xp,/nan), nbad)
  ;Put back initial NaN values
  if keyword_set(repairnans) then begin
    ;For the moment, this option can only deal with arrays being *extended* everywhere, not contracted
    for i=0L, nbad_init-1L do begin
      ;If bad array element is first input array position
      if bad_init[i] eq 0L then begin
        bad_i = where(x le xp[bad_init[i]], nbad_i)
        if nbad_i ne 0L then f[bad_i] = !values.f_nan
        continue
      endif
      ;If bad array element is last input array position
      if bad_init[i] eq (n_elements(xp)-1L) then begin
        bad_i = where(x ge xp[bad_init[i]], nbad_i)
        if nbad_i ne 0L then f[bad_i] = !values.f_nan
        continue
      endif
    ;If bad array element is anywhere else
    bad_i = where(x ge xp[bad_init[i]-1L] and x le xp[bad_init[i]+1L], nbad_i)
    if nbad_i ne 0L then f[bad_i] = !values.f_nan
    endfor
  endif
  if nbad ne 0L then $
    f[bad] = badvalue[0L]
  return, f
End