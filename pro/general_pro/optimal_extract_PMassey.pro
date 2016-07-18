Pro optimal_extract_PMassey, trace, sky, profile, readnoise, spext, NITER=niter, NOSKY=nosky, ROBUST_SHIFT=robust_shift, $
  RSHIFT_MAXDEL=rshift_maxdel, RSHIFT_NPRO=rshift_npro
  
  forward_function weighted_median, robust_sigma, poly, ladfit
  
  ;If the robust_shift method is specified, create an array of shifted profiles
  if keyword_set(robust_shift) then begin
    ;Maximum absolute shift in pixels
    if ~keyword_set(rshift_maxdel) then $
      rshift_maxdel = 5d0
    if ~keyword_set(rshift_npro) then $
      rshift_npro = 600L
    ;Ensure that rshift_npro is odd, such that "zero shift" exists
    if (rshift_npro mod 2) eq 0L then $
      rshift_npro += 1L
    shift_vec = (dindgen(rshift_npro)/double(rshift_npro-1L)-0.5d0) * rshift_maxdel
    profile1d = median(profile,dim=1)
    nx = n_elements(profile1d)
    profiles = dblarr(nx,rshift_npro) + !values.d_nan
    for ii=0L, rshift_npro-1L do begin & $
      if shift_vec[ii] eq 0d0 then begin & $
        profiles[*,ii] = profile1d & $
        continue & $
      endif & $
      profiles[*,ii] = fshift2(profile1d,shift_vec[ii],/MASK_EDGES,MASK_VALUE=0d0) & $
    endfor
  endif
  
  ;Create a sky-subtracted image
  trace_sky = trace - sky
  
  ;Create a version without the object in it
  trace_empty = trace
  bad = where(abs(trace_sky) gt weighted_median(abs(trace_sky),medval=.8)*3., nbad)
  if nbad ne 0 then trace_empty[bad] = !values.f_nan
  
  ;Determine pixel-to-pixel stddev of variations of this empty sky image.
  ; This will serve as the error on the sky determination
  err_sky = robust_sigma(((trace_empty-shift(trace_empty,1,0))[1:*,1:*]+(trace_empty-shift(trace_empty,0,1))[1:*,1:*])/2.)
  
  ;Loop on spatial pixels
  npix = (size(trace))[1]
  xi = lindgen((size(trace))[2])
  spext = findgen(npix)+!values.f_nan
  if ~keyword_set(robust_shift) then nj = 1L else $
    nj = rshift_npro
  spext_shifts_2D = findgen(npix,nj)+!values.f_nan
  for i=0L, npix-1L do begin
    
    spext_shifts = dblarr(nj)+!values.d_nan
    for j=0L, nj-1L do begin
      
      if keyword_set(robust_shift) then begin
        profile_i = profiles[*,j]
      endif else begin
        profile_i = reform(profile[i,*])
      endelse
      
      tracei = reform(trace_sky[i,*])
      ;Fit the sky *again*
      gg = where(tracei lt weighted_median(tracei,medval=.6)*2., ngg)
      if keyword_set(nosky) or ngg eq 0L then begin
        skyi = fltarr(n_elements(xi))
      endif else begin
        skyi = poly(xi,ladfit(xi[gg],tracei[gg]))
      endelse
      tracei -= skyi
      ;Determine variance
      vari = readnoise^2+reform(trace[i,*])+skyi+err_sky^2+reform(sky[i,*])
      weights = reform(profile_i)^2/vari
      bad = where(~finite(tracei*weights),nbad)
      if nbad ne 0L then begin
        tracei[bad] = !values.f_nan
        weights[bad] = !values.f_nan
      endif
      spexti = total(tracei*weights,/nan)/total(weights,/nan)
      ;If NITER is specified, do iterative treatment of bad pixels
      if keyword_set(niter) then begin
        for j=0L, niter-1L do begin
          varij = readnoise^2+spexti*reform(profile_i)+skyi+err_sky^2+reform(sky[i,*])
          weightsj = reform(profile_i)^2/varij
          bad = where(~finite(tracei*weightsj),nbad)
          if nbad ne 0L then begin
            tracei[bad] = !values.f_nan
            weightsj[bad] = !values.f_nan
          endif
          spexti = total(tracei*weightsj,/nan)/total(weightsj,/nan)
        endfor;Loop on cosmic ray rejection
      endif
      spext_shifts[j] = spexti
      spext_shifts_2D[i,j] = spexti
    endfor;Loop on "robust" shifted profiles 
    if nj eq 1L then $
      spext[i] = spext_shifts[0L]
    if nj ne 1L then $
      spext[i] = max(spext_shifts,/NAN)
  endfor;Loop on spatial pixels
  
  if keyword_set(robust_shift) then begin
    nfilt = 15L
    void = max(spext_shifts_2d,dim=2,wmax,/NAN)
    shift_max_ind = reform((array_indices(spext_shifts_2d,wmax))[1,*])
    shift_max_ind = long(median(double(shift_max_ind),nfilt))
    shift_max = shift_vec[shift_max_ind]
    i_shift_max_ind = lindgen(n_elements(shift_max_ind))
    fitc = ladfit(i_shift_max_ind[(nfilt):(-nfilt-1L)],double(shift_max_ind[(nfilt):(-nfilt-1L)]))
    shift_max_ind_fitted = long(round(poly(i_shift_max_ind,fitc)))
    for i=0L, npix-1L do $
      spext[i] = spext_shifts_2d[i,shift_max_ind_fitted[i]]
    shift_max_fitted = shift_vec[shift_max_ind_fitted]
    ;if total(abs(shift_max_fitted)) ne 0 then stop
  endif
  
End