Function find_1d_peaks, sp_in, FLUX_THRESHOLD=flux_threshold, SIG_THRESHOLD=sig_threshold, DISPLAY=display
  ;Finds peaks in a lamp calibration 1D spectrum
  if ~keyword_set(sp_in) then $
    message, 'You must enter a calibration spectrum (1D) !'
  if ~keyword_set(flux_threshold) then $
    flux_threshold = .01
  if ~keyword_set(sig_threshold) then $
    sig_threshold = 3.
  
  ;Normalize spectrum
  sp = sp_in
  sp /= max(sp,/nan)
  
  ;Find all local maxima
  locmax = where(sp[1:-2] gt (shift(sp,1))[1:-2] and sp[1:-2] gt (shift(sp,-1))[1:-2], nlocmax) + 1L
  if nlocmax eq 0 then message, ' Cannot find a single local maximum !'
  
  ;Mask lines
  clevel = where(sp lt flux_threshold, nclevel)
  if nclevel eq 0 then message, ' Cannot find a region without lines !'
  
  ;Compute standard deviation where lines are masked
  std = stddev(sp[clevel],/nan)
  
  ;Clip all local maxima below the N_SIGMA threshold
  bad = where(sp[locmax] lt std*sig_threshold, nbad)
  if nbad eq nlocmax then message, ' All local maxima are below the N_SIGMA threshold !'
  if nbad ne 0L then remove, bad, locmax
  nlocmax = n_elements(locmax)
  
  ;Display results if needed
  if keyword_set(display) then begin
    plot, sp
    oplot,locmax,sp[locmax],color=255,/ps
  endif
  
  ;Return peak positions
  return, locmax
  
End