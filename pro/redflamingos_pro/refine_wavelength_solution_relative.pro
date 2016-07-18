Function refine_wavelength_solution_relative, lambda_obs, spectrum_obs_in, lambda_ref, spectrum_ref_in, $
  USE_LINES=use_lines, WID=wid, POLYPAR=polypar, NODISPLAY=nodisplay, MINAMP=minamp, SLOPE=slope, $
  MINWID=minwid, MAXWID=maxwid, dmaxfit=dmaxfit, maxdxref=maxdxref, REMOVECONT=removecont
  ;Refines a wavelength solution by comparing the absorption lines of two spectra at a given set of wavelengths
  ;; By default, telluric absorption lines will be used, unless "USE_LINES" is specified.
  ;LAMBDA1, SPETRUM1 is  
  
  forward_function savitzky_golay, windows_for_reduc, sortself, weighted_median, mpfitpeak, remove, sortmultiple, poly_fit, supersmooth, rvb_hex
  
  if ~keyword_set(nodisplay) then begin
    windows_for_reduc, /MORE, /LONGTOP
    device, /decomposed
  endif
  
  if keyword_set(removecont) then begin
    spectrum_obs = spectrum_obs_in/savitzky_golay(spectrum_obs_in,80L)
    spectrum_ref = spectrum_ref_in/savitzky_golay(spectrum_ref_in,80L)
  endif else begin
    spectrum_obs = spectrum_obs_in
    spectrum_ref = spectrum_ref_in
  endelse
  
  if ~keyword_set(wid) then $
    wid = 8L
  if ~keyword_set(minamp) then $
    minamp = !values.f_nan
  if ~keyword_set(minwid) then $
    minwid = 4
  if ~keyword_set(maxwid) then $
    maxwid = 25
  
  if ~keyword_set(use_lines) then $
    use_lines = [1.133,1.394,1.431,1.47,1.82125,1.9205,1.957,2.005,2.368,2.4315];,1.399;,1.381;,1.4448;,1.4185;,1.45;,1.8456;,1.906;1.9,;,1.885;,2.055;,1.767;,1.268
  sortself, use_lines, /uniq

  ;Remove H lines outside of range
  bad = where(use_lines lt min(lambda_obs,/nan) or use_lines gt max(lambda_obs,/nan) or use_lines lt min(lambda_ref,/nan) or use_lines gt max(lambda_ref,/nan), nbad)
  if nbad eq n_elements(use_lines) then $
    message, ' All lines are outside of spectral range !'
  if nbad ne 0L then remove, bad, use_lines

  ;Determine pixel position of each line
  npix_obs = n_elements(lambda_obs)
  npix_ref = n_elements(lambda_ref)
  linepos_obs = interpol2(findgen(npix_obs),lambda_obs,use_lines)
  linepos_ref = interpol2(findgen(npix_ref),lambda_ref,use_lines)
  bad = where(~finite(linepos_obs) or ~finite(linepos_ref), nbad)
  if nbad eq n_elements(use_lines) then $
    message, ' All pixel positions of lines are NaNs !'
  if nbad ne 0L then remove, bad, use_lines

  ;Create a smoothed version of the spectra
  mspectrum_obs = spectrum_obs
  mspectrum_ref = spectrum_ref
  ;mspectrum_obs = median(spectrum_obs,3)
  ;mspectrum_ref = median(spectrum_ref,3)

  ;Loop on lines to perform the fits on the "observed" spectrum
  nlines = n_elements(use_lines)
  pars = fltarr(5L,nlines,2L)+!values.f_nan
  if ~keyword_set(nodisplay) then wset, 2
  spectrum = spectrum_obs
  mspectrum = mspectrum_obs
  npix = npix_obs
  linepos = linepos_obs
  lambda = lambda_obs
  for i=0L, nlines-1L do begin
    gg = lindgen(wid+1)-wid/2 + round(linepos[i])
    bad = where(gg lt 0 or gg gt (npix-1L), nbad)
    if nbad eq n_elements(gg) then $
      message, 'Window was totally outside of range ! (This should not happen)'
    if nbad ne 0L then remove, bad, gg
    ngg = n_elements(gg)
    
    ;Guess parameters of line profile
    guess_amp = weighted_median(spectrum[gg],medval=.1) - weighted_median(spectrum[gg],medval=.9)
    midg = median(gg)
    min_amp = (min(spectrum[gg],/nan) - max(spectrum[gg],/nan))*1.5
    max_amp = 0.
    ;void = min(abs(spectrum[gg] - weighted_median(spectrum[gg],medval=.05)),wmin)
    void = min(spectrum[gg],wmin,/nan)
    guess_pos = lambda[gg[wmin]]
    
    if wmin lt 1L or wmin gt (ngg-2L) then begin
      message, ' Skipped a line on the edge', /continue
      continue
    endif
    min_pos0 = (lambda[(gg[wmin]-1L)>0L] < lambda[(gg[wmin]+1L)<(npix_obs-1L)])
    max_pos0 = (lambda[(gg[wmin]-1L)>0L] > lambda[(gg[wmin]+1L)<(npix_obs-1L)])
    min_pos = guess_pos - abs(max_pos0-min_pos0)/4.
    max_pos = guess_pos + abs(max_pos0-min_pos0)/4.
    ;min_pos = (lambda[min(gg,/nan)+2L] < lambda[max(gg,/nan)-2L])
    ;max_pos = (lambda[min(gg,/nan)+2L] > lambda[max(gg,/nan)-2L])
    ;guess_pos >= min_pos
    ;guess_pos <= max_pos
    min_wid_pix = minwid
    max_wid_pix = maxwid
    guess_wid_pix = 10.^((alog10(maxwid)+alog10(minwid))/2.)
    guess_wid = abs(lambda[(round(midg)+guess_wid_pix/2)<(npix-1)]-lambda[(round(midg)-guess_wid_pix/2)>0])
    min_wid = abs(lambda[(round(midg)+min_wid_pix/2)<(npix-1)]-lambda[(round(midg)-min_wid_pix/2)>0])
    max_wid = abs(lambda[(round(midg)+max_wid_pix/2)<(npix-1)]-lambda[(round(midg)-max_wid_pix/2)>0])
    if keyword_set(slope) then begin
      guess_slope = (median(spectrum[gg[-5:*]])-median(spectrum[gg[0:4]]))/(lambda[gg[-1L]] - lambda[gg[0L]])
      min_slope = -2.*abs(guess_slope)
      max_slope = 2.*abs(guess_slope)
    endif else begin
      guess_slope = 0.
      min_slope = -1.
      max_slope = 1.
    endelse
    guess_continuum = weighted_median(spectrum[gg],medval=.9) - median(lambda[gg])*guess_slope
    min_continuum = min(spectrum[gg],/nan) - median(lambda[gg])*guess_slope
    max_continuum = max(spectrum[gg],/nan)*1.5 - median(lambda[gg])*guess_slope
    p0 = [guess_amp, guess_pos, guess_wid, guess_continuum,guess_slope]

    ;Fit a line profile
    parinfo = replicate({VALUE:!values.f_nan,FIXED:0B,LIMITED:[1B,1B],LIMITS:[0.,0.]},n_elements(p0))
    parinfo.VALUE = p0
    parinfo[0].LIMITS = [min_amp,max_amp]
    parinfo[1].LIMITS = [min_pos,max_pos]
    parinfo[2].LIMITS = [min_wid,max_wid]
    parinfo[3].LIMITS = [min_continuum,max_continuum]
    parinfo[4].LIMITS = [min_slope,max_slope]
    if ~keyword_set(slope) then $
      parinfo[4].FIXED = 1
    par = !NULL
    yfit = mpfitpeak(lambda[gg],mspectrum[gg], par, NTERMS=n_elements(p0), WEIGHTS=1D, /QUIET, ESTIMATES=parinfo.VALUE, PARINFO=parinfo, STATUS=status, ERRMSG=errmsg, /LORENTZ)
    if errmsg ne '' then begin
      message, ' One peak fit has failed !', /continue
      continue
      stop
    endif
    ;print,abs(par[1]-guess_pos),maxdxref,abs(par[0]),minamp
    
;    if abs(par[0]) lt minamp then $
;      continue
    
    ;Save parameters
    pars[*,i,0L] = par
    
    if ~keyword_set(nodisplay) then begin
      ;Plot the best fit
      plot, lambda[gg], mspectrum[gg], yrange=[min(spectrum[gg],/nan),max(spectrum[gg],/nan)], ystyle=1
      ;oplot, lambda[gg], mpfitpeak_lorentz(lambda[gg],p0), color=rvb_hex(100,255,100)
      oplot, lambda[gg], yfit, color=255
      stop
    endif
  endfor
  
  ;Loop on lines to perform the fits on the "reference" spectrum
  if ~keyword_set(nodisplay) then wset, 2
  spectrum = spectrum_ref
  mspectrum = mspectrum_ref
  npix = npix_ref
  linepos = linepos_ref
  lambda = lambda_ref
  for i=0L, nlines-1L do begin
    gg = lindgen(wid+1)-wid/2 + round(linepos[i])
    bad = where(gg lt 0 or gg gt (npix-1L), nbad)
    if nbad eq n_elements(gg) then $
      message, 'Window was totally outside of range ! (This should not happen)'
    if nbad ne 0L then remove, bad, gg
    ngg = n_elements(gg)
    
    ;Guess parameters of line profile
    guess_amp = weighted_median(spectrum[gg],medval=.1) - weighted_median(spectrum[gg],medval=.9)
    midg = median(gg)
    min_amp = (min(spectrum[gg],/nan) - max(spectrum[gg],/nan))*1.5
    max_amp = 0.
    void = min(spectrum[gg],wmin,/nan)
    ;void = min(abs(spectrum[gg] - weighted_median(spectrum[gg],medval=.05)),wmin)
    guess_pos = lambda[gg[wmin]]
    if wmin lt 1L or wmin gt (ngg-2L) then continue
    min_pos0 = (lambda[(gg[wmin]-1L)>0L] < lambda[(gg[wmin]+1L)<(npix_obs-1L)])
    max_pos0 = (lambda[(gg[wmin]-1L)>0L] > lambda[(gg[wmin]+1L)<(npix_obs-1L)])
    min_pos = guess_pos - abs(max_pos0-min_pos0)/4.
    max_pos = guess_pos + abs(max_pos0-min_pos0)/4.
    ;min_pos = (lambda[min(gg,/nan)+2L] < lambda[max(gg,/nan)-2L])
    ;max_pos = (lambda[min(gg,/nan)+2L] > lambda[max(gg,/nan)-2L])
    ;guess_pos >= min_pos
    ;guess_pos <= max_pos
    min_wid_pix = minwid
    max_wid_pix = maxwid
    guess_wid_pix = 10.^((alog10(maxwid)+alog10(minwid))/2.)
    guess_wid = abs(lambda[(round(midg)+guess_wid_pix/2)<(npix-1)]-lambda[(round(midg)-guess_wid_pix/2)>0])
    min_wid = abs(lambda[(round(midg)+min_wid_pix/2)<(npix-1)]-lambda[(round(midg)-min_wid_pix/2)>0])
    max_wid = abs(lambda[(round(midg)+max_wid_pix/2)<(npix-1)]-lambda[(round(midg)-max_wid_pix/2)>0])
    if keyword_set(slope) then begin
      guess_slope = (median(spectrum[gg[-5:*]])-median(spectrum[gg[0:4]]))/(lambda[gg[-1L]] - lambda[gg[0L]])
      min_slope = -2.*abs(guess_slope)
      max_slope = 2.*abs(guess_slope)
    endif else begin
      guess_slope = 0.
      min_slope = -1.
      max_slope = 1.
    endelse
    guess_continuum = weighted_median(spectrum[gg],medval=.9) - median(lambda[gg])*guess_slope
    min_continuum = min(spectrum[gg],/nan) - median(lambda[gg])*guess_slope
    max_continuum = max(spectrum[gg],/nan)*1.5 - median(lambda[gg])*guess_slope
    p0 = [guess_amp, guess_pos, guess_wid, guess_continuum,guess_slope]

    ;Fit a line profile
    parinfo = replicate({VALUE:!values.f_nan,FIXED:0B,LIMITED:[1B,1B],LIMITS:[0.,0.]},n_elements(p0))
    parinfo.VALUE = p0
    parinfo[0].LIMITS = [min_amp,max_amp]
    parinfo[1].LIMITS = [min_pos,max_pos]
    parinfo[2].LIMITS = [min_wid,max_wid]
    parinfo[3].LIMITS = [min_continuum,max_continuum]
    parinfo[4].LIMITS = [min_slope,max_slope]
    if ~keyword_set(slope) then $
      parinfo[4].FIXED = 1
    par = !NULL
    yfit = mpfitpeak(lambda[gg],mspectrum[gg], par, NTERMS=n_elements(p0), WEIGHTS=1D, /QUIET, ESTIMATES=parinfo.VALUE, PARINFO=parinfo, STATUS=status, ERRMSG=errmsg, /LORENTZ)
    if errmsg ne '' then stop
    
    if keyword_set(maxdxref) then $
      if abs(par[1]-guess_pos) gt maxdxref then $
        continue
    ;if abs(par[0]) lt minamp then $
    ;  continue
    
    ;Save parameters
    pars[*,i,1L] = par
    
    if ~keyword_set(nodisplay) then begin
      ;Plot the best fit
      plot, lambda[gg], mspectrum[gg], yrange=[min(spectrum[gg],/nan),max(spectrum[gg],/nan)], ystyle=1
      ;oplot, lambda[gg], mpfitpeak_lorentz(lambda[gg],p0), color=rvb_hex(100,255,100)
      oplot, lambda[gg], yfit, color=255
      stop
    endif
  endfor
  
  ;Fit a 2nd-order relation to the differences
  pos_obs = reform(pars[1,*,0L])
  pos_ref = reform(pars[1,*,1L])
  bad = where(~finite(pos_obs) or ~finite(pos_ref), nbad)
  if nbad eq n_elements(pos_obs) then message, ' All line fits have failed !'
  if nbad ne 0 then remove, bad, pos_obs, pos_ref
  sortmultiple, sort(pos_ref), pos_ref, pos_obs
  polypar = reform(poly_fit(pos_obs,pos_ref,2))
  
  if keyword_set(dmaxfit) then begin
    bad = where(abs(poly(pos_obs,polypar)-pos_ref) gt dmaxfit, nbad)
    if nbad ne 0 and (n_elements(pos_obs)-nbad) gt 2L then remove, bad, pos_obs, pos_ref
    polypar = reform(poly_fit(pos_obs,pos_ref,2L))
  endif
  
  if ~keyword_set(nodisplay) then begin
    xv = findgen(2000)/1999.*(max(pos_obs,/nan)*1.01-min(pos_obs,/nan)*.99)+min(pos_obs,/nan)*.99
    wset, 0
    plot, xv, poly(xv,polypar)
    oplot, pos_obs, pos_ref, /ps, color=255,symsize=2
    wset, 1
    plot, [xv,pos_obs], [poly(xv,polypar)-xv,pos_ref-pos_obs], /ps, /nodata
    oplot, xv, poly(xv,polypar)-xv
    oplot, pos_obs, pos_ref-pos_obs, /ps, color=255,symsize=2
    stop
  
    wset, 0
    plot, poly(lambda_obs,polypar), spectrum_obs, xrange=[min(lambda_obs,/nan),mean(lambda_obs,/nan)]
    oplot, pos_ref, interpol2(spectrum_obs,poly(lambda_obs,polypar),pos_ref), /ps, color=255,symsize=3
    wset, 1
    plot, poly(lambda_obs,polypar), spectrum_obs, xrange=[mean(lambda_obs,/nan),max(lambda_obs,/nan)]
    oplot, pos_ref, interpol2(spectrum_obs,poly(lambda_obs,polypar),pos_ref), /ps, color=255,symsize=3
    stop
    
    wset, 0
    plot, poly(lambda_obs,polypar), spectrum_obs/supersmooth(spectrum_obs,50), xrange=[min(lambda_obs,/nan),mean(lambda_obs,/nan)]
    oplot, lambda_ref, spectrum_ref/supersmooth(spectrum_ref,50), color=255
    oplot, pos_ref, pos_ref*0+1, /ps, color=rvb_hex(100,255,100),symsize=3
    wset, 1
    plot, poly(lambda_obs,polypar), spectrum_obs/supersmooth(spectrum_obs,50), xrange=[mean(lambda_obs,/nan),max(lambda_obs,/nan)]
    oplot, lambda_ref, spectrum_ref/supersmooth(spectrum_ref,50), color=255
    oplot, pos_ref, pos_ref*0+1, /ps, color=rvb_hex(100,255,100),symsize=3
    stop
    
    windows_for_reduc, /MORE
  endif
  
  return, poly(lambda_obs,polypar)
End