Function refine_wavelength_solution, lambda, spectrum, USE_LINES=use_lines, WID=wid, POLYPAR=polypar, $
  NODISPLAY=nodisplay
  ;Refines a wavelength solution using the H lines of an A star.
  
  forward_function windows_for_reduc, gpath, sortself, readcol, deletevar, sortmultiple, unmatch, interpol2, weighted_median, mpfitpeak, remove, poly_fit, poly
  
  if ~keyword_set(nodisplay) then $
    windows_for_reduc, /MORE, /LONGTOP
  
  if ~keyword_set(wid) then $
    wid = 25L
  
  ;File containing cataloged H lines
  ;file_H = '/Users/gagne/Dropbox/IDL/IDL_Library/05-BASS/E_Spectral_Reduction/SpeXTool_dup/data/HI.dat'
  file_H = gpath('HI_FILE')
  if file_H eq '' then $
    file_H = '/Users/gagne/Documents/IDL/IDL_Library/General_Astronomy/05-BASS/E_Spectral_Reduction/SpeXTool/data/HI.dat'
  
  ;Choose which H lines to use
  if ~keyword_set(use_lines) then $
    use_lines = ['Pa !7e','Pa !7d','Pa !7c','Pa !7b','Br 15','Br 13','Br 12','Br 11','Br 10','Br !7c','Pf 20'];'Pa 09',
  sortself, use_lines, /uniq
  
  ;Read H lines
  readcol, file_H, wav, variety, num, format='D,A,A', skipline=3, /silent
  name = variety+' '+num
  deletevar, variety, num
  sortmultiple, sort(name), name, wav
  sortmultiple, uniq(name), name, wav
  
  ;Cross-match cataloged lines with those that should be used and only keep those
  unmatch, name, use_lines, ia, ib
  if n_elements(ia) eq 1 and ia[0] eq -1 then $
    message, ' No H lines which were set to be used were found in the catalog !' 
  if n_elements(ia) ne n_elements(use_lines) then $
    message, ' Some H lines which were set to be used were not found in the catalog !', /continue 
  name = name[ia]
  wav = wav[ia]
  
  ;Remove H lines outside of range
  bad = where(wav lt min(lambda,/nan) or wav gt max(lambda,/nan), nbad)
  if nbad eq n_elements(wav) then $
    message, ' All H lines are outside of spectral range !'
  if nbad ne 0L then remove, bad, wav, name
  
  ;Determine pixel position of each line
  npix = n_elements(lambda)
  linepos = interpol2(findgen(npix),lambda,wav)
  bad = where(~finite(linepos), nbad)
  if nbad eq n_elements(wav) then $
    message, ' All pixel positions of H lines are NaNs !'
  if nbad ne 0L then remove, bad, wav, name
  
  ;Create a smoothed version of the spectrum
  mspectrum = median(spectrum,3)
  
  ;Loop on lines to perform the fits
  nlines = n_elements(name)
  pars = fltarr(5L,nlines)+!values.f_nan
  if ~keyword_set(nodisplay) then $
    wset, 2
  for i=0L, nlines-1L do begin
    gg = lindgen(wid+1)-wid/2 + round(linepos[i])
    bad = where(gg lt 0 or gg gt (npix-1L), nbad)
    if nbad eq n_elements(gg) then $
      message, 'Window was totally outside of range ! (This should not happen)'
    if nbad ne 0L then remove, bad, gg 
    
    ;Guess parameters of line profile
    guess_amp = weighted_median(spectrum[gg],medval=.1) - weighted_median(spectrum[gg],medval=.9)
    midg = median(gg)
    min_amp = (min(spectrum[gg],/nan) - max(spectrum[gg],/nan))*1.5
    max_amp = 0.
    void = min(abs(spectrum[gg] - weighted_median(spectrum[gg],medval=.1)),wmin)
    guess_pos = lambda[gg[wmin]]
    min_pos = (lambda[min(gg,/nan)+2L] < lambda[max(gg,/nan)-2L])
    max_pos = (lambda[min(gg,/nan)+2L] > lambda[max(gg,/nan)-2L])
    guess_pos >= min_pos
    guess_pos <= max_pos
    guess_wid_pix = 10L
    min_wid_pix = 4L
    max_wid_pix = wid
    guess_wid = abs(lambda[(round(midg)+guess_wid_pix/2)<(npix-1)]-lambda[(round(midg)-guess_wid_pix/2)>0])
    min_wid = abs(lambda[(round(midg)+min_wid_pix/2)<(npix-1)]-lambda[(round(midg)-min_wid_pix/2)>0])
    max_wid = abs(lambda[(round(midg)+max_wid_pix/2)<(npix-1)]-lambda[(round(midg)-max_wid_pix/2)>0])
    guess_slope = (median(spectrum[gg[-5:*]])-median(spectrum[gg[0:4]]))/(lambda[gg[-1L]] - lambda[gg[0L]])
    min_slope = -2.*abs(guess_slope)
    max_slope = 2.*abs(guess_slope)
    guess_continuum = weighted_median(spectrum[gg],medval=.9) - median(lambda[gg])*guess_slope
    min_continuum = min(spectrum[gg],/nan) - median(lambda[gg])*guess_slope
    max_continuum = max(spectrum[gg],/nan) - median(lambda[gg])*guess_slope
    p0 = [guess_amp, guess_pos, guess_wid, guess_continuum,guess_slope]
    
    ;Fit a line profile
    parinfo = replicate({VALUE:!values.f_nan,FIXED:0B,LIMITED:[1B,1B],LIMITS:[0.,0.]},n_elements(p0))
    parinfo.VALUE = p0
    parinfo[0].LIMITS = [min_amp,max_amp]
    parinfo[1].LIMITS = [min_pos,max_pos]
    parinfo[2].LIMITS = [min_wid,max_wid]
    parinfo[3].LIMITS = [min_continuum,max_continuum]
    parinfo[4].LIMITS = [min_slope,max_slope]
    par = !NULL
    yfit = mpfitpeak(lambda[gg],mspectrum[gg], par, NTERMS=n_elements(p0), WEIGHTS=1D, /QUIET, ESTIMATES=parinfo.VALUE, PARINFO=parinfo, STATUS=status, ERRMSG=errmsg, /LORENTZ)
    if errmsg ne '' then stop
    
    if abs(par[0]) lt abs(min_amp/50) then continue
    
    ;Save parameters
    pars[*,i] = par
    
    if ~keyword_set(nodisplay) then begin
      ;Plot the best fit
      plot, lambda[gg], mspectrum[gg], yrange=[min(spectrum[gg],/nan),max(spectrum[gg],/nan)], ystyle=1
      ;oplot, lambda[gg], mpfitpeak_lorentz(lambda[gg],p0), color=rvb_hex(100,255,100)
      oplot, lambda[gg], yfit, color=255
      stop
    endif
  endfor
  
  ;Fit a 2nd-order relation to the differences
  pos = reform(pars[1,*])
  sortmultiple, sort(wav), wav, name, pos
  bad = where(~finite(pos), nbad)
  if nbad ne 0L then remove, bad, wav, name, pos
  polypar = reform(poly_fit(pos,wav,2))
  
  if ~keyword_set(nodisplay) then begin
    xv = findgen(2000)/1999.*(max(pos,/nan)*1.01-min(pos,/nan)*.99)+min(pos,/nan)*.99
    wset, 0
    plot, xv, poly(xv,polypar)
    oplot,pos, wav, /ps, color=255,symsize=2
    wset, 1
    plot, xv, poly(xv,polypar)-xv
    oplot,pos, wav-pos, /ps, color=255,symsize=2
    stop
    
    wset, 0
    plot, poly(lambda,polypar), spectrum, xrange=[min(lambda,/nan),mean(lambda,/nan)]
    oplot, wav, interpol2(spectrum,poly(lambda,polypar),wav), /ps, color=255, symsize=3
    wset, 1
    plot, poly(lambda,polypar), spectrum, xrange=[mean(lambda,/nan),max(lambda,/nan)]
    oplot, wav, interpol2(spectrum,poly(lambda,polypar),wav), /ps, color=255, symsize=3
    stop
    
    windows_for_reduc, /MORE
  endif
  
  return, poly(lambda,polypar)
End