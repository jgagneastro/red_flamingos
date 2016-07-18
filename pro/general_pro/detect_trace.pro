Function detect_trace, im_in, SMOOTH=smooth, DISPLAY=display, $
  PEAK_VALUE=peak_value, FIT_PARAMS=fit_params, STOP=stop, SKY_IND=sky_ind, $
  PEAK_WIDTH=peak_width, PEAK_INDICE=peak_indice, SUBIND=subind, SKY_DEG=sky_deg, $
  SKY_FIT=sky_fit, BINARY=binary, UNRESOLVED=unresolved, FIBERSCRAMBLER=fiberscrambler, $
  EMPIRIC_PROFILE=empiric_profile, EMP_TRACE=emp_trace, YFIT=yfit, EMP_SKY=emp_sky
  ;On doit fournir une pose A en entree
  ;SMOOTH donne la taille du filtre pour adoucir la trace. 0 = pas de smooth.
  ;DISPLAY sert a afficher les etapes. Mettre un string pour s'en servir comme
  ; titre du plot.
  
  forward_function mpfitfun, rvb_hex, robust_sigma, mpfitfun, rvb_hex, robust_sigma, weighted_median, onemoffat_sky
  
  if keyword_set(unresolved) and ~keyword_set(binary) then binary = 1
  
  ;Avoid altering input variables
  diff = im_in
  ny = (size(diff))[2]
  xind = findgen(ny)
  
  ;Vertical median to get spatial profile
  trace = median(diff,dim=1)
  
  ;Smooth trace to facilitate detection
  if keyword_set(smooth) then $
    smtrace = smooth(trace,smooth) else $
    smtrace = trace
  
  ;Create a sky-subtracted version of the trace
  gg = where(finite(xind) and finite(smtrace))
  skypar = ladfit(xind[gg],smtrace[gg])
  skysmtrace = smtrace - poly(xind,skypar)
  
  ;Determine extremal position and its value
  nshift = 3L
  imax = where(smtrace[nshift:ny-nshift-1L] ge (shift(smtrace,nshift))[nshift:ny-nshift-1L] and smtrace[nshift:ny-nshift-1L] ge (shift(smtrace,-nshift))[nshift:ny-nshift-1L])+nshift
  ;imax = where(smtrace[3:ny-4L] ge (shift(smtrace,3))[3:ny-4L] and smtrace[3:ny-4L] ge (shift(smtrace,-3))[3:ny-4L])+3L
  void = max(smtrace[imax],wmax)
  wmax = imax[wmax]
  
  ;Estimate width of peak
  lowsig_pos = where(skysmtrace lt skysmtrace[wmax]/2.)
  delpos = (min(abs(lowsig_pos - wmax))) < 8
  
  ;Estimate sky level
  zerolevel = skypar[0]
  
  ;Detect second component if /BINARY is set
  if keyword_set(binary) then begin
    
    ;If the secondary is marginally resolved, then start estimate with the same peak position and width as the primary
    if keyword_set(unresolved) then begin
      
      wmaxsec = wmax
      delpossec = delpos
      
    endif else begin
      
      ;Also detect local minima
      imin = where(smtrace[nshift:ny-nshift-1L] le (shift(smtrace,nshift))[nshift:ny-nshift-1L] and smtrace[nshift:ny-nshift-1L] le (shift(smtrace,-nshift))[nshift:ny-nshift-1L])+nshift

      ;Reject all local maxima if no local minimum are
      ; found between it and the global maximum
      nmax = n_elements(imax)
      bad = [-1]
      for i=0L, nmax-1L do begin
        if wmax gt imax[i] then $
          good = where(imin lt wmax and imin gt imax[i], ng) else $
          good = where(imin gt wmax and imin lt imax[i], ng)
        if ng eq 0 then bad = [bad,i]
      endfor
      
      imax_sec = imax
      if n_elements(bad) ge (n_elements(imax_sec)+1) then $
        message, 'Could not find a satisfying local maximum for the companion !'
      if n_elements(bad) ne 1L then $
        remove, bad[1:*], imax_sec
      
      ;Choose the brightest local maxima in those remaining
      void = max(skysmtrace[imax_sec],/nan,wmaxsec)
      ;void = min(abs(imax_sec-wmax),wmaxsec)
      wmaxsec = imax_sec[wmaxsec]
      
      ;Estimate width of peak
      lowsig_pos = where(skysmtrace lt skysmtrace[wmaxsec]/2.)
      delpossec = (min(abs(lowsig_pos-wmaxsec))) < 8
    endelse
    
  endif
  
  ;Fit a Moffat function in the spatial profile
  ;Indices corresponding to the central peak position
  if keyword_set(binary) then subind = [2,6] else subind = [2]
  
  ;Parameters estimates
  moffat_index_estimates = 1.
  if keyword_set(binary) then $
    p0 = [zerolevel,skysmtrace[wmax],wmax,delpos,moffat_index_estimates,(skysmtrace[wmaxsec] - interpol(nmoffat_sky(xind,[0.,skysmtrace[wmax],wmax,delpos,1.]),xind,wmaxsec))>0,wmaxsec,delpossec,moffat_index_estimates] else $
    p0 = [zerolevel,skysmtrace[wmax],wmax,delpos,moffat_index_estimates]
  if keyword_set(unresolved) then p0[5] = p0[1]/2.
  if keyword_set(sky_deg) then p0 = [p0, skypar[1]]
  if n_elements(sky_deg) gt 1L then p0 = [p0, fltarr(n_elements(sky_deg)-1L)]
  
  ;Create parameters limits and start values
  parinfo = replicate({VALUE:!values.f_nan,FIXED:0B,LIMITED:[0B,0B],LIMITS:[0.,0.]},n_elements(p0))
  parinfo.VALUE = p0
  ;Accept only positive trace amplitudes
  parinfo[subind-1L].LIMITED = [1B,0B]
  parinfo[subind-1L].LIMITS = [0.,0.]
  ;Limit X position to 2* the estimated width
  parinfo[subind].LIMITED = [1B,1B]
  parinfo[subind[0]].LIMITS = parinfo[subind[0]].VALUE + [-1,1]*2*parinfo[subind[0]+1L].VALUE
  if keyword_set(binary) then $
    parinfo[subind[1]].LIMITS = parinfo[subind[1]].VALUE + [-1,1]*2*parinfo[subind[1]+1L].VALUE
  ;Limit LSF width from 1 to 12 pixels
  parinfo[subind+1L].LIMITED = [1B,1B]
  parinfo[subind+1L].LIMITS = [1.,12.]
  ;Limit Moffat indices to sane values (high values are degenerate with width, values lower than .5 are ill-defined)
  parinfo[subind+2L].LIMITED = bytarr(2,n_elements(subind))+1B
  parinfo[subind+2L].LIMITS[0] = 0.51
  parinfo[subind+2L].LIMITS[1] = 4.5
  ;Extend acceptable width if the width estimate is very large
  gg = where(parinfo[subind+1].VALUE gt (parinfo[subind+1L].LIMITS)[1], ngg)
  if ngg ne 0L then $
    for j=0L, ngg-1L do $
      parinfo[subind[gg[j]]+1L].LIMITS = [(parinfo[subind[gg[j]]+1L].LIMITS)[0],parinfo[subind[gg[j]]+1].VALUE*1.1]
  ;Apply the fit
  gg = where(finite(xind) and finite(trace))
  par = mpfitfun('nmoffat_sky',xind[gg],trace[gg],0., WEIGHTS=1D, YFIT=yfit, /QUIET, PARINFO=parinfo, STATUS=status, ERRMSG=errmsg)
  if keyword_set(sky_deg) then $
    SKY_FIT = par[(5L+fix(keyword_set(binary))*4L):*]
  
  ;Component A is ALWAYS that on the left
  if keyword_set(binary) then $
    if par[subind[0]] gt par[subind[1]] then begin
      par0 = par
      par[1:4] = par0[5:8]
      par[5:8] = par0[1:4]
    endif
  
  ;Create an empiric profile if needed
  if keyword_set(empiric_profile) then begin
    ;Create an empiric trace profile
    emp_trace = trace
    ;Subtract sky
    emp_sky = poly(xind,skypar)
    emp_trace -= emp_sky
    
    ;Kill NaNs in the profile
    bad = where(~finite(emp_trace), nbad)
    if nbad ne 0L then emp_trace[bad] = (median(emp_trace,5L))[bad]
    bad = where(~finite(emp_trace), nbad)
    if nbad ne 0L then emp_trace[bad] = 0.
    
    ;Empiric estimate of Moffat width form fitting
    moff_width = par[subind[0L]+1L]/sqrt(2*par[subind[0]+2L]+2.266065)
    
    ;Use Moffat fit to put everything at ± N sigma to zero
    lims_profile = round(par[subind[0]]+[-1,1]*10*(moff_width>par[subind[0]+1]))
    bad = where(xind lt lims_profile[0] or xind gt lims_profile[1], nbad)
    if nbad ne 0L then emp_trace[bad] = 0.
    
;    ;Smooth the profile slightly
;    emp_trace = smooth(emp_trace,2)
    
    ;Normalize profile
    emp_trace /= total(emp_trace,/nan,/double)
    
    ;Put row position to center of mass
    par[subind[0]] = total(xind*emp_trace^2,/nan)/total(emp_trace^2,/nan)
  endif
  
  ;If Fiber Scrambler is present, fit a top-hat function instead
  if keyword_set(fiberscrambler) and ~keyword_set(empiric_profile) then begin
    parinfos = [replicate({VALUE:!values.f_nan,FIXED:0B,LIMITED:[0B,0B],LIMITS:[0.,0.]},2L), parinfo]
    subind += 2L
    parinfos[0].VALUE = 13.
    parinfos[0].LIMITED = [1B,1B]
    parinfos[0].LIMITS = [0.,30.]
    parinfos[1].VALUE = 3.
    parinfos[1].LIMITED = [1B,1B]
    parinfos[1].LIMITS = [1.,5.]
    parinfos[2:*].VALUE = par
    parinfos[subind-1L].VALUE = max(smtrace,/nan)
    parinfos[subind+1L].LIMITS[0] = 1d-5
    parinfos[subind+1L].VALUE /= 20.
    parinfos[subind+2L].VALUE = 1.
    par = mpfitfun('nmoffat_sky_scrambler',xind,trace,0., WEIGHTS=1D, YFIT=yfit, /QUIET, PARINFO=parinfos, STATUS=status, ERRMSG=errmsg)
    if errmsg ne '' then stop
    ;plot,trace
    ;oplot,nmoffat_sky_scrambler(xind,parinfos.value),color=255
    ;oplot,nmoffat_sky_scrambler(xind,pars),color=rvb_hex(100,255,100)
  endif
  
  ;Display results
  if keyword_set(display) then begin
    if isa(display,'string') then title = display
    plot, trace, TITLE=title, yrange=yrange
    oplot, imax, trace[imax], color=255, /ps
    if keyword_set(empiric_profile) then $
      oplot, xind, emp_trace/max(emp_trace,/nan)*max(trace,/nan), color=rvb_hex(100,255,100) else $
      oplot, xind, yfit, color=rvb_hex(100,255,100)
    oplot, [wmax], [p0[1]+poly(wmax,skypar)], color=255, psym=4,symsize=3, thick=3
    if keyword_set(binary) then $
      oplot, [wmaxsec], [p0[5]+poly(wmax,skypar)], color=rvb_hex(100,100,255), psym=4,symsize=3, thick=3
  endif
  if keyword_set(stop) then stop
  
  ;Find sky regions
  stats = abs(trace/robust_sigma(trace,/nan))
  wsky = where(stats lt weighted_median(stats,medval=.3), nwsky)
  ;Max. 25 points de ciel
  nwsky <= 25L
  if nwsky ne 0 then begin
    ;On choisit les regions de ciel autour des traces
    wsky = wsky[sort(abs(par[subind[0]]-wsky))]
    sky_ind = wsky[0L:(nwsky-1L)]
  endif
  
  ;On affiche la region choisie pour le ciel
  if nwsky ne 0 and keyword_set(display) then $
    oplot, xind[sky_ind], trace[sky_ind], color=rvb_hex(100,100,255), /ps
  
  ;On renvoie les positions des pics
  fit_params = par
  peak_value = par[subind-1L]
  peak_width = par[subind+1L]
  peak_indice = par[subind+2L]
  return, par[subind]
  
End