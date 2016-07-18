Function best_lamp_shift, lambda_ref, lamp_ref, lamp_obs, poly_estim, YFIT=yfit, XFIT=xfit, XVEC=xvec, MONTECARLO=montecarlo, $
  CONSTELLATIONS=constellations, METRIC=metric, LOW_CLIP=low_clip, PEAK_CHI=peak_chi, BAD_REDUCE=bad_reduce, REFINE=refine, INDFIT=indfit, POLYPAR=polypar
  
  forward_function killnan, weighted_median, build_calib_constellation, supersmooth, find_1d_peaks, mpfitfun, build_calib_spectrum, refine_wavelength_solution_relative
  common lampshiftamoeba, pixels, lam, sp, weights, norm
  
  ;Import variables
  lam = killnan(lambda_ref,0.)
  sp = killnan(lamp_ref,0.)
  spo = killnan(lamp_obs,0.)
  p0 = poly_estim
  
  if ~keyword_set(low_clip) then $
    low_clip = .01
  
  ;Ensure that both spectra have their base at zero
  sp -= weighted_median(sp,medval=.1)
  spo -= weighted_median(spo,medval=.1)
  
  ;Ensure that both spectra are normalized the same way
  sp /= max(sp,/nan)
  spo /= max(sp,/nan)
  ;sp /= weighted_median(sp,medval=.99)
  ;spo /= weighted_median(spo,medval=.99)
  
  if keyword_set(constellations) then begin
    
    if ~keyword_set(metric) then metric = [1.,sqrt(3.)]
    metric = [0.,1.]
    
    ;Compute "constellation coordinates" for each pair of 3 lines, which are independent upon linear transformations
    obs_constellations = build_calib_constellation(lindgen(n_elements(spo)), spo, nsec=10, nuse=50, SIG_THRESHOLD=3., FLUX_THRESHOLD=.01)
    ref_constellations = build_calib_constellation(lambda_ref, sp, nsec=10, nuse=50, SIG_THRESHOLD=3., FLUX_THRESHOLD=.01)
    
    obs_constellations = obs_constellations[sort(obs_constellations[*,0]),*]
    ref_constellations = ref_constellations[sort(ref_constellations[*,0]),*]
    
    ;Compute the cross-match probabilityn density function (PDF) for each combination of constellations
    nobs_constellations = (size(obs_constellations))[1]
    nref_constellations = (size(ref_constellations))[1]
    pdf_matrix = fltarr(nobs_constellations,nref_constellations)+!values.f_nan
    for i=0L, nobs_constellations-1L do $
      for j=0L, nref_constellations-1L do $
        pdf_matrix[i,j] = -total((reform(obs_constellations[i,1:4])-reform(ref_constellations[j,1:4]))^2*metric[[0,0,0,1]]^2)
    
    void=max(exp(pdf_matrix/.01),dim=2,wmax)
    wmax2 = array_indices(pdf_matrix,wmax)
    plot,wmax2[1,*],/ps
    ;vf,exp(pdf_matrix/.01)>.95
    stop
    ;Doesn't work very well ..!
  endif
  
  ;Define abscissa vector
  nx = n_elements(spo)
  if keyword_set(xvec) then $
    x = xvec else $
    x = lindgen(nx)
  
  if ~keyword_set(bad_reduce) then $
    bad_reduce = 10.
  
  ;Use a tapering of the edges
  taper_array = fltarr(nx)+1.
  taper_wid = long(double(nx)/35d0)
  if taper_wid ge 1L then begin
    taper_array[0L:taper_wid-1L] = findgen(taper_wid)/float(taper_wid-1L)
    taper_array[-taper_wid:*] = reverse(findgen(taper_wid)/float(taper_wid-1L))
    ;Smooth the shit out of it
    taper_array = supersmooth(supersmooth(taper_array,taper_wid),taper_wid)
    bad = where(spo lt low_clip, nbad)
    if nbad ne 0 then taper_array[bad] /= bad_reduce
    bad = where(spo lt 0., nbad)
    if nbad ne 0 then begin
      spo[bad] = 0
      taper_array[bad] = 0
    endif
  endif
  
  ;Define arguments that will be passed to the fitting subroutine 
  norm = 1
  functargs = {lam:lam,sp:sp,NORM:norm}
  bounds = .97
  if keyword_set(montecarlo) then begin
    if montecarlo eq 1 then montecarlo = 40
    nMonte = montecarlo
    nMontej = nMonte*2L
    factorsi = 10.^((findgen(nMonte)/float(nMonte-1L)-.5)/60.)
    factorsj = 10.^((findgen(nMontej)/float(nMontej-1L)-.5)/30.)
    factorsi[0] = 1.
    factorsj[0] = 1.
    redchis = fltarr(nMonte,nMontej)+!values.f_nan
    redpeaks = fltarr(nMonte,nMontej)+!values.f_nan
    peaks_o = find_1d_peaks(spo, FLUX_THRESHOLD=.01, SIG_THRESHOLD=3.)
    pars = fltarr(nMonte,nMontej,n_elements(poly_estim))+!values.f_nan
    weights = taper_array
    for i=0L, nMonte-1L do begin
      for j=0L, nMontej-1L do begin
        parinfo = replicate({value:0d0, fixed:0, limited:[1,1], limits:[0d0,0], mpmaxstep:0, step:0}, n_elements(poly_estim)+1)
        parinfo.VALUE = [poly_estim*[factorsj[j],factorsi[i],fltarr(n_elements(poly_estim)-2L)+1.],1.]
        ;Limits are a multiplicative factor 1/bound or bound; however make sure that they're in increasing order
        parinfo.LIMITS = (([bounds,1./bounds]#(fltarr(n_elements(parinfo))+1.))^((fltarr(2)+1.)#sign(parinfo.VALUE)))*((fltarr(2)+1.)#(parinfo.VALUE))
        parinfo[2:-2].FIXED = 1
        parinfo[-1].LIMITS=[0.05,20.]
        if n_elements(poly_estim) gt 3L then begin
          parinfo.MPMAXSTEP = [.1, 1d-4, 1d-11, replicate(1d-11,n_elements(poly_estim)-3L), .3]
          parinfo.STEP = [.1, 1d-4, 1d-11, replicate(1d-11,n_elements(poly_estim)-3L), .3]/10.
        endif else begin
          parinfo.MPMAXSTEP = [.1, 1d-4, 1d-11, .3]
          parinfo.STEP = [.1, 1d-4, 1d-11, .3]/10.
        endelse
        par = mpfitfun('build_calib_spectrum', x, spo, PARINFO=parinfo, WEIGHTS=weights, YFIT=yfit, STATUS=status, ERRMSG=errmsg, FUNCTARGS=functargs, /QUIET)
        ;par = mpfitfun('build_calib_spectrum', x, spo, 0., poly_estim*[1.,factorsi[i],factorsj[j]], WEIGHTS=taper_array, YFIT=yfit, STATUS=status, ERRMSG=errmsg, FUNCTARGS=functargs, /QUIET)
        if errmsg ne '' then continue
        spmod = build_calib_spectrum(x,par,lam=lam,sp=sp,NORM=functargs.NORM)
        veci = (spo-spmod)^2
        ggi = lindgen(n_elements(veci))
        ;ggi = where(veci ge .1 * max(veci,/nan), nggi)
        redchis[i,j] = mean(veci[ggi]*taper_array[ggi],/nan)
        pars[i,j,*] = par[0:-2]
        if keyword_set(peak_chi) then begin
          peaks_ij = find_1d_peaks(spmod, FLUX_THRESHOLD=.01, SIG_THRESHOLD=3.)
          peakdist_2d = (peaks_o#(fltarr(n_elements(peaks_ij))+1.)-(fltarr(n_elements(peaks_o))+1.)#peaks_ij)^2
          redpeaks[i,j] = total(min(peakdist_2d,dim=2,/nan))
        endif
      endfor
    endfor
    
    if keyword_set(peak_chi) then begin
      if total(finite(redpeaks)) eq 0L then message, ' All reduced chi^2 are NaNs !'
      void = min(redpeaks,/nan,wmin)
      wmin2 = array_indices(redpeaks,wmin)
    endif else begin
      if total(finite(redchis)) eq 0L then message, ' All reduced chi^2 are NaNs !'
      void = min(redchis,/nan,wmin)
      wmin2 = array_indices(redchis,wmin)
    endelse
    par = reform(pars[wmin2[0L],wmin2[1L],*])
    ;for i=0L, (size(pars))[1]-1L do begin & for j=0L, (size(pars))[2]-1L do begin & plot,lamp_obs^indfit & oplot,build_calib_spectrum(x,reform(pars[i,j,*]),NORM=functargs.NORM,lam=lam,sp=sp)^indfit,color=255 & wait, 0.1 & endfor & endfor
  endif else begin
    parinfo = replicate({value:0d0, fixed:0, limited:[1,1], limits:[0d0,0]}, n_elements(poly_estim))
    parinfo.VALUE = poly_estim
    parinfo[2].FIXED = 1
    ;Limits are a multiplicative factor 1/bound or bound; however make sure that they're in increasing order
    parinfo.LIMITS = (([bounds,1./bounds]#(fltarr(n_elements(parinfo))+1.))^((fltarr(2)+1.)#sign(parinfo.VALUE)))*((fltarr(2)+1.)#(parinfo.VALUE))
    par = mpfitfun('build_calib_spectrum', x, spo, PARINFO=parinfo, WEIGHTS=taper_array, YFIT=yfit, STATUS=status, ERRMSG=errmsg, FUNCTARGS=functargs, /QUIET,NORM=functargs.NORM)
    ;par = mpfitfun('build_calib_spectrum', x, spo, 0., poly_estim, WEIGHTS=taper_array, YFIT=yfit, STATUS=status, ERRMSG=errmsg, FUNCTARGS=functargs, /QUIET)
    if errmsg ne '' then stop
  endelse
  
  if keyword_set(refine) then begin
    lamp_ref2 = build_calib_spectrum(x,par,lam=lam,sp=sp)^indfit
    lamp_obs2 = lamp_obs^indfit
    
    lamp_ref2 /= max(lamp_ref2,/nan)
    lamp_obs2 /= max(lamp_obs2,/nan)
    flux_threshold = .03
    sig_threshold = 5.
    peaks_ref = find_1d_peaks(lamp_ref2,flux_threshold=flux_threshold,sig_threshold=sig_threshold)
    bad = where(lamp_ref2[peaks_ref] lt flux_threshold, nbad)
    if nbad ne 0. then remove, bad, peaks_ref
    ;plot,lamp_ref2 & oplot,peaks_ref,lamp_ref2[peaks_ref],/ps,color=255,symsize=2
    
;    ;Start with a crude cross-match
;    nlag = 51L
;    lag = (findgen(nlag) - float(nlag-1L)/2.)
;    cc = c_correlate(killnan(lamp_obs2),killnan(lamp_ref2),lag)
;    ;On determine le point où la corrélation est maximale
;    maxi = max(cc,ii,/nan)
;    stop
    lag=0 & ii=0
    
    refinex = refine_wavelength_solution_relative(x,1-lamp_obs2,shift(x,lag[ii]),1-lamp_ref2,use_lines=x[peaks_ref],MINAMP=flux_threshold,/slope,wid=12,minwid=0.01,maxwid=2.,$
      dmaxfit=.4,POLYPAR=polypar,MAXDXREF=1.,/nodisplay)
    refinex = shift(refinex,lag[ii])
    
    ;plot,x,poly(x,par) & oplot,refinex,poly(x,par),color=255
    ;plot,x,poly(x,par)-x & oplot,refinex,poly(x,par)-x,color=255
  endif
  
  xfit = x
  return, par 
End