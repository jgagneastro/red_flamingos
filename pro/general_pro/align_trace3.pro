Pro align_trace3, im_in, calA, calB, calC, NSIGMA=nsigma, NDEG=ndeg, SMOOTH=smooth, BORDERCUT=bordercut, SMPTS=smpts, NITER=niter, STOP=stop, $
  MINPIXSHIFT=minpixshift, NP=np, SIZELAG=sizelag, HORIZONTAL_SMOOTH=horizontal_smooth, MEDIANSMOOTH=mediansmooth, BOX_SIZE=box_size, DEBUG=debug

  forward_function reverse, remove, robust_sigma, c_correlate, poly_fit, interpol2, fc_correlate, weighted_median

  ;Fait pour un seul CCD a la fois
  if ~keyword_set(nsigma) then nsigma = 3.
  if minpixshift eq !NULL then minpixshift = 0.
  ;BORDERCUT : Le % qu'on coupe de chaque bord
  if ~keyword_set(bordercut) then bordercut = .01
  if ~keyword_set(niter) then niter = 2L
  if ~keyword_set(ndeg) then ndeg = 2L
  if np eq !NULL then np = 5L
  if ~keyword_set(sizelag) then $
    sizelag = 120.
  if keyword_set(stop) then stop
  if horizontal_smooth eq !NULL then horizontal_smooth = 10L 
  if mediansmooth eq !NULL then mediansmooth = 5L
  if box_size eq !NULL then box_size = 30L
  
  im = im_in
  nx = (size(im))[1]
  nanpos = where(~finite(im), nnanpos)
  
  ;Remove vertical structures
  ny = (size(im))[2]
  im -= median(im,dim=2)#make_array(ny,value=1d0,/double)
  
  ;Do some horizontal_smoothing
  if keyword_set(horizontal_smooth) then $
    for i=0L, ny-1L do $
      im[*,i] = median(im[*,i],horizontal_smooth)
  if keyword_set(smooth) then im = smooth(median(im,smooth),smooth)
  if keyword_set(mediansmooth) then im = median(im,mediansmooth)
  
  imder = (shift(im,0,0)-shift(im,0,-1))/2.
  if keyword_set(horizontal_smooth) then $
    for i=0L, ny-1L do $
      imder[*,i] = median(imder[*,i],horizontal_smooth)
  if keyword_set(smooth) then imder = smooth(median(imder,smooth),smooth)
  if keyword_set(mediansmooth) then imder = median(imder,mediansmooth)
  
  if keyword_set(debug) then begin
    vf, imder
    stop
  endif
  
  im0 = im
  im = imder
  
  ver_smooth_trace = im0
  medmax = fltarr(nx)+!values.d_nan
  medmin = fltarr(nx)+!values.d_nan
  for i=0L, nx-1L do begin & $
    ver_smooth_trace[i,*] = median(im0[(i-box_size/2)>0:(i+box_size/2)<(nx-1),*],dim=1) & $
    medmax[i] = max(ver_smooth_trace[i,*],/nan) & $
    medmin[i] = max(-ver_smooth_trace[i,*],/nan) & $
  endfor
  intensity_profile = smooth((median(medmax,10)<median(medmin,10)),box_size,/nan)
  max_intensity = max(intensity_profile,/nan,wmax_intensity)
  
  for ite=0, niter-1L do begin
    
    sz = size(im)
    ;trace = median(im,dim=1)
    ;trace -= median(trace)
    ;trace >= 0
    trace = total(im[(wmax_intensity-box_size/2)>0:(wmax_intensity+box_size/2)<(nx-1),*],1,/nan)
    trace -= median(trace)
    trace /= weighted_median(abs(trace),medval=.97)
    
    pos = fltarr(nx) + !values.f_nan
    lag = findgen(sizelag)+1.
    lag = [-reverse(lag),0.,lag]
    bad = where(abs(lag) gt (size(im_in))[1], nbad)
    if nbad ne 0 then remove, bad, lag

    x0 = findgen(n_elements(lag))-float(n_elements(lag))/2.

    ngood = 0L
    ccs_value = dblarr(nx)+!values.d_nan
    for k=0L,nx-1L do begin
      medim = median(im[(k-box_size/2)>0:(k+box_size/2)<(nx-1),*],dim=1)
      if abs(max(abs(medim),/nan)/median(abs(medim)-min(abs(medim),/nan))) le nsigma then continue
      if abs(max(abs(im[k,*]),/nan))/median(abs(im[k,*])-min(abs(im[k,*]),/nan)) le nsigma then continue
      if abs(max(abs(medim),/nan)/robust_sigma(medim,/nan)) le nsigma then continue

      medim -= median(medim)
      medim /= weighted_median(abs(medim),medval=.97)
      ;medim >= 0

      ccs = fc_correlate(medim,trace,lag, /INTEGERS)
      ;ccs = c_correlate(medim,trace,lag)
      ccs *= exp(-x0^2/(2.*sizelag^2)*3)

      ;On fait un fit polynomial de 2e degré dans la fonction de corrélation autour du maximum
      ccs_value[k] = max(ccs,ii)
      if ii gt n_elements(lag)-np-1L or ii lt np then continue
      ngood += 1
      fit = poly_fit(lag[ii-np:ii+np],ccs[ii-np:ii+np],2.)
      ;fit = poly_fit(lag,ccs,2.)

      pos[k] = -.5*fit[1]/fit[2]
    endfor

    if keyword_set(stop) then stop
    xx = findgen(nx)
    yy = pos
    g = where(finite(yy),ng)
    if ng eq 0 then continue

    if bordercut gt 0 then begin
      if ng lt 2*long(bordercut*ng)+10L then message, ' Could not straighten the trace !'
      if long(bordercut*ng) ge 1 then begin
        remove, ng-lindgen(long(bordercut*ng))-1L, g
        remove, lindgen(long(bordercut*ng)), g
      endif
      ng = n_elements(g)
    endif

    if keyword_set(smpts) then $
      yy[g] = smooth(median(yy[g],smpts),smpts)

    if keyword_set(stop) then stop

    if max(abs(yy[g]),/nan) lt minpixshift then break

    ;fit = poly_fit(xx[g],yy[g],ndeg)
    nsmooth_ccs = 100L
    ccs_value_relative = ccs_value/smooth(ccs_value,nsmooth_ccs,/nan)
    ccs_value_relative[0:nsmooth_ccs/2L] = 1d-12
    ccs_value_relative[(-(nsmooth_ccs/2L)-1L):-1L] = 1d-12
    fit = poly_fit(xx[g],yy[g],ndeg,measure_errors=1./ccs_value_relative[g]^2)
    
    if keyword_set(debug) then begin
      plot,xx,poly(xx,fit) & oplot,xx[g],yy[g],color=255,/ps
      stop
      plot,xx[g],yy[g],/ps & oplot,xx,poly(xx,fit), color=255
      stop
    endif
    
    ;Repeat with a robust fit
    bad = where( (yy[g]-poly(xx[g],fit))^2 gt 2., nbad)
    if nbad ne 0L then remove, bad, g
    ng = n_elements(g)
    fit = poly_fit(xx[g],yy[g],ndeg,measure_errors=1./ccs_value_relative[g]^2)
    
    if keyword_set(debug) then begin
      plot,xx,poly(xx,fit) & oplot,xx[g],yy[g],color=255,/ps
      stop
      plot,xx[g],yy[g],/ps & oplot,xx,poly(xx,fit), color=255
      stop
    endif
    
    if keyword_set(stop) then stop
    polyr = poly(lindgen(nx),fit)
    polyr -= median(polyr)
    for k=0L,nx-1L do begin

      if keyword_set(calA) then $
        calA[k,*] = interpol2(calA[k,*],findgen(sz[2]),findgen(sz[2])-polyr[k], badvalue=!values.d_nan)
      if keyword_set(calB) then $
        calB[k,*] = interpol2(calB[k,*],findgen(sz[2]),findgen(sz[2])-polyr[k], badvalue=!values.d_nan)
      if keyword_set(calC) then $
        calC[k,*] = interpol2(calC[k,*],findgen(sz[2]),findgen(sz[2])-polyr[k], badvalue=!values.d_nan)
      im[k,*] = interpol2(im[k,*],findgen(sz[2]),findgen(sz[2])-polyr[k], badvalue=!values.d_nan)
      im_in[k,*] = interpol2(im_in[k,*],findgen(sz[2]),findgen(sz[2])-polyr[k], badvalue=!values.d_nan)

      if nnanpos ne 0L then begin
        im[nanpos] = !values.d_nan
        im_in[nanpos] = !values.d_nan
      endif

      ;      if keyword_set(calA) then $
      ;        calA[k,*] = interpolate(calA[k,*],findgen(sz[2])-polyr[k],cubic=-.6,missing=!values.f_nan)
      ;      if keyword_set(calB) then $
      ;        calB[k,*] = interpolate(calB[k,*],findgen(sz[2])-polyr[k],cubic=-.6,missing=!values.f_nan)
      ;      im[k,*] = interpolate(im[k,*],findgen(sz[2])-polyr[k],cubic=-.6,missing=!values.f_nan)
      ;      im_in[k,*] = interpolate(im_in[k,*],findgen(sz[2])-polyr[k],cubic=-.6,missing=!values.f_nan)
    endfor
    if keyword_set(stop) then stop
  endfor

End