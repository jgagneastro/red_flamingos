Function fringing_fit_1d, lambda_in, fringing_in, ESTIMATED_PERIOD=estimated_period, SINUS=sinus, YFIT=yfit, ERR=err, WEIGHTS=weights_in, FACTOR=factor
  
  forward_function interpol2, supersmooth, mpfitfun, cshell_fringing_1d, cshell_fringing_1d2, remove
  
  ;Interpolate the fringing on a regular wavelength grid
  if ~keyword_set(factor) then factor = 3.
  nlam = n_elements(lambda_in) * factor
  lambda = dindgen(nlam)/double(nlam-1L) * (max(lambda_in,/nan) - min(lambda_in,/nan)) + min(lambda_in,/nan)
  fringing = interpol2(fringing_in, lambda_in, lambda, /repairnan, badval=!values.f_nan)
  if keyword_set(weights_in) then $
    weights = interpol2(weights_in, lambda_in, lambda, badval=0.)
  smooth_fringing = supersmooth(fringing,4.*factor)
  
  ;Fit the 1D fringing with a full interference profile (with Finesse parameter)
  amp0 = (max(smooth_fringing,/nan)-min(smooth_fringing,/nan))/2.
  level0 = (max(smooth_fringing,/nan)+min(smooth_fringing,/nan))/2.
  if ~keyword_set(estimated_period) then begin
    ;Find all local maxima
    ggmax = where(smooth_fringing[1L:-2L] gt (shift(smooth_fringing,1L))[1L:-2L] and smooth_fringing[1L:-2L] gt (shift(smooth_fringing,-1L))[1L:-2L], nggmax) + 1L
    ;Find all local minima
    ggmin = where(smooth_fringing[1L:-2L] lt (shift(smooth_fringing,1L))[1L:-2L] and smooth_fringing[1L:-2L] lt (shift(smooth_fringing,-1L))[1L:-2L], nggmin) + 1L
    ;Determine median period
    estimated_period = median([ggmax[1L:*]-ggmax[0:-2L],ggmin[1L:*]-ggmin[0:-2L]]) * 2
  endif
  p0 = [amp0,0.01,(!dpi/estimated_period),0.,level0,0.]
  if keyword_set(sinus) then begin
    p0[2] *= 2.
    remove, 1, p0
  endif
  par = mpfitfun('cshell_fringing_1d'+(['2',''])[keyword_set(sinus)],lambda,fringing,1.,p0,YFIT=yfit_hires,status=status,err=err,/quiet,/nan, WEIGHTS=weights)
  par2 = par
  par2[4L-keyword_set(sinus)] = 1.
  par2[5L-keyword_set(sinus)] = 0.
  yfit_hires = cshell_fringing_1d(lambda,par2)
  yfit = interpol2(yfit_hires, lambda, lambda_in, badval=1.); / (par[4L-keyword_set(sinus)] + lambda*par[5L-keyword_set(sinus)])
  return, par
End