Function build_calib_spectrum, pixels, params, LAM=lam, SP=sp, NORM=norm
  
  forward_function interpol2, killnan
  
  ;Number of params and pixels
  npar = n_elements(params)
  if keyword_set(norm) then npar -= 1L
  npix = n_elements(pixels)
  
  ;Build a wavelength vector from the parameters
  wv_vector = total((make_array(npix,value=1.,/float)#params)*(pixels#make_array(npar,value=1.,/float))^(make_array(npix,value=1.,/float)#findgen(npar)),2,/nan)
  
  ;Interpolate the reference calibration on the pixel array
  out_sp = interpol2(sp,lam,wv_vector,badval=0.)
  out_sp = killnan(out_sp,0.)
  
  if keyword_set(norm) then out_sp *= params[-1L] else out_sp /= median(out_sp) 
  
  return, out_sp
End