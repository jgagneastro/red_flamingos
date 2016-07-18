Function optextract_1d, datav, ABSTHRESH=absthresh, COEFFV=coeffv, NSIGMAS=nsigmas, FUNC=func, MASKV=maskv, $
  PROFILE=profile, GAIN=gain, THRESH=thresh, READNOISE=readnoise, ERROR=error, VTHRESH=vthresh, OLD=old, $
  SKY=sky
  
  ;Method of Horne et al. (1986) ; http://adsabs.harvard.edu/abs/1986PASP...98..609H
  ;datav is a 1D slice of the image containing the data
  ;skyvarv is a 1D slice of the image containing the variance of the sky
  ;multv is a 1D slice of the image containing the profile
  ;maskv is a 1D slice of the image containing the mask
  
  forward_function extractfunc
  
  nx = (size(datav))[1]           ; size of vector
  xvals = dindgen(nx)
  if not keyword_set(thresh) then thresh = 5.
  if not keyword_set(gain) then gain = 1.
  if not keyword_set(readnoise) then readnoise = 0.
  if ~keyword_set(func) then func = 'extractfunc'
  if ~keyword_set(maskv) then maskv = dblarr(n_elements(datav))+1.
  crv = dblarr(n_elements(datav)) ; row of residuals
  multv = profile
  
  ; set the error threshold to be the greatest of 6 pixel, 10% of the
  ; total pixels, or the given percentage good of pixels passed in
  q = gain
  varv = abs(datav) / q + readnoise^2
  if keyword_set(sky) then $
    skyvarv = sky else $
    skyvarv = replicate(readnoise^2,n_elements(datav))
  
  if ~keyword_set(vthresh) then $
    vthresh= 20.
  coeffv  = 0                     ; initialize coeffv
  funcdone = 1                    ; initialize for while loop
  funccount = 0                   ; number of loops
  
  ; MAIN LOOP
  ; on each iteration first check to make sure there is enough good
  ; pixels left to fit the data.  Next pull out the good pixels and pass
  ; them to the function. Calculate the residuals of each pixel, using
  ; either the difference between actual and estimated, or if sigma
  ; rejection is used square the difference and divide by the variance
  ; of the pixel. If the user requests a summary plot, plot the
  ; estimated versus actual and the residual.  If needed, update the
  ; variance to reflect the new estimation. Reject the pixel with the
  ; largest residual larger then the threshold.  If no bad pixels are
  ; found, exit the loop.
  
  while (funcdone eq 1) do begin
    funcdone = 0
    goodvals = where(maskv[xvals] eq 1)
    if (n_elements(goodvals) lt thresh) then begin
      fiteval = call_function(func, dindgen(nx), datav, varv, multv*maskv, 1, coeffv, variance)
      error = sqrt(variance)
      nsigmas = sqrt(crv)
      return, fiteval
    endif
    fitx = xvals[goodvals]        ; xvals for good pixels
    fitdata = datav[fitx]         ; data for good pixels
    fitvar = varv[fitx]           ; variance of good pixels
    fitmult = multv[fitx]         ; multiplier for good pixels
    
    if ~keyword_set(old) then begin
      multv >= 0d0
      multv /= total(multv, /nan, /double)
    endif
    
    est = call_function(func, fitx, fitdata, fitvar, fitmult, 0, coeffv)
    
    If keyword_set(absthresh) then begin
      ;crv[fitx] =  abs(fitdata / fitmult - est) / abs(est)
      crv[fitx] =  abs(fitdata - fitmult * est)
    endif else begin
      crv[fitx] = (fitdata - fitmult * est)^2 / (fitvar > 1.e-6) ; no divide by 0
    endelse
    
    badpix = where(crv[fitx] gt vthresh) ; get bad locations
    if (badpix[0] ne -1) then begin
      badx = fitx[badpix]
      maxpos = where(crv[badx] eq  max(crv[badx])) ; only elimate max pixel
      maxx = badx[maxpos]
      funccount = funccount + n_elements(maxx) ; add count for bad pix
      maskv[maxx] = 0b             ; mask bad pix
      funcdone = 1                 ; set so sequence loops
    endif
    varv[fitx] = abs(fitmult * est ) / Q + skyvarv[fitx]
  endwhile
  nsigmas = sqrt(crv)
  
  if ~keyword_set(old) then begin
    multv >= 0d0
    multv /= total(multv, /nan, /double)
  endif
  
  fiteval = call_function(func, dindgen(nx), datav, varv, multv*maskv, 1, coeffv, variance)
  varv = abs(multv * fiteval ) / Q + skyvarv  ;get var for all pixels
  error = sqrt(variance)
  
  return, fiteval
  
End