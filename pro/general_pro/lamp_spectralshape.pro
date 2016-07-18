Function lamp_spectralshape, flat_in, NWID=nwid, NOPATCH=nopatch
  ;Correct for the spectral distribution of a lamp
  ;Now we will correct flat's illumination for the spectral signature of the lamp
  ;NWID : Width of the median filter used for smoothing the spectral distribution of the lamp
  
  ;Determine the X dimension of each trace
  flat = flat_in
  nx = (size(flat))[1]
  ny = (size(flat))[2]
  if ~keyword_set(nwid) then nwid = 64L
  
  ;Find the spectral distribution of the lamp
  med = median(flat,dim=2)
  med = smooth(median(med,nwid),nwid)
  
  ;Remove the spectral structure of the lamp
  flat /= (med#make_array(ny,value=1.,/float))
  
  ;Remove bad pixels
  if ~keyword_set(nopatch) then begin
    bad = where(flat le 0.005, nbad)
    if nbad ne 0 then flat[bad] = !values.d_nan
  endif
  
  return, flat
  
End