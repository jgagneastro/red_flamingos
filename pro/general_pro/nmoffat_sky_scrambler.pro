Function nmoffat_sky_scrambler, x, A
  ; Fits N Moffat functions + sky
  ; A[0] : Width of the top-hat region
  ; A[1] : Smooth Factor
  ; A[2] : Zero point (vertical shift)
  ; A[3] : Amplitude
  ; A[4] : Central position
  ; A[5] : Peak width
  ; A[6] : Moffat Index
  ; A[7] : 2nd peak Amplitude
  ; A[8] : 2nd peak Central position
  ; A[9] : 2nd peak Peak width
  ; A[10] : 2nd peak Moffat Index
  ; (...)
  ; A[n] : 1st polynomial degree for sky fitting
  ; A[n+1] : 2nd polynomial degree for sky fitting
  ; etc.
  forward_function interpol ,frac
  ;nmof = floor(float(n_elements(A)-3L)/4.)
  nmof = ceil(float(n_elements(A)-3L)/4.)
  if nmof ne 1L then message, 'This code is not meant for binary traces !'
  nsky = n_elements(A)-2L-(nmof*4L+1L)
  ;Start with the sky estimate
  val = A[2]
  ;Add a gmos profiles
  for i=0L, nmof-1L do $
    val += A[4*i+3L]/(((x-A[4*i+4L])/A[4*i+5L])^2 + 1)^A[4*i+6L]
  
  if finite(A[0]) and A[0] gt 0. then begin
    pos = A[4]
    flux = A[3]
    nsmooth = A[1]
    plateau = A[0]
    valint = interpol(val,x,x+frac(pos))
    valint_left = valint[0:floor(pos)]
    valint_right = valint[floor(pos):*]
    x_left =  x[0:floor(pos)]
    x_right =  x[floor(pos):*]
    midpl = replicate(flux,floor(plateau)+1L)
    full_profile = [interpol(valint_left,x_left,x_left-(1.-frac(plateau))/2.),midpl,interpol(valint_right,x_right,x_right+(1.-frac(plateau))/2.)]
    full_x = [x, (x[-1]-x[-2])*(1.+findgen(n_elements(midpl)+1L))+x[-1]]
    val = interpol(full_profile,full_x,x+(double(floor(plateau)+1L))/2.)
    if A[1] ge 1L and round(A[1]) lt (ceil(pos-plateau/2.*1.2)-ceil(pos-plateau/2.*1.2)+1L) then $
      val[floor(pos-plateau/2.*1.2):ceil(pos+plateau/2.*1.2)] = smooth(val[floor(pos-plateau/2.*1.2):ceil(pos+plateau/2.*1.2)],round(A[1]))
  endif
  
  isky0 = n_elements(a)-nsky
  ;Add an N-order polynomial
  if nsky gt 0L then $
    val += total(A[isky0:*]#(dblarr(n_elements(x))+1D) * ((dblarr(nsky)+1D)#x)^((dindgen(nsky)+1D)#(dblarr(n_elements(x))+1D)), 1, /nan)
  return, val
End