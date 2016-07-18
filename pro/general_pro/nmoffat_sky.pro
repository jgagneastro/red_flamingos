Function nmoffat_sky, x, A
  ; Fits N Moffat functions + sky
  ; A[0] : Zero point (vertical shift)
  ; A[1] : Amplitude
  ; A[2] : Central position
  ; A[3] : Peak width
  ; A[4] : Moffat Index
  ; A[5] : 2nd peak Amplitude
  ; A[6] : 2nd peak Central position
  ; A[7] : 2nd peak Peak width
  ; A[8] : 2nd peak Moffat Index
  ; (...)
  ; A[n] : 1st polynomial degree for sky fitting 
  ; A[n+1] : 2nd polynomial degree for sky fitting
  ; etc.
  nmof = floor(float(n_elements(A)-1L)/4.)
  nsky = n_elements(A)-(nmof*4L+1L)
  ;Start with the sky estimate
  val = A[0]
  ;Add N gmos profiles
  for i=0L, nmof-1L do $
    val += A[4*i+1L]/(((x-A[4*i+2L])/A[4*i+3L])^2 + 1)^A[4*i+4L]
  isky0 = n_elements(a)-nsky
  ;Add an N-order polynomial
  if nsky gt 0L then $
    val += total(A[isky0:*]#(dblarr(n_elements(x))+1D) * ((dblarr(nsky)+1D)#x)^((dindgen(nsky)+1D)#(dblarr(n_elements(x))+1D)), 1, /nan)
  return, val
End