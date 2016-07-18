Function onemoffat_sky, x, A
  ; Fits one Moffat function + sky
  ; A[0] : Zero point (vertical shift)
  ; A[1] : Amplitude
  ; A[2] : Central position
  ; A[3] : Peak width
  ; A[4] : Moffat Index
  ; A[5] : 1st polynomial degree for sky fitting 
  ; A[6] : 2nd polynomial degree for sky fitting
  ; A[7] : (...)
  ; etc.
  val = A[0]
  nsky = n_elements(A)-5L
  if n_elements(A) lt 5L then message, ' Wrong number of arguments for '+strtrim(nmof,2)+' onemoffat_sky fits !'
  val += A[1L]/(((x-A[2L])/A[3L])^2 + 1)^A[4L]
  if nsky gt 0L then $ 
    val += total(A[5:*]#(dblarr(n_elements(x))+1D) * ((dblarr(nsky)+1D)#x)^((dindgen(nsky)+1D)#(dblarr(n_elements(x))+1D)), 1, /nan)
  ;val += total((dblarr(n_elements(x))+1D)#A[5:*] * (x#(dblarr(nsky)+1D))^((dblarr(n_elements(x))+1D)#(dindgen(nsky)+1D)), 2, /nan)
  return, val
End