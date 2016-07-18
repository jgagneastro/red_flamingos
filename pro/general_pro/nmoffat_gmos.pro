Function nmoffat_gmos, x, A
; Fait fitter N gaussiennes d'un coup.
; A[0] : Niveau zero
; A[1] : Valeur maximale de la 1ere gaussienne
; A[2] : x0 de la 1ere gaussienne
; A[3] : σ de la 1e gaussienne
; A[4] : Indice de Moffat de la 1e gaussienne
; A[5] : Valeur maximale de la 2e gaussienne
; A[6] : x0 de la 2e gaussienne
; A[7] : σ de la 2e gaussienne
; A[8] : Indice de Moffat de la 2e gaussienne
; etc.
  val = A[0]
  nmof = ceil((n_elements(A)-1.)/4.)
  if n_elements(A) lt (nmof*4+1L) then message, ' Wrong number of arguments for '+strtrim(nmof,2)+' gaussian fits !'
  for i=0L, nmof-1L do $
    val += A[4*i+1L]/(((x-A[4*i+2L])/A[4*i+3L])^2 + 1)^A[4*i+4L]
  return, val
End