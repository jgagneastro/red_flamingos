Function cshell_fringing_1D, x, A
  ; Fits a fringing solution on a 2D cshell image
  ; A[0] : Amplitude
  ; A[1] : Period
  ; A[2] : Phase
  ; A[3] : Constant
  ; A[4] : Slope
  ; A[5] : Power
  if n_elements(A) ge 6L then begin
    retval = (A[0])*((sin(x*A[1]+A[2])+1.)^A[5]-1.)+A[3]
  endif else begin
    retval = A[0]*sin(x*A[1]+A[2])+A[3]
  endelse
  if n_elements(A) ge 5L then $
    retval += A[4]*x
  return, retval
End