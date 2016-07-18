Function cshell_fringing_1D2, x, A
  ; Fits a fringing solution on a 2D cshell image
  ; A[0] : Amplitude
  ; A[1] : Finesse
  ; A[2] : Period
  ; A[3] : Phase
  ; A[4] : Zero-level
  ; A[5] : Slope
  
  sinf = sin(x*A[2]+A[3])
  hrange = (1.-1./(1.+A[1]))/2.
  interf = (1./(1.+A[1]*sinf^2)-1+hrange)/hrange
  retval = A[0]*interf+A[4]
  if n_elements(a) gt 5L then $
    retval += A[5]*x
  return, retval
End