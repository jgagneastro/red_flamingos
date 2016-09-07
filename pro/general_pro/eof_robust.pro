Function eof_robust, file
  compile_opt hidden
 ; return, eof(file)
  ;S'il y  a une erreur, on attend une demi-seconde et on recommence
  catch, error_status
  if error_status ne 0 then begin
    error_status = 0
    wait, 0.5
    return, eof_robust(file)
  endif
  return, eof(file)
End