Function killnan,DATAIN,remplace
  
  compile_opt hidden
  
  DATA = DATAIN
  mauvais = where(finite(data) eq 0)
  if total(size(remplace)) eq 0 then remplace = 0
  
  if mauvais[0] eq (-1) then begin
  	return,data
  endif else begin
  	DATA[mauvais] = remplace
  	return,DATA
  endelse
  
End