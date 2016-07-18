Function rvb_hex, RR, VV, BB, totalv,RVB=RVB
  ;Prend en entrée des valeurs R,V,B sur un total de 255 ou "totalv", puis les convertit en BBVVRR hexadécimal
  
  on_error, 2
  compile_opt hidden
  forward_function to_hex
  
  ;Couleur blanche
  if n_params() eq 1 then begin
    if size(RR,/type) eq 7 then begin
      RR = strlowcase(strtrim(RR,2))
      if RR eq 'blanc' or RR eq 'b' or RR eq 'white' or RR eq 'w' then return, 16777215L
    endif
  endif
  
  if keyword_set(rvb) then begin
    RR = rvb[0]
    VV = rvb[1]
    BB = rvb[2]
  endif
  
  if n_params(invec) lt 3 and ~keyword_set(rvb) then $
    message, 'Les trois premiers arguments doivent etre les 3 composantes RVB '
    
  if ~keyword_set(totalv) then totalv = 256. 
  
  str1 = strupcase(strmid(to_hex(fix((BB/float(totalv)*256. mod 256))),0,2))
  str2 = strupcase(strmid(to_hex(fix((VV/float(totalv)*256. mod 256))),0,2))
  str3 = strupcase(strmid(to_hex(fix((RR/float(totalv)*256. mod 256))),0,2))
  if strlen(str1) then str1 = '0'+str1
  if strlen(str2) then str2 = '0'+str2
  if strlen(str3) then str3 = '0'+str3
  void = execute('ll = '''+str1+str2+str3+'''x')
  
  return, ll
  
End