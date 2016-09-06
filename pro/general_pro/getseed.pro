Function getseed
;  Cette fonction sert a obtenir un seed quelconque base sur le temps du systime, qui change tres rapidement
;  en fonction du temps.
  compile_opt hidden ;Evite de toujours afficher 'Compiled module'
  x = double(abs(float(strmid(systime(),11,2)+strmid(systime(),14,2)+strmid(systime(),17,2)))*(systime(/seconds)-1216243200d0))
  return, double(abs(x*cos(sin(x)*sqrt(x)))/100.)/1d9
End