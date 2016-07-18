function frac, x
  ;Cette fonction retourne la partie fractionnaire d'un nombre en précision double.
  compile_opt hidden ;Evite de toujours afficher 'Compiled module'
  pos = where(x ge 0)
  neg = where(x lt 0)
  rep=x
  if pos(0) ne -1 then rep(pos) = (x(pos)-floor(x(pos)))
  if neg(0) ne -1 then rep(neg) = (abs(x(neg))-floor(abs(x(neg))))
  return, rep
End