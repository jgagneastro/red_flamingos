Pro sortself, vect, UNIQUE=unique
  ;Permet de classer directement un vecteur en ordre croissant.
  ;/unique permet d'enlever les repetitions en meme temps
  compile_opt hidden
  on_error, 2
  forward_function uniq
  vect = vect[sort(vect)]
  if keyword_set(unique) then vect = vect[uniq(vect)]
End