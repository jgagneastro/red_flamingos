Function AB_moffat_fixed_sep, x, A
  ; Fait fitter N gaussiennes d'un coup.
  ; Tous les indices de Moffat sont fixes
  ; A[0] : Niveau zero
  ; A[1] : Position de la primaire de gauche
  ; A[2] : Separation primaire / secondaire
  ; A[3] : Taille de dither
  ; A[4] : Indice de Moffat
  ;
  ; A[5] : Valeur maximale de la 1ere gaussienne
  ; A[6] : σ de la 1e gaussienne
  ;
  ; A[7] : Ratio de flux 2e/1e gaussiennes
  ; A[8] : σ de la 2e gaussienne
  ;
  ; A[9] : Ratio de flux positif / negatif ou l'inverse (trace de droite / trace de gauche; 3e / 1e)
  ; A[10] : σ de la 3eme gaussienne
  ;
  ; A[11] : Ratio de flux 4e/3e gaussiennes
  ; A[12] : σ de la 4eme gaussienne
  ;Note : 4/2 = (4/3) / (2/1) * (3/1) = A[11] / A[7] * A[9]
  
  val = A[0] + A[5L]/(((x-A[1L])/A[6L])^2 + 1)^A[4L] + $ ;Niveau 0 + Profil A gauche
    A[7L]*A[5L]/(((x-(A[1L]+A[2L]))/A[8L])^2 + 1)^A[4L] + $;Profil B gauche
    A[9L]*A[5L]/(((x-(A[1L]+A[3L]))/A[10L])^2 + 1)^A[4L] + $;Profil A droite
    A[11L]*A[9L]*A[5L]/(((x-(A[1L]+A[2L]+A[3L]))/A[12L])^2 + 1)^A[4L];Profil B droite
  
  return, val
End