Function offset_spectra, spectrum_or_cal, cal_or_calref, calref_or_lamref, spectrum2, $
  spectrum3, spectrum4, spectrum5, NITER=niter, N_COR_SEC=n_cor_sec, NLAG=nlag, $
  LAMBDA=lambda, STOP=stop, FINAL_FIT=final_fit, MASK=mask, NOWARNING=nowarning, $
  FIT_OFFSETS=fit_offsets
  ;Prend en entree un spectre, sa calibration et une calibration de reference (avec son
  ; vecteur de longueurs d'onde), puis retourne :
  ; - Le spectre deplace pour concorder avec le vecteur de longueurs d'onde (par defaut).
  ;   Dans ce cas, spectrum2, spectrum3 et spectrum4 seront aussi decales
  ; - Le vecteur de longueurs d'onde deplace pour concorder avec le spectre (si /LAMBDA est
  ;     utilise).
  ;LAMBDA est un mot-cle desuet
  
  forward_function poly_fit, killnan, c_correlate
  
  if keyword_set(lambda) then begin
    calibrations = spectrum_or_cal
    calref = cal_or_calref
    lamref = calref_or_lamref
  endif else begin
    spectrum = spectrum_or_cal
    calibrations = cal_or_calref
    calref = calref_or_lamref
  endelse
  
  if total(calref) eq 0. then $
    if ~keyword_set(nowarning) then message, ' CALREF contains only zeroes !' else return, spectrum_or_cal
  
  if total(finite(calref)) eq 0. then $
    if ~keyword_set(nowarning) then message, ' CALREF contains only NaNs !' else return, spectrum_or_cal
  
  if total(calibrations) eq 0. then $
    if ~keyword_set(nowarning) then message, ' CALIBRATIONS contains only zeroes !' else return, spectrum_or_cal
    
  if total(finite(calibrations)) eq 0. then $
    if ~keyword_set(nowarning) then message, ' CALIBRATIONS contains only NaNs !' else return, spectrum_or_cal
  
  ;Nombre d'iterations pour chercher la meilleure correlation
  if ~keyword_set(niter) then niter = 3L
  
  ;Nombre de sections en lesquelles on brisera les spectres pour que la correction
  ; soit d'un ordre superieur a 0 (c'est a dire un etirement possible en plus d'un deplacement)
  ; Il faut s'assurer d'avoir des features dans chaque section !
  if ~keyword_set(n_cor_sec) then $
    n_cor_sec = 2L
  
  ;Nombre de decalages testes
  if ~keyword_set(nlag) then $
    nlag = 801L
  
  if keyword_set(mask) and n_elements(mask) ne nlag then $
    message, ' MASK must have NLAG elements !'
  
  ;Nombre de points a gauche et a droite du maximum de correlation
  ; pour faire le fit
  np = 5L
  
  ;Degre du fit pour trouver precisement le maximum de la correlation
  fit_deg = 2L < n_cor_sec
  
  nx = (size(calibrations))[1]
  nel_per_sec = ceil(float(nx)/float(n_cor_sec))
  lag = (findgen(nlag) - float(nlag-1L)/2.)/2.
  
  fit_offsets = fltarr(niter, fit_deg)
  calibrations_decal = calibrations
  for i=0L, niter-1L do begin
    
    ;Vecteurs pour stocker les resultats des correlations de chaque section
    offsets = fltarr(n_cor_sec)+!values.f_nan
    offsets_x = fltarr(n_cor_sec)+!values.f_nan
    cen_pix = fltarr(n_cor_sec)+!values.f_nan
    cen_wvl = fltarr(n_cor_sec)+!values.f_nan
    
    for j=0L, n_cor_sec-1L do begin
      
      ;Indices de la section consideree
      ind = (lindgen(nx))[j*nel_per_sec:((j+1L)*nel_per_sec-1L)<(nx-1L)]
      
      ;C'est le vecteur correspondant au 2e input qui sera decale
      cc = c_correlate(abs(killnan(calref[ind])),abs(killnan(calibrations_decal[ind])),lag)
      if keyword_set(mask) then cc *= mask
      
      ;On determine le point où la corrélation est maximale
      cc[0:np-1L] = !values.f_nan
      cc[-np:*] = !values.f_nan
      maxi = max(cc,ii,/nan)
      if keyword_set(stop) then stop
      if max(finite(cc)) eq 0 then begin
        if ~keyword_set(nowarning) then message, ' Correlation has yielded only NaN values !' else $
          return, spectrum_or_cal
      endif
      if ii le np or ii ge nlag-1L-np then begin
        if ~keyword_set(nowarning) then message, ' Correlation of spectra has failed because maximal correlation was on the edge of the arrays !' else $
          return, spectrum_or_cal
      endif
      
      ;On fait un fit polynomial de 2e degré dans la fonction de corrélation autour du maximum
      fit = poly_fit(lag[ii-np:ii+np],cc[ii-np:ii+np],2.)
      
      ;On détermine la position exacte du maximum (dérivée nulle)
      ;if fit_deg ne 2 then message, ' Ici la formule doit etre changee pour trouver le point maximal !'
      offset = -.5*fit[1]/fit[2]
      
      ;On stocke le resultat
      cen_pix[j] = median(ind)
      offsets[j] = offset
    endfor
    ;293.592. Min g = 216
    ;77.5920 ; plot,x & oplot,shift(ar_ref,78),color=255
    ;On fait un fit lineaire pour l'offset en fonction de la region
;Je ne comprends plus a quoi servaient les lignes suivantes
;    if min(offsets eq offsets[0]) eq 1 then begin
;      fit = fltarr(fit_deg)
;      fit[0] += offsets[0]
;      yft = offsets
;    endif else $
      fit = poly_fit(cen_pix,offsets,(fit_deg-1L),yfit=yfit)
    fit_offsets[i,*] = fit
    
    ;On effectue cet offset
    if max(abs(yfit)) ge 0.03 then begin
      g = where(finite(calibrations_decal), ng)
      if ng ne 0 then begin
        calibrations_off = calibrations_decal * !values.f_nan
        ;temp = spline((lindgen(nx)-(fit[0]+lindgen(nx)*fit[1]))[g],calibrations_decal[g],(lindgen(nx))[g])
        xx = lindgen(nx)
        if fit_deg eq 1 then begin
          vstart = (xx-fit[0])[g]
        endif else begin
          vstart = (xx - total((replicate(1,nx)#fit)*(xx#replicate(1,fit_deg))^(replicate(1,nx)#findgen(fit_deg)),2))[g]
          ;vstart = (lindgen(nx)-(final_fit[0]+lindgen(nx)*final_fit[1]))[g]
        endelse
        ;vstart = (lindgen(nx)-(fit[0]+lindgen(nx)*fit[1]))[g]
        vout = (lindgen(nx))[g]
        temp = interpol(calibrations_decal[g],vstart,vout,/nan)
        bad = where(vout gt max(vstart,/nan) or vout lt min(vstart,/nan), nbad)
        if nbad ne 0 then temp[bad] = !values.f_nan
        calibrations_off[g] = temp
        ;test_pad = spline((lindgen(nx)-(fit[0]+lindgen(nx)*fit[1]))[g],(findgen(nx))[g],(lindgen(nx))[g])
;        bad = where(test_pad lt 0. or test_pad gt (nx-1.), nbad)
;        if keyword_set(stop) then stop
;        if nbad ne 0 then calibrations_off[g[bad]] = !values.f_nan
        calibrations_decal = calibrations_off
      endif
    endif
  endfor
  final_fit = total(fit_offsets,1,/nan)
  
  if keyword_set(lambda) then begin
    lambda_return = lamref * !values.f_nan
    g = where(finite(lamref),ng)
    if ng eq 0 then message, ' All the reference wavelength are NaNs !'
    xx = lindgen(nx)
    if fit_deg eq 1 then begin
      vstart = (xx-final_fit[0])[g]
    endif else begin
      vstart = (xx - total((replicate(1,nx)#final_fit)*(xx#replicate(1,fit_deg))^(replicate(1,nx)#findgen(fit_deg)),2))[g]
      ;vstart = (lindgen(nx)-(final_fit[0]+lindgen(nx)*final_fit[1]))[g]
    endelse
    ;vstart = (lindgen(nx)+(final_fit[0]+lindgen(nx)*final_fit[1]))[g]
    vout = (lindgen(nx))[g]
    temp = interpol(lamref[g],vstart,vout,/nan)
    ;temp = spline((lindgen(nx)+(final_fit[0]+lindgen(nx)*final_fit[1]))[g],lamref[g],(lindgen(nx))[g])
    bad = where(vout gt max(vstart,/nan) or vout lt min(vstart,/nan), nbad)
    if nbad ne 0 then temp[bad] = !values.f_nan
    lambda_return[g] = temp
    return, lambda_return
  endif
  
  ;On effectue les offsets finaux
  spectrum_off = spectrum * !values.f_nan
  g = where(finite(calibrations), ng)
  if ng ne 0 then begin
    xx = lindgen(nx)
    if fit_deg eq 1 then begin
      vstart = (xx-final_fit[0])[g]
    endif else begin
      vstart = (xx - total((replicate(1,nx)#final_fit)*(xx#replicate(1,fit_deg))^(replicate(1,nx)#findgen(fit_deg)),2))[g]
      ;vstart = (lindgen(nx)-(final_fit[0]+lindgen(nx)*final_fit[1]))[g]
    endelse
    vout = (lindgen(nx))[g]
    temp = interpol(spectrum[g],vstart,vout,/nan)
    bad = where(vout gt max(vstart,/nan) or vout lt min(vstart,/nan), nbad)
    if nbad ne 0 then temp[bad] = !values.f_nan
    spectrum_off[g] = temp
    if keyword_set(spectrum2) then begin
      temp = interpol(spectrum2[g],vstart,vout,/nan)
      if nbad ne 0 then temp[bad] = !values.f_nan
      spectrum2[g] = temp
    endif
    if keyword_set(spectrum3) then begin
      temp = interpol(spectrum3[g],vstart,vout,/nan)
      if nbad ne 0 then temp[bad] = !values.f_nan
      spectrum3[g] = temp
    endif
    if keyword_set(spectrum4) then begin
      temp = interpol(spectrum4[g],vstart,vout,/nan)
      if nbad ne 0 then temp[bad] = !values.f_nan
      spectrum4[g] = temp
    endif
    if keyword_set(spectrum5) then begin
      temp = interpol(spectrum5[g],vstart,vout,/nan)
      if nbad ne 0 then temp[bad] = !values.f_nan
      spectrum5[g] = temp
    endif
  endif
  return, spectrum_off
End