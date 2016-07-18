Function extract_traces_AB, diff, SMOOTH=smooth, DISPLAY=display, $
  BINARY=binary, SKY=sky, UNRESOLVED=unresolved, POST=post, SUBIND=subind, $
  FIT_PARAMS=fit_params, STOP=stop, ONLY_AMPLITUDE=only_amplitude, XRANGE=xrange, $
  PARAMS_SKY=params_sky, TELLURIC=telluric, OPTIMAL=optimal, READ_NOISE=read_noise, $
  CR_ITER=cr_iter, FIXED_SEP=fixed_sep, TITLE=title, DIFFS=diffs, $
  ENVELOPES=envelopes, GAIN=gain, ERROR=error, $
  SIGMA_THRESHOLD_OPTIMAL=sigma_threshold_optimal, EXTRACTION_PROFILES=extraction_profiles
  ;POST : On utilise SUBIND et FIT_PARAMS au lieu de detecter les traces
  ;TELLURIC : Cette option soustrait les telluriques en utilisant la région du ciel (pratique pour de gros CCDs par exemple, où la soustraction tellurique précédente n'est pas nécessairement efficace autour de la trace de science)
  ;OPTIMAL : Pour faire une "optimal extraction". Dans ce cas il faut donner READ_NOISE et SK_IMAGE (obtenue avec im1 < im2)
  ;Nombre d'iterations pour enlever les rayons cosmiques

  forward_function detect_traces_AB, mpfitfun, extract_trace, optextract_1d, myprocvect

  nx = (size(diff))[1]
  ny = (size(diff))[2]
  xind = findgen(ny)

  if keyword_set(optimal) and ~keyword_set(sigma_threshold_optimal) then $
    sigma_threshold_optimal = 7.

  if keyword_set(post) then begin
    pos_traces = fit_params[subind]
    peak_values = fit_params[subind-1L]
    peak_widths = fit_params[subind+1L]
    peak_indices = fit_params[subind+2L]
  endif else begin
    ;Here I replaced diff with "diff - median thing" (J. Gagne April 28, 2015)
    diffsub = diff-median(diff,dim=2)#make_array((size(diff))[2],value=1d0,/double)
    pos_traces = detect_traces_AB(diffsub, SMOOTH=smooth, DISPLAY=display, $
      BINARY=binary, PEAK_VALUES=peak_values, SKY_IND=sky_ind, $
      PEAK_WIDTHS=peak_widths, PEAK_INDICES=peak_indices, $
      FIT_PARAMS=fit_params, SUBIND=subind, STOP=stop, UNRESOLVED=unresolved,$
      ONLY_AMPLITUDE=only_amplitude, XRANGE=xrange, FIXED_SEP=fixed_sep, TITLE=title, DIFFS=diffs)
  endelse
  n_traces = n_elements(pos_traces)
  
  ;Remove tellurics precisely around the trace, using sky positions near the trace.
  if keyword_set(telluric) then begin
    
    message, ' Must determine width and position of traces from the trace itself, not the fitting parameters (use abs value of trace)'
    
    ;Build indices around each trace
    indices = ptrarr(n_traces,/ALLOCATE)
    medregs = ptrarr(n_traces,/ALLOCATE)
    madregs = ptrarr(n_traces,/ALLOCATE)
    for i=0L, n_traces-1L do begin
      width = moff_widths[i];fit_params[subind[i]+1L]
      ind0 = (long(pos_traces[i]-6.*width) > 0L)
      ind1 = (long(pos_traces[i]+6.*width) < (ny-1L) )
      ind = lindgen(ind1-ind0+1L)+ind0
      ;Remove indices that overlap one of the traces
      ;sigma/sqrt(2*v-2.266065)
      for j=0L, n_traces-1L do begin
        ;ssig = moff_widths[j];fit_params[subind[j]+1L]/sqrt(2*fit_params[subind[j]+2L]+2.266065)
        bad = where(abs(ind - pos_traces[j]) lt 3.*moff_widths[j], nbad)
        if nbad ne 0 then remove, bad, ind
      endfor
      *indices[i] = ind
      *medregs[i] = median(diff[*,ind],dim=2L)
      if keyword_set(sk_image) then $
        *madregs[i] = median(abs(diff[*,ind]-median(diff[*,ind],dim=2L)#make_array(n_elements(ind),value=1.,/float)),dim=2L)
    endfor

    ;Then, correct 6-sigma around each trace with the results, carefully avoiding to correct twice the same spots
    corrected_indices = []
    for i=0L, n_traces-1L do begin
      if *medregs[i] eq !NULL then continue
      ind0 = (long(pos_traces[i]-8.*width) > 0L)
      ind1 = (long(pos_traces[i]+8.*width) < (ny-1L) )
      ind = lindgen(ind1-ind0+1L)+ind0
      if n_elements(corrected_indices) ne 0 then begin
        match,corrected_indices,ind,ia,ib
        if n_elements(ib) eq n_elements(ind) then continue
        if ib[0] ne -1 then remove, ib, ind
      endif
      corrected_indices = [corrected_indices,ind]
      sortself, corrected_indices, /uniq
      nind = n_elements(ind)
      diff[*,ind] -= *medregs[i]#make_array(nind,value=1.,/float)
      ;Propagate this error on the sky image, if any
      if keyword_set(sk_image) then $
        sk_image[*,ind] += ( *madregs[i]#make_array(nind,value=1.,/float) )^2
    endfor
  endif
  
  if keyword_set(binary) then begin
    
    ;Méthode spéciale pour extraire les binaires non résolues
    
    ;On fait fitter 4 fonctions de Moffat d'un coup dans chaque position de la trace,
    ; mais seulement les amplitudes et le point zéro ont le droit de varier
    ;Parameter inputs
    parinfo = replicate({VALUE:!values.f_nan,FIXED:0B,LIMITED:[0B,0B],LIMITS:[0.,0.]},n_elements(fit_params))
    parinfo.VALUE = fit_params
    parinfo.FIXED = 1B
    parinfo[[0L,1L,5L,7L,11L]].FIXED = 0B
    ;On limite les valeurs possibles du shift des traces
    parinfo[1L].LIMITED = bytarr(2)+1B
    parinfo[1L].LIMITS[0L] = fit_params[1L] - fit_params[6]
    parinfo[1L].LIMITS[1L] = fit_params[1L] + fit_params[6]
    
    ;Create envelopes
    envelopes = AB_moffat_indiv_envelopes(xind,fit_params, /NORM)
    
    ;Pour chaque position spectrale on va faire le fit
    print, ' Using special reduction for unresolved binaries... '
    nys = make_array(ny,value=1.,/float)
    modeli = finite(diff)*0.
    spectra = dblarr(nx,n_traces)*!values.d_nan
    for ks=0L, nx-1L do begin
      ;On évite d'altérer la structure parinfo
      parinfoi = parinfo
      ;On fait le fit
      pari = fit_params + !values.f_nan
      g = where(finite(xind*reform(diff[ks,*])), ng)
      if ng eq 0L then continue
      pari = mpfitfun('AB_moffat_fixed_sep', xind[g], diff[ks,g], 0., WEIGHTS=1D, /QUIET, PARINFO=parinfoi, STATUS=status, /NAN, YFIT=yfit)
      ;If the fit failed, don't store the results
      if status le 0L then continue
      ;plot, xind[g], diff[ks,g]
      ;oplot, xind[g], yfit[g], color=255
      ;wait,0.2
      modeli = AB_moffat_indiv_envelopes(xind[g],pari)
      spectra[ks,*] = total(modeli*envelopes[g,*],1,/nan)
    endfor
    
  endif else begin
  
    ;On a estime numeriquement une relation entre la largeur d'une courbe normale et de celle de Moffat
    moff_widths = fit_params[subind+1L]/sqrt(2*fit_params[subind+2L]+2.266065)
  
    ;On construit des enveloppes de Moffat normalisées à 1 pour chaque trace
    envelopes = dblarr(ny,n_traces)+!values.d_nan
    norms = ( gamma(peak_indices-.5)/gamma(peak_indices)*sqrt(!pi)*peak_widths )
    ;Si la fonction Gamma tend vers l'infini, on doit utiliser la limite.
    bad = where(~finite(norms), nbad)
    if nbad ne 0 then begin
      for i=0L, nbad-1L do begin
        if gamma(peak_indices[bad[i]]) eq !values.f_infinity then begin
          norms[bad[i]] = 1./(peak_indices[bad[i]])^.5*sqrt(!pi)*peak_widths[bad[i]]
        endif else message, ' The Moffat function for the trace could not be properly normalized !'
      endfor
    endif
  
    bad = where(peak_indices lt .5, nbad)
    if nbad ne 0 then begin
      norms[bad] = 1.
      message, ' Some Moffat functions could not be properly normalized !', /continue
    endif
    for i=0L, n_traces-1L do $
      envelopes[*,i] = sign(peak_values[i])*1d0/(((xind-pos_traces[i])/peak_widths[i])^(2d0) + 1d0)^peak_indices[i] / norms[i]
  
    ;On coupe les regions ou l'enveloppe est a moins de 1%, ou a 2-sigma
    for i=0L, n_traces-1L do begin
      if sign(peak_values[i]) eq 1. then $
      bad = where(envelopes[*,i] lt 1./norms[i]*.01 or abs(xind-pos_traces[i]) gt 2.*peak_widths[i], nbad) $
    else $
      bad = where(envelopes[*,i] gt -1./norms[i]*.01 or abs(xind-pos_traces[i]) gt 2.*peak_widths[i], nbad)
      if nbad ne 0 then envelopes[bad,i] = 0.
    endfor
  
    ;On construit les enveloppes en 2D
    envelopes_2D = fltarr(nx,ny,n_traces)+!values.f_nan
    for i=0L, n_traces-1L do $
      envelopes_2D[*,*,i] = make_array(nx,value=1.,/float)#envelopes[*,i]
  
    ;On extrait les spectres
    spectra = fltarr(nx,n_traces)+!values.f_nan
    sky = fltarr(nx,n_traces)+!values.f_nan
    
    ;Optimal extraction
    if optimal eq !NULL then optimal = 0
    if optimal eq 2 then begin
      
      ;Start with a standard extraction
      spectra_estim = spectra*0.
      for i=0L, n_traces-1L do $
        spectra_estim[*,i] = sign(peak_values[i])*total(envelopes_2D[*,*,i]*diff,2,/nan)/total(envelopes_2D[*,*,i],2,/nan)
      
      ;Create a smooth 2D profile
      profile2d = dblarr(nx,ny,n_traces)
      residual2d = dblarr(nx,ny,n_traces)
      maskprofile2d = dblarr(nx,ny,n_traces)
      var2d = dblarr(nx,ny,n_traces)
      for i=0L, n_traces-1L do begin & $
        env_i = total(envelopes_2d[*,*,i],1) & $
        gg = where(abs(env_i) gt .01*max(env_i,/nan), ngg) & $
        if ngg eq 0L then continue & $
        for j=0L, ngg-1L do begin & $
          data_j = reform(diff[*,gg[j]])*sign(peak_values[i]) & $
          var_j = abs(data_j) / gain + read_noise^2 & $
          maskv = !NULL & $
          crv = !NULL & $
          errflag = 0 & $
          profile2d[*,gg[j],i] = myprocvect(data_j, varv=var_j, crv=crv, Q=gain, v0=read_noise^2, $
            skyvarv=replicate(read_noise^2,nx),maskv=maskv, multv=replicate(mean(supersmooth(spectra_estim[*,i],12),/nan),nx)) & $
          if keyword_set(maskv) then $
            maskprofile2d[*,gg[j],i] = maskv & $
        endfor & $
      endfor
      
      error = fltarr(nx,n_traces)+!values.f_nan & $
      mask2d = fltarr(nx,ny,n_traces)+!values.f_nan & $
      nsigmas2d = fltarr(nx,ny,n_traces)+!values.f_nan & $
      for i=0L, n_traces-1L do begin & $
        for j=0L, nx-1L do begin & $
          env_j = sign(peak_values[i])*reform(envelopes_2d[j,*,i]) & $
          gg = where(env_j gt 0.01*max(env_j,/nan), ngg) & $
          if ngg eq 0 then continue & $
          env_j = profile2d[j,gg,i] & $
          env_j >= 0. & $
          env_j /= total(env_j,/nan) & $
          data_j = reform(diff[j,gg])*sign(peak_values[i]) & $
          maskv = !NULL & $
          nsigmas = !NULL & $
          spectra[j,i] = optextract_1d(data_j, PROFILE=env_j, GAIN=gain, $
            THRESH=ngg/2, READNOISE=read_noise, ERROR=error_j, MASKV=maskv, $
            NSIGMAS=nsigmas, VTHRESH=sigma_threshold_optimal^2) & $
          if keyword_set(maskv) then $
            mask2d[j,gg,i] = maskv & $
          if keyword_set(nsigmas) then $
            nsigmas2d[j,gg,i] = nsigmas & $
          if keyword_set(error_j) then $
            error[j,i] = error_j & $
        endfor & $
      endfor
      
      extraction_profiles = profile2d
      for i=0L, n_traces-1L do $
        extraction_profiles[*,*,i] *= sign(peak_values[i])
      
      ;Replace NaNs
      bad = where(~finite(spectra), nbad)
      if nbad ne 0L then spectra[bad] = spectra_estim[bad]
      
;      ;Compare
;      wset,0
;      plot,spectra[*,0]/median(spectra[*,0])
;      oplot,spectra_estim[*,0]/median(spectra_estim[*,0]),color=255
;      oplot,spectra[*,0]/median(spectra[*,0])
;      oplot, error[*,0]/median(spectra[*,0]), color=rvb_hex(100,255,100)
;      wset, 1
;      plot, spectra[*,0]/error[*,0]
;      wset,2
;      plot,spectra[*,1]/median(spectra[*,1]);, xrange=[800,1300]
;      oplot,spectra_estim[*,1]/median(spectra_estim[*,1]),color=255
;      oplot,spectra[*,1]/median(spectra[*,1])
;      oplot, error[*,1]/median(spectra[*,1]), color=rvb_hex(100,255,100)
;      stop
      
    endif
    if optimal eq 1 then begin
      sky_2d = finite(diff)*0d0
      for i=0L, n_traces-1L do begin
        optimal_extract_PMassey, diff, sky_2d, envelopes_2D[*,*,i], read_noise, spext, ROBUST_SHIFT=robust_shift
        spectra[*,i] = spext
        ;spectra[*,i] = sign(peak_values[i])*total(envelopes_2D[*,*,i]*diff,2,/nan)/total(envelopes_2D[*,*,i],2,/nan)
      endfor
    endif
    if optimal eq 0 then begin
      
      if keyword_set(post) and keyword_set(extraction_profiles) then begin
        envelopes_2d = extraction_profiles
        envelopes_2d >= 0.
        for i=0L, n_traces-1L do $
          for j=0L, nx-1L do $
            envelopes_2d[j,*,i] /= total(envelopes_2d[j,*,i],/nan)
      endif else $
        extraction_profiles = envelopes_2d
      
      ;Standard Method (does not account for barely resolved binaries)
      for i=0L, n_traces-1L do $
        spectra[*,i] = sign(peak_values[i])*total(envelopes_2D[*,*,i]*diff,2,/nan)/total(envelopes_2D[*,*,i],2,/nan)
      
    endif
  endelse
  
  ;On extrait le ciel
  if keyword_set(sky_ind) then begin

    pos_sky = fltarr(n_traces)+!values.f_nan
    for i=0L, n_traces-1L do begin
      if keyword_set(binary) then width = max(fit_params[[6,8,10,12]],/nan) else $
        width = fit_params[subind[i]+1L]/sqrt(2*fit_params[subind[i]+2L]+2.266065)
      ind0 = (long(pos_traces[i]-15.*width) > 0L)
      ind1 = (long(pos_traces[i]+15.*width) < (ny-1L) )
      ind = lindgen(ind1-ind0+1L)+ind0
      ;Remove indices that overlap one of the traces
      for j=0L, n_traces-1L do begin
        bad = where(abs(ind - pos_traces[j]) lt 3.*width, nbad)
        if nbad ne 0 then remove, bad, ind
      endfor
      ;On choisit la position de ciel la plus éloignée de la trace
      void = max(abs(pos_traces[i]-ind),wmin)
      pos_sky[i] = ind[wmin]
    endfor
    
    sky = fltarr(nx,n_traces)+!values.f_nan
    if keyword_set(binary) then begin
      ;On cree un ciel "median" et on l'extrait
      fit_params_sky = fit_params
      fit_params_sky[[1L,2L,3L]] = [pos_sky[0],pos_sky[1]-pos_sky[0],pos_sky[2]-pos_sky[0]]
      sky_envelopes = AB_moffat_indiv_envelopes(dindgen(ny),fit_params_sky,/NORM)
      for ci=0L, n_traces-1L do $
        sky[*,ci] = total((make_array(nx,value=1d,/double)#sky_envelopes[*,ci])*diff,2,/nan,/double)
    endif else begin
      params_sky = fit_params
      params_sky[subind] = pos_sky
      for i=0L, n_traces-1L do begin
        sky_parami = [params_sky[0],params_sky[subind[i]-1L],params_sky[subind[i]],params_sky[subind[i]+1],params_sky[subind[i]+2]]
        sky[*,i] = extract_trace(diff, /POST, /NOSKY, DISPLAY=display, $
          FIT_PARAMS=sky_parami, SUBIND=2L, FIXED_SEP=fixed_sep, TITLE=title)
      endfor
    endelse
  endif
  
  return, spectra

End