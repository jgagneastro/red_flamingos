;diff contains the science image. The spectral trace must be horizontal !
;SMOOTH is used to smooth the image before detecting the trace (does *not* smooth the data)
Function extract_trace, diff_in, SMOOTH=smooth, DISPLAY=display, $
  IMSKY=sky, POST=post, SUBIND=subind, OPTIMAL=optimal, $
  NOSKY=nosky, FIT_PARAMS=fit_params, PARAMS_SKY=params_sky, $
  FIXED_SEP=fixed_sep, TITLE=title, SK_IMAGE=sk_image, CR_ITER=cr_iter, $
  MEDIAN=median, SKY_DEG=sky_deg, LOOPFIT=loopfit, SKY_2D=sky_2d, $
  BINARY=binary, UNRESOLVED=unresolved, READ_NOISE=read_noise, FIBERSCRAMBLER=fiberscrambler, $
  SPEC_2D=spec_2d, EMPIRIC_PROFILE=empiric_profile, EMP_TRACE=emp_trace, SKY_OUTPUT=sky_output, $
  ROBUST_SHIFT=robust_shift, GAIN=gain, SIGMA_THRESHOLD_OPTIMAL=sigma_threshold_optimal, $
  MASK_OPTIMAL_2D=mask_optimal_2d, EMP_SKY=emp_sky, HORNE_PROFILE=horne_profile, $
  CORRELATIONPROFILE=correlationprofile, ENVELOPES_2D=envelopes_2d, SKY_FIT=sky_fit, KEEPBOTH=keepboth, CHI2S=chi2s, NPTS=npts
  ;POST : On utilise SUBIND et FIT_PARAMS au lieu de detecter les traces
  
  forward_function detect_trace, gamma, nmoffat_sky_scrambler, nmoffat_gmos, optimal_extract_PMassey, mpfitfun, nmoffat_sky, extract_trace 
  
  if keyword_set(unresolved) and ~keyword_set(binary) then binary = 1
  diff = diff_in
  
  ;If read_noise is not defined, then set it as zero
  if ~keyword_set(read_noise) then $
    read_noise = 0.
  ;Set the default number of iterations for the optimal extraction algorithm
  if ~keyword_set(cr_iter) then $
    cr_iter = 6L
  if ~keyword_set(sky_fit) then sky_fit = !NULL
  
  ;Number of spectral pixels
  nx = (size(diff))[1]
  ;Number of spatial pixels
  ny = (size(diff))[2]
  xind = findgen(ny)
  
  if keyword_set(post) then begin
    pos_trace = fit_params[subind]
    peak_value = fit_params[subind-1L]
    peak_width = fit_params[subind+1L]
    peak_indice = fit_params[subind+2L]
  endif else begin
    pos_trace = detect_trace(diff, SMOOTH=smooth, DISPLAY=display, $
      PEAK_VALUE=peak_value, SKY_IND=sky_ind, $
      PEAK_WIDTH=peak_width, PEAK_INDICE=peak_indice, $
      FIT_PARAMS=fit_params, SUBIND=subind, SKY_DEG=sky_deg, $
      SKY_FIT=sky_fit, BINARY=binary, UNRESOLVED=0, $
      FIBERSCRAMBLER=fiberscrambler, EMPIRIC_PROFILE=empiric_profile, $
      EMP_TRACE=emp_trace, EMP_SKY=emp_sky)
    ind = 0.
  endelse
  
  spectrum = fltarr(nx,1L+keyword_set(binary)) + !values.f_nan
  spectrumold = fltarr(nx,1L+keyword_set(binary)) + !values.f_nan
  
  ;Built Moffat envelopes with integral normalized to 1
  norm = ( gamma(peak_indice-.5)/gamma(peak_indice)*sqrt(!pi)*peak_width )
  ;Si la fonction Gamma tend vers l'infini, on doit utiliser la limite.
  bad = where(~finite(norm),nbad)
  if nbad ne 0 then begin
    bbad = where(gamma(peak_indice) eq !values.f_infinity, nbbad)
    if nbbad ne 0 then $
      norm[bad[bbad]] = 1./(peak_indice[bad[bbad]])^.5*sqrt(!pi)*peak_width[bad[bbad]]
    if nbbad ne nbad then message, ' Some Moffat functions for the trace could not be properly normalized !' 
  endif
  bad = where(peak_indice lt .5, nbad)
  if nbad ne 0 then begin
    norm[bad] = 1.
    message, ' Some Moffat functions for the trace could not be properly normalized !', /continue
  endif
  
  if keyword_set(emp_trace) then begin
    envelopes = emp_trace
  endif else begin
    envelopes = fltarr(ny,1L+long(keyword_set(binary)))
    ;envelopes[*,0] = nmoffat_gmos(xind,[0.,1.,fit_params[2:4]]) / norm[0L]
    if keyword_set(fiberscrambler) then begin
      envelopes[*,0] = nmoffat_sky_scrambler(xind,[fit_params[0:1],0.,fit_params[3:-2],0.])
      envelopes[*,0] /= total(envelopes[*,0],/nan,/double)
    endif else $
      envelopes[*,0] = nmoffat_gmos(xind,[0.,1.,fit_params[2:4]]) / norm[0L]
    if keyword_set(binary) then $
      envelopes[*,1] = nmoffat_gmos(xind,[0.,1.,fit_params[6:8]]) / norm[1L]
    
    ;envelope = sign(peak_value[0L])*1./(((xind-pos_trace[0L])/peak_width[0L])^2 + 1)^peak_indice[0L] / norm[0L]
    if min(peak_width,/nan) lt 0. then $
      message, ' Some extraction envelopes have negative widths !'
    if total(envelopes[*,0]) eq 0. then $
      message, ' The extraction envelope contains only zeroes !'
    if keyword_set(binary) then $
      if total(envelopes[*,1]) eq 0. then $
        message, ' The secondary extraction envelope contains only zeroes !'
      
    ;On coupe les regions ou l'enveloppe est a moins de 1%
    ;Cut regions where envelopes fall below 1%
    bad = where(envelopes[*,0] lt 1./norm[0]*.01, nbad)
    if nbad ne 0L then envelopes[bad,0] = 0.
    if keyword_set(binary) then begin
      bad = where(envelopes[*,1] lt 1./norm[1]*.01, nbad)
      if nbad ne 0L then envelopes[bad,1] = 0.
    endif
  endelse
  
  ;Remove sky polynomial in the data
  if keyword_set(sky_deg) then begin
    y = findgen(ny)
    sky_1D = fit_params[0] + total(sky_fit#(dblarr(ny)+1D) * ((dblarr(sky_deg)+1D)#y)^((dindgen(sky_deg)+1D)#(dblarr(ny)+1D)), 1, /nan)
    sky_2D = (fltarr(nx)+1.)#sky_1D
  endif else sky_2D = fit_params[0] + float(finite(diff))*0.
  
  ;Remove sky from diff image
  diffsky = diff - sky_2D
  
  ;J. Gagne 2015 Feb. 30. Ensure envelopes are positive
  envelopes >= 0.
  
  ;Normal extraction
  if ~keyword_set(unresolved) then begin
    
    ;Build 2D extraction envelopes
    envelopes_2D = fltarr(nx,ny,1L+long(keyword_set(binary)))+!values.f_nan
    envelopes_2D[*,*,0] = make_array(nx,value=1.,/float)#envelopes[*,0]
    if keyword_set(binary) then $
      envelopes_2D[*,*,1] = make_array(nx,value=1.,/float)#envelopes[*,1]
    
    if keyword_set(optimal) then begin
      
      thresh_profile_count = .08
      
      if keyword_set(horne_profile) then begin
        
        ;Start with a standard extraction
        ;spectrum_estim = sign(peak_value[0L])*total(envelopes_2D*diff,2,/nan)/total(envelopes_2D,2,/nan)
        ;Blaze correction
        ;sm_spectrum_estim = correct_cshell_blaze(supersmooth(spectrum_estim,3), NDEG=1, NSMOOTH=20L)
        
        ;Create a smooth 2D profile
        profile2d = dblarr(nx,ny)
        residual2d = dblarr(nx,ny)
        i = !NULL
        ;if keyword_set(mask_optimal_2d) then $
        ;  maskprofile2d = mask_optimal_2d else $
        maskprofile2d = dblarr(nx,ny)+1
        var2d = dblarr(nx,ny)
;        bad = where(gg[0:-2]+1 ne (shift(gg,-1))[0:-2], nbad)+1
;        if nbad ne 0L then remove, bad, gg
;        if gg[-1L] ne (gg[-2L]+1) then remove, n_elements(gg)-1, gg
;        if gg[0L] ne (gg[1L]-1) then remove, 0L, gg
;        ngg = n_elements(gg)
        
        ;Eliminating the obvious bad pixels
        diff2 = diff
        bad = where(~finite(diff2), nbad)
        if nbad ne 0L then begin & $
          diff2m = diff2*!values.f_nan & $
          for ji=0L, ny-1L do $
            diff2m[*,ji] = median(diff2[*,ji],12L) & $
          diff2[bad] = diff2m[bad] & $
        endif
        bad = where(~finite(diff2), nbad)
        if nbad ne 0L then diff2[bad] = median(diff2)
        bad = where(abs(diff2)/median(diff2) gt 50. or diff2/abs(peak_value[0L]) gt 5., nbad)
        if nbad ne 0L then diff2[bad] = !values.f_nan
        bad = where(~finite(diff2), nbad)
        if nbad ne 0L then begin & $
          diff2m = diff2*!values.f_nan & $
          for ji=0L, ny-1L do $
            diff2m[*,ji] = median(diff2[*,ji],12L) & $
          diff2[bad] = diff2m[bad] & $
        endif
        bad = where(~finite(diff2), nbad)
        if nbad ne 0L then diff2[bad] = median(diff2)
        diff2m = diff2*!values.f_nan
        for ji=0L, ny-1L do $
          diff2m[*,ji] = median(diff2[*,ji],12L)
        bad = where((diff2-diff2m)/diff2m gt 10., nbad)
        if nbad ne 0L then diff2[bad] = !values.f_nan
        bad = where(~finite(diff2), nbad)
        if nbad ne 0L then begin & $
          diff2m = diff2*!values.f_nan & $
          for ji=0L, ny-1L do $
            diff2m[*,ji] = median(diff2[*,ji],12L) & $
          diff2[bad] = diff2m[bad] & $
        endif
        bad = where(~finite(diff2), nbad)
        if nbad ne 0L then diff2[bad] = median(diff2)
        
        ;Remove the sky
        diff2 -= (dblarr(nx)+1d0)#emp_sky
        
;        ;Building the extraction profiles
;        for j=0L, ngg-1L do begin & $
;          data_j = reform(diff2[*,gg[j]])*sign(peak_value[0L]) & $
;          var_j = abs(data_j) / gain + read_noise^2 & $
;          maskv = !NULL & $
;          ;maskv = maskprofile2d[*,gg[j]] & $
;          crv = !NULL & $
;          errflag = 0 & $
;          profile2d[*,gg[j]] = myprocvect(data_j, varv=var_j, crv=crv, Q=gain, v0=read_noise^2, $
;            skyvarv=replicate(read_noise^2,nx),maskv=maskv, multv=replicate(median(spectrum_estim),nx)) & $;, multv=sm_spectrum_estim) & $
;          if keyword_set(maskv) then $
;            maskprofile2d[*,gg[j]] = maskv & $
;        endfor
        
        ;Remove NaNs
;        bad = where(~finite(profile2d), nbad)
;        if nbad ne 0L then profile2d[bad] = 0
        
        ;mask2d = fltarr(nx,ny)+!values.f_nan
        
        ;Recreate an empirical profile from "diff2"
        diff3 = diff2
        diff3 >= 0.
        diff3 /= (total(diff3,2)#(dblarr(ny)+1d0))
        nsm = 35
        xx = lindgen(nx)
        profile2d *= !values.d_nan
        for j=0L, ny-1L do begin & $
          hyper_sm = supersmooth(diff3[*,j],nsm) & $
          hyper_sm[0:nsm-1L] = !values.d_nan & $
          hyper_sm[-nsm-2L:*] = !values.d_nan & $
          gf = where(finite(hyper_sm), ngf) & $
          ;First approx
          coeff = ladfit(xx[gf],hyper_sm[gf]) & $
          profile2d[*,j] = poly(xx,coeff) & $
        endfor
        
        for j=0L, nx-1L do begin
        
          gg = where(profile2d[j,*] gt thresh_profile_count*max(profile2d[j,*],/nan), ngg); and median(envelopes_2d,dim=1) ne 0
    
          ;Remove satellites
          ;Where is the center expected to be ?
          void = max(total(envelopes_2d,1),wmax)
          ;Identify which gg this is
          gw = (where(gg eq wmax, ngw))[0L]
          if ngw eq 0 then message, ' The maximum of the diff2 image seems to happen far from where it''s supposed to !'
          ;Find out where the monotonic sequence is broken
          ibroken = where((lindgen(ngg)-gw+wmax) ne gg, nbroken)
          if nbroken ne 0L then begin
            ;Find out which broken index is the first one AFTER the expected peak
            iright = (where(gg[ibroken] gt wmax, nright))[0L]
            ;Find out which broken index is the first one BEFORE the expected peak
            ileft = (where(gg[ibroken] lt wmax, nleft))[-1L]
            ;Restrain to central gg's
            if nright ne 0L then $
              gg = gg[0L:(ibroken[iright]-1L)]
            if nleft ne 0L then $
              gg = gg[(ibroken[ileft]+1L):*]
          endif
          ngg = n_elements(gg)
          
          gg2 = lindgen(nx)
          remove, gg, gg2
          
          profile2d[j,gg2] = 0.
        endfor
        
        profile2d >= 0.
        profile2d /= (total(profile2d,2)#(dblarr(ny)+1d0))
        
        if keyword_set(correlationprofile) then begin
          addinterp = 1d0 & $
          lag = dindgen(13)-6
          np = 2
          lags = dblarr(nx)+!values.d_nan & $
          for i=0L, nx-1L do begin & $
            ;Compute cross-correlation
            cc = c_correlate(reform(envelopes_2d[i,*]),reform(profile2d[i,*]),lag) & $
            ;Find approximate maximum
            void = max(cc,/nan,wmax) & $
              
            ;Refine it
            fit = reform(poly_fit(lag[wmax-np:wmax+np],cc[wmax-np:wmax+np],2.,yfit=g))
            ;  ;fit2 = mpfitfun('poly',offset[ii-np:ii+np],cc[ii-np:ii+np],0.,[fit,0.],WEIGHTS=1D, STATUS=status, /NAN)
            ;
            ;  ;On détermine la position exacte du maximum (dérivée nulle)
            ;  ;if fit_deg ne 2 then message, ' Ici la formule doit etre changee pour trouver le point maximal !'
            ;  bestoffset = -.5*fit[1]/fit[2]
             
             lags[i] = -.5*fit[1]/fit[2] & $ 
            ;lags[i] = lag[wmax] & $
          endfor
          
          if max(abs(lags),/nan) ne 0. then begin
            ;Do a robust polynomial fit
            xx = dindgen(nx)
            coeff = ladfit(xx,lags)
            if total(abs(coeff)) ne 0. then begin
              ;Re-interpolate extraction profile
              yy = dindgen(ny)
              for i=0L, nx-1L do $
                envelopes_2d[i,*] = interpol2(reform(envelopes_2d[i,*]),yy,yy-poly(i,coeff),badval=0d0)
            endif; else keep the envelopes_2d
          endif; else keep the envelopes_2d
          
        endif else begin
          envelopes_2d = profile2d
        endelse
        
      endif else $
        profile2d = envelopes_2d[*,*,0]
      
      
      if optimal eq 1 then begin
        optimal_extract_PMassey, diff, sky_2d, envelopes_2D[*,*,0], read_noise, spext, ROBUST_SHIFT=robust_shift
        spectrum[*,0] = spext
        if keyword_set(binary) then begin
          optimal_extract_PMassey, diff, sky_2d, envelopes_2D[*,*,1], read_noise, spext2, ROBUST_SHIFT=robust_shift
          spectrum[*,1] = spext2
        endif
      endif
      if optimal eq 2 then begin
        
        error = fltarr(nx)+!values.f_nan
        
        if keyword_set(mask_optimal_2d) then $
          mask2d = mask_optimal_2d else $
          mask2d = dblarr(nx,ny)+1
        
        nsigmas2d = fltarr(nx,ny)+!values.f_nan
        normf = fltarr(nx)+!values.f_nan
        for j=0L, nx-1L do begin & $
          env_j = sign(peak_value[0L])*reform(profile2d[j,*]) & $
          ;env_j = sign(peak_value[0L])*reform(envelopes_2d[j,*]) & $
          env_j >= 0. & $
          gg = where(env_j gt thresh_profile_count*max(env_j,/nan), ngg) & $; and median(envelopes_2d,dim=1) ne 0
          if ngg eq 0 then continue & $
          
          ;Remove satellites
          ;Where is the center expected to be ?
          void = max(total(envelopes_2d,1),wmax)
          ;Identify which gg this is
          gw = (where(gg eq wmax, ngw))[0L]
          if ngw eq 0 then message, ' The maximum of the diff2 image seems to happen far from where it''s supposed to !'
          ;Find out where the monotonic sequence is broken
          ibroken = where((lindgen(ngg)-gw+wmax) ne gg, nbroken)
          if nbroken ne 0L then begin
            ;Find out which broken index is the first one AFTER the expected peak
            iright = (where(gg[ibroken] gt wmax, nright))[0L]
            ;Find out which broken index is the first one BEFORE the expected peak
            ileft = (where(gg[ibroken] lt wmax, nleft))[-1L]
            ;Restrain to central gg's
            if nright ne 0L then $
              gg = gg[0L:(ibroken[iright]-1L)]
            if nleft ne 0L then $
              gg = gg[(ibroken[ileft]+1L):*]
          endif
          ngg = n_elements(gg)
;          bad = where(gg[0:-2]+1 ne (shift(gg,-1))[0:-2], nbad)+1 & $
;          if nbad eq n_elements(gg) then continue & $
;          if nbad ne 0L then remove, bad, gg & $
;          if n_elements(gg) eq 1L then continue & $
;          if gg[-1L] ne (gg[-2L]+1) then remove, n_elements(gg)-1, gg & $
;          if n_elements(gg) eq 1L then continue & $
;          if gg[0L] ne (gg[1L]-1) then remove, 0L, gg & $
;          if n_elements(gg) eq 1L then continue & $
;          ngg = n_elements(gg) & $
          
          env_j = reform(profile2d[j,gg]) & $
          env_j >= 0. & $
          env_j /= total(env_j,/nan) & $
          data_j = reform(diff[j,gg])*sign(peak_value[0L]) & $
          ;maskv = !NULL
          maskv = reform(mask2d[j,gg]) & $
          ;if min(maskv) eq 0 then stop
          nsigmas = !NULL & $
          
          ;Remove NaNs in data
          bad = where(~finite(data_j), nbad) & $
          if nbad eq ngg then continue & $
          if nbad ne 0L then begin & $
            data_j[bad] = 0. & $
            maskv[bad] = 0. & $
          endif & $
          spectrum[j] = optextract_1d(data_j, PROFILE=env_j, GAIN=gain, $
            THRESH=ngg/2, READNOISE=read_noise, ERROR=error_j, MASKV=maskv, $
            NSIGMAS=nsigmas, VTHRESH=sigma_threshold_optimal^2) & $
;          spectrumold[j] = optextract_1d(data_j, PROFILE=env_j, GAIN=gain, $
;            THRESH=ngg/2, READNOISE=read_noise, ERROR=error_j, MASKV=maskv, $
;            NSIGMAS=nsigmas, VTHRESH=sigma_threshold_optimal^2, /OLD) & $
          ;Enforce normalization similar to OPTIMAL=1
          normf[j] = total(env_j^2,/NAN) & $
          if keyword_set(maskv) then $
            mask2d[j,gg] = maskv & $
          if keyword_set(nsigmas) then $
            nsigmas2d[j,gg] = nsigmas & $
          if keyword_set(error_j) then $
            error[j] = error_j & $
        endfor
        
;        wset,2
;        plot,spectrum_estim/median(spectrum_estim)
;        oplot, spectrumold/median(spectrumold), color=255
;        oplot, spectrum/median(spectrum), col=rvb_hex(100,255,100)
;        stop
;        saveimage, 'SPCOMP_'+curcompdate(/path)+'_'+strtrim(long(getseed()*1d4),2)+'.png', /PNG
        
        ;Enforce normalization similar to OPTIMAL=1
        spectrum *= mean(normf,/nan)
        
        extraction_profiles = profile2d * sign(peak_value[0L])
        
        ;Replace NaNs (NOPE)
        bad = where(~finite(spectrum), nbad)
        if nbad ne 0L then stop
      endif
    endif else begin
      ;Normal extraction (non-optimal)
      spectrum[*,0] = total(envelopes_2D[*,*,0]^2*diffsky,2,/nan)/total(envelopes_2D[*,*,0]^2,2,/nan)
      if keyword_set(binary) then spectrum[*,1] = total(envelopes_2D[*,*,1]^2*diffsky,2,/nan)/total(envelopes_2D[*,*,1]^2,2,/nan)
    endelse
    
    ;x_extobjopt,diff,lindgen(nx)#(make_array(ny,value=1L,/long)),finite(diff)*0.+1./stddev(diffsky,/nan)^2,sky_2D,[109,69],fin
    ;x_extobjopt,transpose(diff),transpose(lindgen(nx)#(make_array(ny,value=1L,/long))),transpose(finite(diff)*0.+stddev(diffsky,/nan)),transpose(sky_2D),[109,69],fin
    
  endif else begin
    ;Extraction for barely resolved binaries
    parinfo = replicate({VALUE:!values.f_nan,FIXED:1B,LIMITED:[0B,0B],LIMITS:[0.,0.]},n_elements(fit_params));-n_elements(sky_fit))
    ;Amplitudes can vary
    parinfo[subind-1L].FIXED = 0B
    parinfo[subind-1L].LIMITED[0] = 1
    parinfo[subind-1L].LIMITS[0] = 0.
    ;Width can vary by 15%
    var = .15
    parinfo[subind+1L].FIXED = 0B
    parinfo[subind+1L].LIMITED = bytarr(2,n_elements(subind))+1B
    parinfo[subind+1L].LIMITS[0] = (fit_params[subind+1L]*(1.-var))>0.
    parinfo[subind+1L].LIMITS[1] = fit_params[subind+1L]*(1.+var)
    ;Moffat index can vary by 5%
    var = .05
    parinfo[subind+2L].FIXED = 0B
    parinfo[subind+2L].LIMITED = bytarr(2,n_elements(subind))+1B
    parinfo[subind+2L].LIMITS[0] = (fit_params[subind+2L]*(1.-var))>0.51
    parinfo[subind+2L].LIMITS[1] = fit_params[subind+2L]*(1.+var)
    ;Zero level can vary
    parinfo[0L].FIXED = 0B
    tmp = fit_params[0:n_elements(fit_params)-1L]
    ;if keyword_set(sky_fit) then tmp = [tmp, 0.]
    parinfo.VALUE = double(tmp)
    ;Sky slope can vary
    parinfo[-1L].FIXED = 0B
;    parinfo[-1].LIMITED = 1B
    ;Sky slope and zero levels can't vary
    if keyword_set(sky_fit) then begin
      meddiff = median(diffsky,dim=1)
      mslope = abs(max(meddiff,/nan)-min(meddiff,/nan))/abs(max(xind,/nan)-min(xind,/nan))
      parinfo[-1].LIMITS[0] = -mslope/4.
      parinfo[-1].LIMITS[1] = mslope/4.
      ;parinfo[0L].FIXED = 1B
      ;parinfo[-1L].FIXED = 1B
      parinfo[0].VALUE = double(fit_params[0])
      parinfo[-1].VALUE = double(fit_params[-1])
    endif
    
    ;Make the double Moffat fit on each spatial column 
    print, ' Using special reduction for unresolved binaries... '
    nys = make_array(ny,value=1.,/float)
    chi2s = dblarr(nx)+!values.d_nan
    npts = lonarr(nx)
    yfits = dblarr(ny,nx)+!values.d_nan
    allpars = dblarr(nx,n_elements(parinfo))+!values.d_nan
    for ks=0L, nx-1L do begin
      if total(finite(diffsky[ks,*])) eq 0 then continue
      ;Avoid changing the parinfo structure
      parinfoi = parinfo
      ;Make the fit
      pari = fit_params + !values.f_nan
      gg = where(finite(diffsky[ks,*]), ngg)
      npts[ks] = ngg
      pari = mpfitfun('nmoffat_sky', xind[gg], diffsky[ks,gg], 0., WEIGHTS=1D, PARINFO=parinfoi, STATUS=status, /NAN, yfit=yfit, /QUIET)
      ;If the fit failed, don't store results
      if status le 0L then continue
      chi2s[ks] = total((diffsky[ks,gg]-yfit)^2,/nan)
      yfits[gg,ks] = yfit
      ;MODELI is the analog of DIFF when using regular extraction
      allpars[ks,*] = pari
      modeli0 = nmoffat_sky(xind,[0.,pari[1:4]])
      modeli1 = nmoffat_sky(xind,[0.,pari[5:8]])
      spectrum[ks,0L] = total(modeli0)
      spectrum[ks,1L] = total(modeli1)
      ;spectrum[ks,0L] = total(modeli0*envelopes[*,0],/nan)/total(envelopes[*,0],/nan)
      ;spectrum[ks,1L] = total(modeli1*envelopes[*,1],/nan)/total(envelopes[*,1],/nan)
    endfor
    ;Identify bad fits
    chi2s0 = chi2s
    spectrum0 = spectrum
    bad = where(chi2s gt robust_sigma(chi2s,/nan)*30d0, nbad)
    for kbad=0L, nbad-1L do begin
      ks = bad[kbad]
      if total(finite(diffsky[ks,*])) eq 0 then continue
      ;Avoid changing the parinfo structure
      parinfoi = parinfo
      ;Make the fit
      pari = fit_params + !values.f_nan
      gg = where(finite(diffsky[ks,*]), ngg)
      ;pari = mpfitfun('nmoffat_sky', xind[gg], diffsky[ks,gg], 0., WEIGHTS=1D, PARINFO=parinfoi, STATUS=status, /NAN, yfit=yfit, /QUIET)
      ;If the fit failed, don't store results
      if status le 0L then continue
      badk = where(abs(diffsky[ks,gg]-yfits[gg,ks]) gt 5d0*robust_sigma(diffsky[ks,gg]-yfits[gg,ks],/nan), nbadk)
      if nbadk eq 0 then continue
      if nbadk eq ngg then continue
      remove, badk, gg
      ngg = n_elements(gg)
      npts[ks] = ngg
      pari = mpfitfun('nmoffat_sky', xind[gg], diffsky[ks,gg], 0., WEIGHTS=1D, PARINFO=parinfoi, STATUS=status, /NAN, yfit=yfit, /QUIET)
      chi2s[ks] = total((diffsky[ks,gg]-yfit)^2,/nan)
      ;MODELI is the analog of DIFF when using regular extraction
      allpars[ks,*] = pari
      modeli0 = nmoffat_sky(xind,[0.,pari[1:4]])
      modeli1 = nmoffat_sky(xind,[0.,pari[5:8]])
      spectrum[ks,0L] = total(modeli0)
      spectrum[ks,1L] = total(modeli1)
    endfor
    ;if nbad ge 3 then stop
  endelse
  
  ;Extract sky around the trace
  if keyword_set(sky_ind) and ~keyword_set(nosky) then begin
    g = where(sky_ind-pos_trace[0L] ge 0, ng)
    if ng eq 0 then g = where(sky_ind-pos_trace[0L] lt 0)
    pos_sky = pos_trace + median((sky_ind-pos_trace[0L])[g])
    params_sky = fit_params
    params_sky[subind] = pos_sky
    sky = extract_trace(diff, /POST, /NOSKY, DISPLAY=display, $
      FIT_PARAMS=params_sky, SUBIND=subind, READ_NOISE=read_noise)
    sky_output = sky
  endif
  
  if keyword_set(keepboth) then begin
    fit_paramsout = fit_params
    return, spectrum
  endif
  
  ;Choose which binary component to keep
  if keyword_set(binary) then begin
    if binary eq 1L then begin
      spectrum = spectrum[*,0]
      fit_paramsout = fit_params[0:4L]
      if n_elements(fit_params) ge 10L then fit_paramsout=[fit_paramsout,fit_params[9:*]]
      subind = subind[0L]
    endif
    if binary eq 2L then begin
      spectrum = spectrum[*,1]
      fit_paramsout = [fit_params[0L],fit_params[5:8L]]
      if n_elements(fit_params) ge 10L then fit_paramsout=[fit_paramsout,fit_params[9:*]]
      subind = subind[1L]
    endif
    fit_params = fit_paramsout
  endif
  
  return, spectrum
  
End