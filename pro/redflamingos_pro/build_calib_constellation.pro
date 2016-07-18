Function build_calib_constellation, abscissa_in, sp_in, NSEC=nsec, NUSE=nuse, FLUX_THRESHOLD=flux_threshold, SIG_THRESHOLD=sig_threshold, $
  CONSTELLATIONS_INDICES=constellations_indices
  ;Builds "constellations" for an 1D arc spectrum, a bit like astrometrynet.
  ;NSEC is the number of nearby secondary lines that will be used
  ;NUSE is the maximal number of primary lines that will be used.
  ;Computer timing will scale as NUSE * (NSEC^2 - NSEC)/2 
  
  if ~keyword_set(sp_in) then $
    message, 'You must enter a calibration spectrum (1D) !'
  if ~keyword_set(abscissa_in) then $
    message, 'You must enter an abscissa array (typically pixels or wavelength) !'
  
  ;Pass parameters
  sp = sp_in
  abscissa = abscissa_in
  
  ;Normalize spectrum
  sp /= max(sp,/nan)
  
  ;Default keywords
  if ~keyword_set(nsec) then nsec = 10L
  if ~keyword_set(nuse) then nuse = 50L
  
  ;Detect all peaks
  peaks = find_1d_peaks(sp, FLUX_THRESHOLD=flux_threshold, SIG_THRESHOLD=sig_threshold)
  
  ;Limit the total number of peaks (choose the brightest first)
  if n_elements(peaks) gt nuse then $
    peaks = (peaks[reverse(sort(sp[peaks]))])[0L:nuse-1L]
  npeaks = n_elements(peaks)
  nsec <= npeaks
  
  ;Loop on peaks to build the constellations
  nconst = npeaks*(nsec^2-nsec)/2L
  constellations_indices = lonarr(nconst,3L)-1L
  l = 0L
  for i=0L, npeaks-1L do begin
    ;Choose the "NSEC" closest lines
    dist = abs(float(peaks[i] - peaks))
    dist[i] = !values.f_nan
    sec_peaks_i = (peaks[sort(dist)])[0:nsec-1L]
    ;Store indices for all non-repetitive combinations
    for j=0L, nsec-1L do begin
      for k=0L, j-1L do begin
        indices = [peaks[i],sec_peaks_i[j],sec_peaks_i[k]]
        indices = indices[sort(indices)]
        constellations_indices[l,*] = indices
        l += 1L
      endfor
    endfor
  endfor
  
  ;Remove repetitions
  ids = strjoin(transpose(strtrim(constellations_indices,2L)),',')
  constellations_indices = constellations_indices[sort(ids),*]
  ids = ids[sort(ids)]
  constellations_indices = constellations_indices[uniq(ids),*]
  ids = ids[uniq(ids)]
  nconst = (size(constellations_indices))[1]
  
  ;Compute coordinates
  constellations_coordinates = fltarr(nconst,5L) + !values.f_nan
  for i=0L, nconst-1L do begin
    ind_i = reform(constellations_indices[i,*])
    ;Center of constellation
    constellations_coordinates[i,0L] = float(abscissa[ind_i[0L]]+abscissa[ind_i[2L]])/2.
    ;Flux of individual lines
    constellations_coordinates[i,1L:3L] = sp[ind_i]
    ;Relative position of the central line
    constellations_coordinates[i,4L] = float(abscissa[ind_i[1L]]-abscissa[ind_i[0L]])/float(abscissa[ind_i[2L]]-abscissa[ind_i[0L]])
  endfor
  
  ;Remove any position with NaNs
  bad = where(~finite(total(constellations_coordinates,2)), nbad, complement=good)
  if nbad eq nconst then message, ' All constellation coordinates contain NaNs !'
  if nbad ne 0L then begin
    constellations_coordinates = constellations_coordinates[good,*]
    constellations_indices = constellations_indices[good,*]
  endif
  
  return, constellations_coordinates
End