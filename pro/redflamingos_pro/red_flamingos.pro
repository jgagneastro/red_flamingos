Pro red_flamingos, RESET=reset, FORCE=force
  
  forward_function gpath, string_replace, folder_check, nbrlist, read_flamingos2, sortself, flamingos2_numtofile, lamp_spectralshape, $
    weighted_median, read_abba_pattern, windows_for_reduc, extract_traces_AB, spectrum_error, tvimage, subtract_calib_sky, offset_spectra, $
    rvb_hex, best_lamp_shift, saveimage, refine_wavelength_solution, strkill, interpol2, supersmooth, fringing_fit_1d, robust_mean, simbad_data, $
    trim, refine_wavelength_solution_relative, folder_check, writespec_ps, align_trace3
  
  ; *** Parameters that need to be set by the user ***
  
  prereduc = 0
  force = 0
  
  ;Should the trace be straightened ? (yes it should)
  straighten_trace = 1
  
  do_screenshots = 1
  screenshots_dir = gpath('redflamingos_screenshots')
  if do_screenshots eq 1 then folder_check, screenshots_dir
  
  ;Are the reduction windows already open ?
  windows_already_open = 0
  
  ;Get path variables
  logfile_dir =  gpath('redflamingos_logdir')
  datadir_base = gpath('redflamingos_raw_data')
  outdir = gpath('redflamingos_reduced_data')
  resources_dir = gpath('redflamingos_idl_resources')
  darks_dir = gpath('redflamingos_darks')
  
  ;Observing log file
  logfile = logfile_dir+'Sample_FLAMINGOSII_Observing_Log.csv'
  
  ;File containing the OSIRIS wavelength-calibrated lamp spectra
  ar_cal_file = resources_dir+'cuar_'+['JH','HK']+'_FLAMINGOS2.sav'
  
  ;Comparison file for a typical extraction (for display purposes only)
  comparison_file = resources_dir+'comparison_'+['JH','HK']+'_FLAMINGOS2.sav'
  
  ;Display the trace detection
  display_diagnostics = 1
  
  ;Are flat fields to be applied ?
  apply_flat_field = 1
  
  ;Are the darks to be applied ? (STRONGLY SUGGESTED, Dark current is enormous in Flamingos II)
  apply_darks = 1
  
  ;Should the spectral structure of the flat lamp be corrected ?
  correct_spectral_structure_lamp = 1
  
  ;Extract & correct tellurics on A and B beams separately ? (this is more tedious, and probably not necessary unless there are some systematics)
  sep_ab_beams = 0
  
  display_wavelength_refinement = 0
  
  read_noise = 11.7 ;e/read
  gain = 4.44 ;e/ADU
  dark_current = 0.5 ;e/s/pixel
  
  ;Use optimal extraction
  optimal = 2
  
  ;Apply a fringing correction ?
  ;0 : No
  ;1 : Telluric only
  ;2 : Science only
  ;3 : Telluric & Science
  do_fringing_correction = 3
  
  ;Should all calibration files be reset ?
  if reset eq !NULL then $
    reset = 0
  
  ;Overwrite existing files by default ?
  if force eq !NULL then $
    force = 0
  
  ;Do a clean-up of remaining tellurics in the A - B trace
  do_corr_telluric_diff = 1
  
  ;Should the wavelength solution be refined with peak detection of Ar lines ?
  do_refine_wavelength_solution = 1
  
  ; *** End of parameters that need to be set by the user ***

  ;If this routine was called by other routines just to get datadir, etc. then hand the information
  ; and return
  if keyword_set(ask_for_params) then begin
    ask_for_params = [logfile, datadir, outdir]
    return
  endif

  ;Avoid math errors
  !except = 0

  ; *** Parameters related to FLAMINGOS-2 ***
  
  ;Number of orders (or number of CCDs)
  norders = 1L
  
  ;0- Start of the horizontal useful region (pixels)
  ;1- End of the horizontal useful region (pixels)
  ;2- Start of the vertical useful region (pixels)
  ;3- End of the vertical useful region (pixels)
  useful_region = [0L,1522L,236L,1755L]
  
  ;If any pixel is more deviant than "nsig_cleanindtrace"-sigma in the individual 2D
  ; spectra, they will be flagged as bad pixels.
  nsig_cleanindtrace = 10.

  ;If any pixel is deviant than "nsig_cleandifftrace"-sigma in the diff (A - B), 2D
  ; spectra, they will be flagged as bad pixels.
  nsig_cleandifftrace = 6.

  ;Pixel width of the vertical (spatial) median filter used to detect bad pixels in 2D spectra
  nvertical_cleantrace = 11L

  ;Pixel width of the 2D median filter used to detect bad pixels in 2D spectra
  n2dfilt_cleantrace = 5L

  ;Pixel width of the horizontal (spectral) median filter used to correct the spectral structure of the flat lamp
  lamp_corr_wid = 64L

  ;Pixel width of the smoothing filter to be applied on the median profile *only* when detecting the trace position
  ; (this is *not* applied on the actual data)
  extract_smooth = 5L

  ;Number of sections in which to break the calibration lamps to apply the shift/strech
  ; for the wavelength solution (at least 2 is needed for a stretch)
  n_cor_sec = 2L

  ;NSIG detection threshold to consider pixels as part of the trace, for trace straightening
  straighten_nsig = 3.

  ;Polynomial degree for trace straightening
  straighten_ndeg = 2

  ;Horizontal fraction of the 2D order not to be considered in trace straightening
  straighten_bordercut = .01

  ;Number of iterations for trace straightening
  straighten_niter = 2L
  
  ;Initial guess of the wavelength solution (polynomial)
  lam_coeff = $;[[1.81492,-0.000662696,6.98771e-09], $ ;JH GRISM
  [[1.87979,-0.000653085,3.01280e-09], $ ;JH GRISM, J band
   [1.88367,-0.000661216,6.79313e-09], $ ;JH GRISM, H band
   [2.46925,-0.000753955,6.44295e-10], $ ;HK GRISM, H band
   [2.47154,-0.000758697,1.51301e-09]]   ; HK GRISM, K band
  
  ;Wavelength (um) positions correspondic to telluric absorption regions,
  ; where residuals are to be masked before fitting a fringing pattern.
  ; This data is not used when no fringing correction is applied.
  nan = !values.f_nan
  fringing_mask = [ [[nan,nan],[nan,nan],[nan,nan]], $;JH Grism
    [[0.,1.4],[1.78,1.95],[2.365,1d3]] ];HK Grism
  
  ;Wavelength (um) regions correspondic to delimitations between distinct
  ; filter configurations, where residuals are to be fitted with a fringing 
  ; pattern. This data is not used when no fringing correction is applied.
  fringing_regions = [ 1.4, $;JH Grism
    1.95 ];HK Grism
  
  ;Now for the regions where to instert a gap in the wavelength solution
  wavelength_sol_regions = [ 1.4, $;JH Grism
    1.82 ];HK Grism
  
  ;Estimated period (in micrometers) for the fringing in each region/filter
  estimated_fringing_period = [.05,$;JH, blue region (NOT TESTED)
    .05,$;JH, red region (NOT TESTED)
    .05,$;HK, blue region (TESTED)
    .04];HK, red region (TESTED)
  
  ;Files containing a correction factor for the slope on the edge of F2 arrays
  ;They must be given for (a) the JH config, then (b) the HK config
  corr_slope_files = '/Users/gagne/Dropbox/Donnees/my_reductions/GEMINI_SOUTH/FLAMINGOS2/slope_corr_flamingos2_'+['JH','HK']+'.sav'
  ;/Users/gagne/Documents/Donnees/my_reductions_old/GEMINI_SOUTH/FLAMINGOS2/slope_corr_flamingos2_HK.sav
  
  ; *** End of parameters related to FLAMINGOS 2 ***

  ;Fix potential problems with input parameters
  if strmid(datadir_base,strlen(datadir_base)-1L) ne path_sep() then datadir_base += path_sep()
  if ~file_test(datadir_base) then message, ' The data directory was not found !'
  if strmid(outdir,strlen(outdir)-1L) ne path_sep() then outdir += path_sep()
  if ~file_test(outdir) then begin
    message, ' The output directory does not exist ! If you enter .continue it will be created. Else type return', /continue
    message, 'Directory that will be created : "'+outdir+'"', /continue
    stop
    folder_check, outdir
  endif
  if ~file_test(logfile) then message, 'The input log file was not found !'

  ;Read logfile
  readcol, logfile, name, dates, fits, program_ids, exp, filter, localtime, airmass, seeing, type, flat_frames, telluric_frames, lamp_frames, $
    darks_pid, darks_date, dark_frames, pattern, reduce, binary, comments, $
    format=strjoin(strarr(19)+'A',','), delimiter=',', /keepspaces, /preserve_null, skipline=2, /silent
  
  ;////Temporary : Avoid binary option/////
  message, ' No binary reduction yet !', /continue
  bad = where(strtrim(binary,2) eq '1', nbad)
  if nbad ne 0L then binary[bad] = '0'
  ;////////////////////////////////////////
  
  bad = where(strtrim(reduce,2) eq '0', nbad)
  if nbad ne 0L then $
    remove, bad, name, dates, fits, program_ids, exp, filter, localtime, airmass, seeing, type, flat_frames, telluric_frames, lamp_frames, $
      darks_pid, darks_date, dark_frames, pattern, reduce, binary, comments
  
  string_replace, name, ';;', ','
  string_replace, fits, ';;', ','
  string_replace, comments, ';;', ','
  string_replace, flat_frames, ';;', ','
  string_replace, telluric_frames, ';;', ','
  string_replace, lamp_frames, ';;', ','
  string_replace, fits, ';', ','
  string_replace, telluric_frames, ';', ','
  string_replace, name, ' ', '_'
  
  telescope = 'Gemini-South'
  instrument = 'Flamingos II'
  
  ;Determine which files are in the science mode.
  to_reduce = where((strlowcase(type) eq 'science' or strlowcase(type) eq 'telluric') and strtrim(reduce,2) eq '1', nred)
  if nred eq 0L then message, 'No files are set for reduction ! They must have the "TYPE" column set to science and the "REDUCE ? (0/1)" column set to 1 !
  
  i0 = 0L
  ;Loop over the observing log entries
  for i=i0, nred-1L do begin
    
    ;Reset some variables
    nx = !NULL
    ny = !NULL
    
    ;Indices of files that need to be reduced
    ii = to_reduce[i]
    
    print, ' > Reducing log entry ['+strtrim(i+1,2)+'/'+strtrim(nred,2)+'] : '+name[ii]+' , '+filter[ii]
    
    ;Choose program ID and date
    program_id = program_ids[ii]
    date = dates[ii]
    
    ;Choose data directory
    datadir = datadir_base+program_id+path_sep()+date+path_sep()
    ;datadir_darks = datadir_base+darks_pid[ii]+path_sep()+darks_date[ii]+path_sep()
    ;if ~file_test(datadir) then $
    ;  message, ' The data directory '+datadir+' does not exist !'
    outdir_i = outdir+program_id+path_sep()+date+path_sep()
    folder_check, outdir_i
    
    ;Determine individual fits numbers
    fits_list_i = nbrlist(fits[ii])
    nfits_i = n_elements(fits_list_i)
    
    ;Determine lists of calibrations
    flats_list_i = nbrlist(flat_frames[ii])
    lamps_list_i = nbrlist(lamp_frames[ii])
    darks_list_i = nbrlist(dark_frames[ii])
    
    if apply_flat_field eq 0L then flats_list_i = !NULL
    if apply_darks eq 0L then darks_list_i = !NULL
    
    ;Read & Create darks
    if darks_list_i eq !NULL then begin
      if apply_darks eq 0L then $
        message, ' The option to use subtract darks was turned off (this is not suggested) !', /continue else $
        message, ' No dark files frame numbers were given for object '+name[ii]+' ! Continuing on, but no darks subtraction will be applied (this is not suggested) !', /continue
      darks = !NULL
      darks_comb_name_i = 'None Applied'
    endif else begin
      darks_comb_name_i = outdir+'DARKS_COMB_'+dark_frames[ii]+'.sav'
      if file_test(darks_comb_name_i) and ~keyword_set(reset) then begin
        restore, darks_comb_name_i;, darks
      endif else begin
        print, '  > Reading & combining darks (this can take a few minutes) ...'
        ;Read and combine dark files
        darks_list_i = file_search(darks_dir+'*')
        nj = n_elements(darks_list_i)
        dark_headers = ptrarr(nj,/allocate)
        exp_times = dblarr(nj)+!values.d_nan
        for j=0L, nj-1L do begin
          ;file = datadir_darks+flamingos2_numtofile(darks_date[ii],darks_list_i[j])
          file = darks_list_i[j]
          if ~file_test(file) then message, ' File "'+file+'" was not found !'
          darkj = read_flamingos2(file,useful_region,HEADER=darkj_header)
          *dark_headers[j] = darkj_header
          nx = (size(darkj))[1]
          ny = (size(darkj))[2]
          if ~keyword_set(dark_cube) then $
            dark_cube = fltarr(nx,ny,nj)+!values.f_nan
          dark_cube[*,*,j] = darkj
          exp_times[j] = sxpar(darkj_header,'EXPTIME')
        endfor
        
        ;Determine unique possible exposure times and combine them together
        u_exp_times = exp_times
        sortself, u_exp_times, /uniq
        nu_exp_times = n_elements(u_exp_times)
        comb_dark_cube = dblarr(nx,ny,nu_exp_times)+!values.d_nan
        for l=0L, nu_exp_times-1L do begin
          good = where(exp_times eq u_exp_times[l], ngood)
          if ngood eq 0L then continue
          comb_dark_cube[*,*,l] = median(dark_cube[*,*,good],dim=3,/even)
        endfor
        
        ;// THIS NEEDS TO BE CHANGED IF EXTRACTION IS DONE WITH MULTIPLE CCDs
        dark_traces = ptrarr(1L,/allocate)
        *dark_traces[0] = comb_dark_cube
        ;//
        
        ;Save the dark
        save, dark_traces, dark_headers, u_exp_times, file=darks_comb_name_i, /compress
        reset = 0
      endelse
    endelse
    ;END : Read & Create Darks
    
    ;Read & Create flats
    if flats_list_i eq !NULL then begin
      if apply_flat_field eq 0L then $
        message, ' The option to use flat field correction was turned off !', /continue else $
        message, ' No flat files frame numbers were given for object '+name[ii]+' ! Continuing on, but no flat field corrections will be applied !', /continue
      flat = !NULL
      flats_comb_name_i = 'None Applied'
    endif else begin
      flats_comb_name_i = outdir_i+'FLATS_COMB_'+flat_frames[ii]+'.sav'
      if file_test(flats_comb_name_i) and ~keyword_set(reset) then begin
        restore, flats_comb_name_i;, flat_traces
      endif else begin
        print, '  > Reading & combining flat files...'
        ;Read and combine flat files
        nj = n_elements(flats_list_i)
        flat_headers = ptrarr(nj,/allocate)
        flat_exp_times = dblarr(nj)+!values.d_nan
        for j=0L, nj-1L do begin
          file = datadir+flamingos2_numtofile(date,flats_list_i[j])
          if ~file_test(file) then message, ' File "'+file+'" was not found !'
          flatj = read_flamingos2(file,useful_region,HEADER=flatj_header)
          *flat_headers[j] = flatj_header
          nx = (size(flatj))[1]
          ny = (size(flatj))[2]
          if ~keyword_set(flat_cube) then $
            flat_cube = fltarr(nx,ny,nj)+!values.f_nan
          flat_cube[*,*,j] = flatj
          flat_exp_times[j] = sxpar(flatj_header,'EXPTIME')
        endfor
        if nj ne 1L then $
          flat = median(flat_cube,dim=3) else $
          flat = flat_cube[*,*,0L]
        
        ;Ensure that all flats have the same exposure time
        if min(flat_exp_times eq median(flat_exp_times),/nan) eq 0L then $
          message, 'All the flat files must have the same exposure time !'
        flat_exp_time = flat_exp_times[0L]
        
        ;Select best-matching dark files and interpolate
        best_dark_pos = interpol(findgen(n_elements(u_exp_times)),u_exp_times,flat_exp_time)
        best_dark_pos >= 0L
        best_dark_pos <= (n_elements(u_exp_times)-1L)
        comb_dark_cube = (*dark_traces[0])
        best_dark = (1. - frac(best_dark_pos)) * comb_dark_cube[*,*,floor(best_dark_pos)] + frac(best_dark_pos) * comb_dark_cube[*,*,ceil(best_dark_pos)]
        
        ;Subtract best dark
        if darks_list_i ne !NULL then $
          flat -= best_dark
        
        ;// THIS NEEDS TO BE CHANGED IF EXTRACTION IS DONE WITH MULTIPLE CCDs
        ;Remove right part of flat because of low SNR
        if filter[ii] eq 'JH' then begin
          x0 = 160L
          x1 = 968L
          y0 = 38L
          y1 = 1518L
          flat[0:x0,*] = 1.
          flat[x1:*,*] = 1.
        endif else begin
          y0 = 55L
          y1 = 1490L
          x0 = 0L
          x1 = nx-1L
        endelse
        flat[*,0:y0] = 1.
        flat[*,y1:*] = 1.
        flat[x0:x1-1L,y0:y1-1L] = lamp_spectralshape(flat[x0:x1-1L,y0:y1-1L], NWID=lamp_corr_wid)
        ;Normalize flat
          flat /= median(flat)
        
        ;Remove large-scale structures
        mf60 = median(flat,60L)
        flat[x0:x1-1L,y0:y1-1L] /= mf60[x0:x1-1L,y0:y1-1L]
        flat /= median(flat)
        
        ;Repair NaNs
        bad = where(~finite(flat),nbad)
        nfilt1 = 10L
        if nbad ne 0 then $
          flat[bad] = (median(flat,nfilt1))[bad]
        ;Repair the remaining NaNs
        bad = where(~finite(flat),nbad)
        if nbad ne 0 then $
          flat[bad] = 1.
        
        ;Store flat in pointer
        flat_traces = ptrarr(1L,/allocate)
        *flat_traces[0] = flat
        ;//
        
        ;Correct the spectral structure of the lamps
        if correct_spectral_structure_lamp eq 1 then begin
          norders = (size(flat_traces))[1]
          for l=0L, norders-1L do begin
            *flat_traces[l] = lamp_spectralshape(*flat_traces[l], NWID=lamp_corr_wid)
            *flat_traces[l] /= median(*flat_traces[l])
          endfor
        endif
        
        ;Save the flat
        save, flat, flat_traces, flat_headers, file=flats_comb_name_i, /compress
      endelse
    endelse
    ;END : Read & Create Flats
    
    ;Read & Create Lamps
    if lamps_list_i eq !NULL then begin
      message, ' No calibration lamp files frame numbers were given for object '+name[ii]+' ! Continuing on, but tellurics will be used instead of Ar lamps !', /continue
      ;Restore the reference tellurics file
      message, 'Not done yet !'
      ;restore, tell_cal_file;, lam_ref, ar_ref
      lamp = !NULL
      lamps_comb_name_i = 'None Applied'
    endif else begin
      ;Restore the reference Ar file
      restore, ar_cal_file[long(filter[ii] eq 'HK')];, lam_ref, cuar_ref
      lam_ref = reform(lam_ref,n_elements(lam_ref))
      cuar_ref = reform(cuar_ref,n_elements(cuar_ref))
      lamps_comb_name_i = outdir_i+'LAMPS_COMB_'+lamp_frames[ii]+'.sav'
      if file_test(lamps_comb_name_i) and ~keyword_set(reset) then begin
        restore, lamps_comb_name_i;, lamp_traces
        if lamp_traces eq !NULL then begin
          message, ' Thrashing corrupted files ...', /continue
          file_move, lamps_comb_name_i, lamps_comb_name_i+'.thrash'
        endif
      endif
      if ~file_test(lamps_comb_name_i) or keyword_set(reset) then begin
        ;Read and combine lamp files
        print, '  > Reading & combining calibration lamp files...'
        nj = n_elements(lamps_list_i)
        lamp_headers = ptrarr(nj,/allocate)
        lamp_exp_times = dblarr(nj)+!values.d_nan
        for j=0L, nj-1L do begin
          file = datadir+flamingos2_numtofile(date,lamps_list_i[j])
          if ~file_test(file) then message, ' The file "'+file+'" was not found !'
          lampj = read_flamingos2(file,useful_region,HEADER=lampj_header)
          *lamp_headers[j] = lampj_header
          if ~keyword_set(nx) then $
            nx = (size(lampj))[1]
          if ~keyword_set(ny) then $
            ny = (size(lampj))[2]
          if ~keyword_set(lamp_cube) then $
            lamp_cube = fltarr(nx,ny,nj)+!values.f_nan
          lamp_cube[*,*,j] = lampj
          lamp_exp_times[j] = sxpar(lampj_header,'EXPTIME')
        endfor
        if nj eq 1L then begin
          lamp = lampj
        endif else begin
          lamp = median(lamp_cube,dim=3)
        endelse
        
        ;Ensure that all lamps have the same exposure time
        if min(lamp_exp_times eq median(lamp_exp_times),/nan) eq 0L then $
          message, 'All the lamps must have the same exposure time !'
        lamp_exp_time = lamp_exp_times[0L]
        
        ;Select best-matching dark files and interpolate
        best_dark_pos = interpol(findgen(n_elements(u_exp_times)),u_exp_times,lamp_exp_time)
        best_dark_pos >= 0L
        best_dark_pos <= (n_elements(u_exp_times)-1L)
        comb_dark_cube = (*dark_traces[0])
        best_dark = (1. - frac(best_dark_pos)) * comb_dark_cube[*,*,floor(best_dark_pos)] + frac(best_dark_pos) * comb_dark_cube[*,*,ceil(best_dark_pos)]
        
        ;Subtract best dark
        if darks_list_i ne !NULL then $
          lamp -= best_dark
        
        ;Divide with flatfield
        if flats_list_i ne !NULL then $
          lamp /= flat
        
        lamp /= weighted_median(lamp,medvalue=.8)
        
        ;// THIS NEEDS TO BE CHANGED IF EXTRACTION IS DONE WITH MULTIPLE CCDs
        lamp_traces = ptrarr(1L,/allocate)
        *lamp_traces[0] = lamp
        ;//
        
        ;Save the lamp
        save, lamp_traces, lamp_headers, file=lamps_comb_name_i, /compress
      endif
    endelse
    
    print, '  > Reading science data...'
    files = strarr(nfits_i)
    hdrs = !NULL
    sci_cube = !NULL
    for j=0L, nfits_i-1L do begin
      
      ;Output files
      filei = outdir+flamingos2_numtofile(date,fits_list_i[j])
      strkill, filei, '.gz'
      strkill, filei, '.fits'
      files[j] = filei
;      if file_test(filei+'.fits') and ~keyword_set(force) then begin
;        print, ' Skipping existing file '+file_basename(files[j])+' !'
;        if strtrim(file_basename(files[j]),2) ne 'S20140601S0147' and $
;          strtrim(file_basename(files[j]),2) ne 'S20140601S0148' then continue
;      endif
      if file_test(filei+'.fits') and ~keyword_set(force) then begin
        print, ' Skipping existing file '+file_basename(files[j])+' !'
        continue
      endif
      
      ;Read data
      file = datadir+flamingos2_numtofile(date,fits_list_i[j])
      if ~file_test(file) then message, ' The file "'+file+'" was not found !'
      imj = read_flamingos2(file,useful_region,HEADER=hdrj)
      
      if ~keyword_set(nx) then begin
        nx = (size(imj))[1]
        ny = (size(imj))[2]
      endif
      
      if ~keyword_set(hdrs) then $
        hdrs = ptrarr(nfits_i,/allocate)
      *hdrs[j] = hdrj
      
      ;Select best-matching dark files and interpolate
      sci_exp_time = sxpar(hdrj,'EXPTIME')
      best_dark_pos = interpol(findgen(n_elements(u_exp_times)),u_exp_times,sci_exp_time)
      best_dark_pos >= 0L
      best_dark_pos <= (n_elements(u_exp_times)-1L)
      comb_dark_cube = (*dark_traces[0])
      best_dark = (1. - frac(best_dark_pos)) * comb_dark_cube[*,*,floor(best_dark_pos)] + frac(best_dark_pos) * comb_dark_cube[*,*,ceil(best_dark_pos)]
      
      ;Apply darks
      if darks_list_i ne !NULL then $
        imj -= best_dark
      
      ;// THIS NEEDS TO BE CHANGED IF EXTRACTION IS DONE WITH MULTIPLE CCDs
      sci_traces = ptrarr(1L,/allocate)
      *sci_traces[0] = imj
      ;//
      
      ;Apply flatfields
      if flats_list_i ne !NULL then $
        for k=0L, norders-1L do $
          *sci_traces[k] /= *flat_traces[k]
      
      ;Clean 2D traces of bad pixels
;      for k=0L, norders-1L do begin
;        im = *sci_traces[k]
;        im = clean_trace(im, NSIG=nsig_cleanindtrace, FILTER_VERTICAL=nvertical_cleantrace, FILTER_2D=n2dfilt_cleantrace,/NOTELL)
;        *sci_traces[k] = im
;      endfor
      
      if ~keyword_set(sx) then begin
        sx = lonarr(norders)-1
        for l=0L, norders-1L do $
          sx[l] = (size(*sci_traces[l]))[1]
      endif
      
      ;Store in a cube
      if ~keyword_set(sci_cube) then $
        sci_cube = ptrarr(nfits_i, norders, /allocate)
      for k=0L, norders-1L do $
        *sci_cube[j,k] = *sci_traces[k]
    endfor
    
    ;Analyze ABBA patterns
    pattern_vec = strtrim(transpose(byte(pattern[ii])),2L)
    if n_elements(pattern_vec) ne nfits_i then $
      message, ' The pattern does not contain the same number of a,b keys as the number of fits files !'
    couples = read_abba_pattern(pattern[ii],repeated)
    n_couples = (size(couples))[2]
    
    ;Create a final pointer that will contain all spectra
    lambda_ext = ptrarr(nfits_i,norders,/allocate)
    spectra_ext = ptrarr(nfits_i,norders,/allocate)
    sky_ext = ptrarr(nfits_i,norders,/allocate)
    calib_ext = ptrarr(nfits_i,norders,/allocate)
    companion_AB = strarr(nfits_i)+'NaN'
    ndeg = (size(lam_coeff))[1]-1
    lam_coeff_corr = dblarr(ndeg+1,norders,nfits_i)+!values.d_nan
    ap_pos = fltarr(nfits_i,norders)+!values.f_nan
    ap_wid = fltarr(nfits_i,norders)+!values.f_nan
    
    ;Extract all AB couples separately
    msg = 0
    for j=0L, n_couples-1L do begin
      ;if min(strtrim(file_basename(files[couples[*,j]]),2) ne 'S20140601S0147') eq 1 and $
      ;  min(strtrim(file_basename(files[couples[*,j]]),2) ne 'S20140601S0148') eq 1 then $
      if min(file_test(files[couples[*,j]]+'.fits')) eq 1 and ~keyword_set(force) then continue
      
      if msg eq 0 then begin
        print, '  > Extracting spectra...'
        msg = 1
      endif
        
      sci_j = sci_cube[reform(couples[*,j]),*]
      ;Extract each order separately
      for k=0L, norders-1L do begin
        ;Compute the "diff" image
        ;diff = *sci_j[0L,k] - *sci_j[1L,k]
        factors = [weighted_median(*sci_j[0L,k],medval=.95),weighted_median(*sci_j[1L,k],medval=.95)]
        factors /= mean(factors)
        diffa = *sci_j[0L,k]/factors[0]
        diffb = *sci_j[1L,k]/factors[1]
        diff = diffa - diffb
        calib_imagei = *lamp_traces[k]
        
        print, '   > Exposures # '+strjoin(strtrim(fits_list_i[reform(couples[*,j])],2),' , ')+' ...'
        
        if ~keyword_set(ny) then $
          ny = (size(diff))[2]
        
        ;stop
        ;lowfreq = bandpass_filter(killnan(diff),0.,0.002)
        ;maxdev = max(abs((median(diff[*,0:500],20))[21:-21,21:-21]))
        ;maxdevhigh = max(abs((median(diff[*,501:*],20))[21:-21,21:-21]))
        ;if maxdev gt 2.5 * maxdevhigh then diff[*,0:500] = !values.f_nan
        
        ;Straighten the trace if it's required
        if straighten_trace eq 1 then $
          align_trace3, diff, calib_imagei, diffa, diffb, NSIGMA=straighten_nsig, NDEG=straighten_ndeg, BORDERCUT=straighten_bordercut, NITER=straighten_niter
        
          ;align_trace2, diff, calib_imagei, NSIGMA=straighten_nsig, NDEG=straighten_ndeg, BORDERCUT=straighten_bordercut, NITER=straighten_niter
        ;if ~keyword_set(okcontinue) then stop
;        diff = *sci_j[0L,k]/factors[0] - *sci_j[1L,k]/factors[1]
;        calib_imagei = *lamp_traces[k]
;        diff[*,0:500] = !values.f_nan
;        straighten_niter=1
;        straighten_ndeg=1
;        align_trace2, diff, calib_imagei, NSIGMA=straighten_nsig, NDEG=straighten_ndeg, BORDERCUT=straighten_bordercut, NITER=straighten_niter
;        diff[*,0:500] = !values.f_nan
        ;okcontinue=1
        
        ;vf, diff
        ;stop
        ;diff = *sci_j[0L,k]/factors[0] - *sci_j[1L,k]/factors[1] & calib_imagei = *lamp_traces[k] & diff[*,0:610] = !values.f_nan  & diff[*,1300:*] = !values.f_nan & align_trace2, diff,ndeg=1,niter=1 & diff[*,1300:*] = !values.f_nan & diff[*,0:610] = !values.f_nan  & vf, diff
        
        ;Clean the diff image of bad pixels
        ;diff = clean_trace(diff, NSIG=nsig_cleandifftrace, FILTER_VERTICAL=nvertical_cleantrace, FILTER_2D=n2dfilt_cleantrace,/NOTELL)
        
        ;Prepare display
        if keyword_set(display_diagnostics) then begin
          if windows_already_open eq 0 then begin
            windows_for_reduc, /MORE
            windows_already_open = 1
          endif
          device,/decomposed
          wset, 0
        endif
        
        ;Extract the raw spectra
;        bad = where((*sci_traces[k])*(*flat_traces[k]) ge 35000d0, nbad)
;        mask = finite(*sci_j[0L,k])*0.+1.
;        if nbad ne 0L then mask[bad] = 0.
;        spectrait = extract_trace3(*sci_j[0L,k]/factors[0], SMOOTH=extract_smooth, SKY=*sci_j[1L,k]/factors[1],GAIN=gain,READNOISE=read_noise)
;        spectrait2 = extract_trace3(diff, GAIN=gain,READNOISE=read_noise)
        spectrai = extract_traces_AB(diff, SMOOTH=extract_smooth, DISPLAY=display_diagnostics, $
          BINARY=(strtrim(binary[ii],2) eq '1'), SKY=skyi, FIT_PARAMS=fit_params, SUBIND=subind, $
          ENVELOPES=envelopes)
        
        ;Do a better subtraction of the telluric lines 
        if do_corr_telluric_diff then begin
          
          ;Combine envelopes to select a proper region of the sky
          env_comb = total(abs(envelopes) / (make_array(ny,value=1.,/float)#total(abs(envelopes),1)), 2, /nan)
          gg = where( env_comb gt max(env_comb,/nan)*.65 )
          yrange = [min(gg,/nan),max(gg,/nan)]
          yrange += [-1,1]*(yrange[1]-yrange[0])*.2
          
          void = max(abs(envelopes[*,0]),gmaxa)
          void = max(abs(envelopes[*,1]),gmaxb)
          box_size = 12L
          min_radius = 12L
           
          ggupa = where(env_comb lt .005 and lindgen(ny) gt (gmaxa+min_radius), nggupa)
          ggupa = ggupa[sort(ggupa)]
          ggupa = ggupa[0:((box_size-1L)<(nggupa-1L))]
          
          ggdowna = where(env_comb lt .005 and lindgen(ny) lt (gmaxa-min_radius), nggdowna)
          ggdowna = ggdowna[reverse(sort(ggdowna))]
          ggdowna = ggdowna[0:((box_size-1L)<(nggdowna-1L))]
          ggdowna = ggdowna[sort(ggdowna)]
          
          ggupb = where(env_comb lt .005 and lindgen(ny) gt (gmaxb+min_radius), nggupb)
          ggupb = ggupb[sort(ggupb)]
          ggupb = ggupb[0:((box_size-1L)<(nggupb-1L))]

          ggdownb = where(env_comb lt .005 and lindgen(ny) lt (gmaxb-min_radius), nggdownb)
          ggdownb = ggdownb[reverse(sort(ggdownb))]
          ggdownb = ggdownb[0:((box_size-1L)<(nggdownb-1L))]
          ggdownb = ggdownb[sort(ggdownb)]
          
          upmeda = median(diffa[*,ggupa],dim=2)
          downmeda = median(diffa[*,ggdowna],dim=2)
          upmedb = median(diffb[*,ggupb],dim=2)
          downmedb = median(diffb[*,ggdownb],dim=2)
          
          upmeda_forimb = median(diffb[*,ggupa],dim=2)
          downmeda_forimb = median(diffb[*,ggdowna],dim=2)
          upmedb_forima = median(diffa[*,ggupb],dim=2)
          downmedb_forima = median(diffa[*,ggdownb],dim=2)
          
          slopea = (upmeda-downmeda) / (mean(double(ggupa),/nan) - mean(double(ggdowna),/nan))
          slopeb = (upmedb-downmedb) / (mean(double(ggupb),/nan) - mean(double(ggdownb),/nan))
          mida = (upmeda + downmeda) / 2d0
          midb = (upmedb + downmedb) / 2d0
          midposa = (mean(double(ggupa),/nan) + mean(double(ggdowna),/nan)) / 2d0
          midposb = (mean(double(ggupb),/nan) + mean(double(ggdownb),/nan)) / 2d0
          
          slopea_forimb = (upmeda_forimb-downmeda_forimb) / (mean(double(ggupa),/nan) - mean(double(ggdowna),/nan))
          slopeb_forima = (upmedb_forima-downmedb_forima) / (mean(double(ggupb),/nan) - mean(double(ggdownb),/nan))
          mida_forimb = (upmeda_forimb + downmeda_forimb) / 2d0
          midb_forima = (upmedb_forima + downmedb_forima) / 2d0
          
          oh_a = (make_array(nx,/double,value=1d0)#(lindgen(ny)-midposa))*(slopea#make_array(ny,/double,value=1d0))+mida#make_array(ny,/double,value=1d0)
          oh_b = (make_array(nx,/double,value=1d0)#(lindgen(ny)-midposb))*(slopeb#make_array(ny,/double,value=1d0))+midb#make_array(ny,/double,value=1d0)
          
          oh_a_forimb = (make_array(nx,/double,value=1d0)#(lindgen(ny)-midposa))*(slopea_forimb#make_array(ny,/double,value=1d0))+mida_forimb#make_array(ny,/double,value=1d0)
          oh_b_forima = (make_array(nx,/double,value=1d0)#(lindgen(ny)-midposb))*(slopeb_forima#make_array(ny,/double,value=1d0))+midb_forima#make_array(ny,/double,value=1d0)
          
          ggapplya = where(abs(lindgen(ny) - gmaxa) lt 2*min_radius and abs(envelopes[*,1]) lt .005, nggapplya)
          ggapplyb = where(abs(lindgen(ny) - gmaxb) lt 2*min_radius and abs(envelopes[*,0]) lt .005, nggapplyb)
          
          ggoutside = lindgen(ny)
          doremove = []
          if nggapplya ne 0L then doremove = [doremove, ggapplya]
          if nggapplya ne 0L then doremove = [doremove, ggapplyb]
          if n_elements(doremove) ne 0 then remove, doremove, ggoutside
          
          ;Clean up diffa, diffb, and regenerate diff
          diffa[*,ggapplya] = diffa[*,ggapplya] - oh_a[*,ggapplya]
          diffa[*,ggapplyb] = diffa[*,ggapplyb] - oh_b_forima[*,ggapplyb]
          diffb[*,ggapplyb] = diffb[*,ggapplyb] - oh_b[*,ggapplyb]
          diffb[*,ggapplya] = diffb[*,ggapplya] - oh_a_forimb[*,ggapplya]
          diffa[*,ggoutside] = !values.d_nan;(diffa[*,ggoutside] - oh_a[*,ggoutside]) / 100.
          diffb[*,ggoutside] = !values.d_nan;(diffb[*,ggoutside] - oh_b[*,ggoutside]) / 100.
          
          ;Clean the diff image of bad pixels
          diff = diffa - diffb
          
          ;diff = clean_trace(diff, NSIG=nsig_cleandifftrace, FILTER_VERTICAL=nvertical_cleantrace, FILTER_2D=n2dfilt_cleantrace,/NOTELL)
          
          
;          RADIUS = 12L
;          gg = where(env_comb lt .005 and (lindgen(ny) lt (yrange[0L]-RADIUS) or lindgen(ny) gt (yrange[1]+RADIUS)), ngg)
;          diffsky = diff[*,gg]
;          
;          ;Clean up tellurics better
;          cleantell = median(diffsky,dim=2)#make_array(ny,value=1.,/float)
;          diff -= cleantell
          
          ;Re-do extraction
          exptimes = [sxpar(*(hdrs[reform(couples[*,j]),*])[0],'EXPTIME'),sxpar(*(hdrs[reform(couples[*,j]),*])[1],'EXPTIME')]
          eff_noise = sqrt(total((exptimes*dark_current)^2) + 2*read_noise^2) / gain
          spectrai = extract_traces_AB(diff, SMOOTH=extract_smooth, DISPLAY=display_diagnostics, $
            BINARY=(strtrim(binary[ii],2) eq '1'), SKY=skyi, FIT_PARAMS=fit_params, SUBIND=subind, $
            ENVELOPES=envelopes, OPTIMAL=optimal, GAIN=gain, READ_NOISE=eff_noise, ERROR=error, $
            EXTRACTION_PROFILES=extraction_profiles,SIGMA_THRESHOLD_OPTIMAL=7.)
          
        endif
        
        ;Re-apply the multiplicative factors
        spectrai[*,0] *= factors[0]
        spectrai[*,1] *= factors[1]
        
        ;Use sky extraction to derive the measurement errors unless optimal extraction was used
        if ~keyword_set(optimal) then optimal = 0
        if optimal eq 2 then begin
          skyi[*,0] = error[*,0]*factors[0]
          skyi[*,1] = error[*,1]*factors[1]
        endif else begin
          skyi[*,0] *= factors[0]
          skyi[*,1] *= factors[1]
          for lli=0L, (size(skyi))[2]-1L do $
            skyi[*,lli] = spectrum_error(findgen(n_elements(skyi[*,lli])), skyi[*,lli])
        endelse
        
        if total(finite(skyi)) eq 0 then $
          message, ' The sky could not be properly extracted ! (See detect_traces_AB.pro)'
        
        ;Display results of extraction
        if keyword_set(display_diagnostics) then begin
          device,/decomposed
          wset, 1
          envelope = total(abs(envelopes),2)
          ndiv = 4L
          nxi = ceil(nx/ndiv)
          bb2 = where(envelope/max(envelope,/nan) ge .01, nbb2)
          imi2 = dblarr(nxi,nbb2*ndiv)+!values.d_nan
          for kk=0L, ndiv-1L do $
            imi2[*,nbb2*kk:nbb2*(kk+1)-1L] = (diff[*,bb2])[nxi*kk:nxi*(kk+1L)-1L,*]
          im0 = fltarr(nxi,nbb2*ndiv,3L)
          im0[*,*,0] = imi2
          im0[*,*,1] = imi2
          im0[*,*,2] = imi2
          im0[*,round((lindgen(ndiv*2-1)+1)*nbb2/2.),1L] = max(imi2,/nan)
          bb = where(~finite(imi2), nbb)
          if nbb ne 0L then begin
            tmp = im0[*,*,0]
            tmp[bb] = max(imi2,/nan)
            im0[*,*,0] = tmp
            tmp = im0[*,*,1]
            tmp[bb] = min(imi2,/nan)
            im0[*,*,1] = tmp
            tmp = im0[*,*,2]
            tmp[bb] = min(imi2,/nan)
            im0[*,*,2] = tmp
          endif
          tvimage,im0,true=3,minvalue=median(im0)-3*robust_sigma(im0,/nan),maxvalue=median(im0)+3*robust_sigma(im0,/nan)
          device,/decomposed
        endif
        
        ;Extract the calibration lamps if they are specified, else use the tellurics
        fit_params[subind-1] = abs(fit_params[subind-1])
        if lamps_list_i ne !NULL then begin
          calibi = extract_traces_AB(calib_imagei, FIT_PARAMS=fit_params, SUBIND=subind[0:1], /POST, EXTRACTION_PROFILES=abs(extraction_profiles))
        endif else begin
          ;Create a "tellurics" image
          calibi = extract_traces_AB(sk_image, FIT_PARAMS=fit_params, SUBIND=subind[0:1], /POST, EXTRACTION_PROFILES=abs(extraction_profiles))
        endelse
        
        makeastop = 1
        ;Shift the calibration to have a minimum at zero
        calibi -= weighted_median(calibi,medval=.02)
        cuar_ref -= weighted_median(cuar_ref,medval=.02)
        
        ;Fix background
        calibi[*,0] = subtract_calib_sky(calibi[*,0])
        calibi[*,1] = subtract_calib_sky(calibi[*,1])
        ;cuar_ref >= 0
        cuar_ref = subtract_calib_sky(cuar_ref)
        
        ;Interpolate the second science spectrum on the first one
        sp2 = skyi[*,1]
        sp3 = calibi[*,1]
        sp1 = offset_spectra(spectrai[*,1], calibi[*,1], calibi[*,0], sp2, sp3, N_COR_SEC=1L, NITER=1)
        spectrai[*,1] = sp1
        skyi[*,1] = sp2
        calibi[*,1] = sp3
        
        ;Remove background structures
        calibi_c = (calibi[*,0] < calibi[*,1])
        
        ;Build an approximate wavelength solution
        xx = findgen(sx[k])+float(useful_region[0L])
        lam_approxk = lam_coeff[0,(filter[ii] eq 'HK')*2] + lam_coeff[1,(filter[ii] eq 'HK')*2] * xx + lam_coeff[2,(filter[ii] eq 'HK')*2] * xx^2
        
        ;Normalize calibrations
        calibi_c /= max(calibi_c,/nan)
        cuar_ref /= max(cuar_ref,/nan)
        
        if keyword_set(display_diagnostics) then begin
          wset, 2
          yrange = [0.,max(spectrai/median(spectrai),/nan)]
          yrange += [0,1]*yrange[1]*.02
          plot, lam_approxk, spectrai[*,0]/median(spectrai), YRANGE=yrange
          oplot, lam_approxk, calibi_c/max(calibi_c,/nan), color=255
          oplot, lam_ref, cuar_ref/max(cuar_ref,/nan),color=rvb_hex(255,100,255)
          oplot, lam_approxk, spectrai[*,0]/median(spectrai)
          oplot, lam_approxk, spectrai[*,1]/median(spectrai), color=rvb_hex(100,255,100)
        endif
        
        ;Do wavelength solution separately within each band
        print, ' Creating a wavelength solution ...'
        separation = wavelength_sol_regions[fix(filter[ii] eq 'HK')]
        lam_bestk = xx * !values.f_nan
        best_coeffs = fltarr((size(lam_coeff))[1],2L)+!values.f_nan
        indfit = 5.
        low_clip = .35
        for soli=0L, 1L do begin
          ccd_ind = where(([1,-1])[soli]*lam_approxk lt ([1,-1])[soli]*separation and lam_approxk lt 2.2 and not (lam_approxk gt 1.77 and lam_approxk lt 1.83), nccd_ind)
          if nccd_ind eq 0L then message, 'Problem separating CCDs'
          ref_ind = where(lam_ref lt 2.2 and not (lam_ref gt 1.77 and lam_ref lt 1.83), nref_ind)
          if nref_ind eq 0L then message, 'Could not find wavelengths < 2.2 microns for the reference Ar spectrum !'
          ccd_ind2 = where(([1,-1])[soli]*lam_approxk lt ([1,-1])[soli]*separation, nccd_ind2)
          if nccd_ind2 eq 0L then message, 'Problem separating CCDs'
          polypar = !NULL
          best_coeffi = best_lamp_shift(lam_ref[ref_ind], cuar_ref[ref_ind]^(1./indfit),(calibi_c[ccd_ind])^(1./indfit), lam_coeff[*,fix(filter[ii] eq 'HK')*2+soli], XFIT=xfit, YFIT=yfit, XVEC=xx[ccd_ind],MONTECARLO=50,LOW_CLIP=low_clip, BAD_REDUCE=5.,refine=do_refine_wavelength_solution, INDFIT=indfit, POLYPAR=polypar)
          if keyword_set(polypar) then $
            xx2 = poly(xx,polypar) else xx2 = xx
          lam_bestk[ccd_ind2] = best_coeffi[0L] + best_coeffi[1L] * xx2[ccd_ind2] + best_coeffi[2L] * xx2[ccd_ind2]^2
          best_coeffs[*,soli] = best_coeffi
        endfor
        best_coeff = mean(best_coeffs,/nan,dim=2)
        
        ;Save the wavelength solution to file
        save, best_coeff, best_coeffs, lam_approxk, lam_bestk, indfit, low_clip, lam_ref, cuar_ref, calibi_c, calibi, file=file_calib
        
        if keyword_set(display_diagnostics) then begin
          wset, 2
          yrange = [0.,max(spectrai/median(spectrai),/nan)]
          yrange += [0,1]*yrange[1]*.02
          plot, lam_bestk, spectrai[*,0]/median(spectrai), YRANGE=yrange
          oplot, lam_bestk, calibi_c/max(calibi_c,/nan), color=255
          oplot, lam_ref, cuar_ref/max(cuar_ref,/nan),color=rvb_hex(255,100,255)
          oplot, lam_bestk, spectrai[*,0]/median(spectrai)
          oplot, lam_bestk, spectrai[*,1]/median(spectrai), color=rvb_hex(100,255,100)
          
          ;Do a zoomed calibration plot
          wset, 3
          indplot = indfit
          plot, (calibi_c/max(calibi_c,/nan))^(1./indplot), lam_bestk, yrange=[([.95,1.82])[fix(filter[ii] eq 'HK')],([1.82,2.28])[fix(filter[ii] eq 'HK')]],ystyle=1,xrange=[0,1.01],xstyle=1
          oplot, (calibi_c/max(calibi_c,/nan))^(1./indplot), lam_bestk, color=255
          oplot, (cuar_ref/max(cuar_ref,/nan))^(1./indplot), lam_ref,color=rvb_hex(100,255,100)
          oplot, (cuar_ref/max(cuar_ref,/nan))^(1./indplot)*0.+low_clip, lam_ref,color=rvb_hex(100,255,255)
          wset, 4
          plot, (calibi_c/max(calibi_c,/nan))^(1./indplot), lam_bestk, yrange=[1.4,1.82],ystyle=1,xrange=[0,1.01],xstyle=1
          oplot, (calibi_c/max(calibi_c,/nan))^(1./indplot), lam_bestk, color=255
          oplot, (cuar_ref/max(cuar_ref,/nan))^(1./indplot), lam_ref,color=rvb_hex(100,255,100)
          oplot, (cuar_ref/max(cuar_ref,/nan))^(1./indplot)*0.+low_clip, lam_ref,color=rvb_hex(100,255,255)
          ;if makeastop eq 1 then stop
          
          if do_screenshots eq 1 then begin
            screenfile = name[ii]+'_'+filter[ii]+'_'+strjoin(strtrim(file_basename(files[reform(couples[*,j])]),2),'_')
            wset,0
            saveimage, screenshots_dir+screenfile+'_DISPLAY0.png', /PNG
            wset,1
            saveimage, screenshots_dir+screenfile+'_DISPLAY1.png', /PNG
            wset,2
            saveimage, screenshots_dir+screenfile+'_DISPLAY2.png', /PNG
            wset,3
            saveimage, screenshots_dir+screenfile+'_DISPLAY3.png', /PNG
            wset,4
            saveimage, screenshots_dir+screenfile+'_DISPLAY4.png', /PNG
          endif
          
          ;Refine wavelength solution in the case of tellurics
          if strlowcase(type[ii]) ne 'science' then begin
            if display_wavelength_refinement then $
              stop
            lam_bestk2 = refine_wavelength_solution(lam_bestk, spectrai[*,0],NODISPLAY=~display_wavelength_refinement)
            lam_bestk_old = lam_bestk
            lam_bestk = lam_bestk2
          endif
          
        endif
        
        ;Save only the un-repeated AB exposures in the output spectra
        tosave = where(repeated[*,j] eq 0, ntosave)
        if ntosave ne 0L then begin
          external_ind = reform(couples[tosave,j])
          for l=0L, ntosave-1L do begin
            *lambda_ext[external_ind[l],k] = lam_bestk
            *spectra_ext[external_ind[l],k] = spectrai[*,tosave[l]]
            *sky_ext[external_ind[l],k] = skyi[*,tosave[l]]
            *calib_ext[external_ind[l],k] = calibi_c;calibi[*,tosave[l]]
            companion_AB[external_ind[l]] = flamingos2_numtofile(date,fits_list_i[couples[long(1-l),j]])
            lam_coeff_corr[*,k,external_ind[l]] = best_coeff
            ;Store the aperture position and width for SpeXTool
            ap_pos[external_ind[l],k] = fit_params[subind[tosave[l]]]
            ap_wid[external_ind[l],k] = fit_params[subind[tosave[l]]+1L]/sqrt(2.*fit_params[subind[tosave[l]]+2L]+2.266065)
          endfor
        endif
      endfor
    endfor
    
    msg = 0
    
    ;Save data in a format that SpeXTools understand
    for j=0L, nfits_i-1L do begin
      if file_test(files[j]+'.fits') and ~keyword_set(force) then continue
      
      if msg eq 0 then begin
        print, '  > Saving raw spectra...'
        msg = 1
      endif
      
      weather = sxpar(*hdrs[j],'RAWIQ',COMMENT=cfocus)+'|'+sxpar(*hdrs[j],'RAWCC',COMMENT=cfocus)+'|'+sxpar(*hdrs[j],'RAWWV',COMMENT=cfocus)+'|'+sxpar(*hdrs[j],'RAWBG',COMMENT=cfocus)
      strkill, weather, '-percentile'
      pixscale = (strtrim(sxpar(*hdrs[j],'FOCUS'),2) eq 'f/16' ? 0.18 : 0.09)
      ctelescop = ''
      cinstrume = ''
      cobserver = ''
      cairmass = ''
      cobsid = ''
      cra = ''
      cdec = ''
      celevatio = ''
      cazimuth = ''
      ccrpa = ''
      cha = ''
      cut = ''
      cdate = ''
      cexptime = ''
      ccoadds = ''
      cgain = ''
      cdetector = ''
      cfilter = ''
      ccamera = ''
      cslit = ''
      cgrating = ''
      cgrattilt = ''
      cmjd = ''
      ccamfocus = ''
      ctelfocus = ''
      cfocus = ''
      cssa = ''
      vals = {$
        origin:'Gemini', $
        telescop:sxpar(*hdrs[j],'TELESCOP',COMMENT=ctelescop), $
        instrume:sxpar(*hdrs[j],'INSTRUME',COMMENT=cinstrume), $
        observer:sxpar(*hdrs[j],'OBSERVER',COMMENT=cobserver), $
        object:strtrim(name[ii],2), $
        airmass:sxpar(*hdrs[j],'AIRMASS',COMMENT=cairmass), $
        prgid:strtrim(program_id,2), $
        obsid:sxpar(*hdrs[j],'OBSID',COMMENT=cobsid), $
        ra:sxpar(*hdrs[j],'RA',COMMENT=cra), $
        dec:sxpar(*hdrs[j],'DEC',COMMENT=cdec), $
        elevatio:sxpar(*hdrs[j],'ELEVATIO',COMMENT=celevatio), $
        azimuth:sxpar(*hdrs[j],'AZIMUTH',COMMENT=cazimuth), $
        crpa:sxpar(*hdrs[j],'CRPA',COMMENT=ccrpa), $
        ha:sxpar(*hdrs[j],'HA',COMMENT=cha), $
        ut:sxpar(*hdrs[j],'TIME-OBS',COMMENT=cut), $
        date:sxpar(*hdrs[j],'DATE-OBS',COMMENT=cdate), $
        exptime:sxpar(*hdrs[j],'EXPTIME',COMMENT=cexptime), $
        coadds:sxpar(*hdrs[j],'COADDS',COMMENT=ccoadds), $
        gain:sxpar(*hdrs[j],'GAIN',COMMENT=cgain), $
        detector:sxpar(*hdrs[j],'DETTYPE',COMMENT=cdetector), $
        filter:sxpar(*hdrs[j],'FILTER',COMMENT=cfilter), $
        camera:sxpar(*hdrs[j],'DETID',COMMENT=ccamera), $
        pixscale:trim(pixscale), $
        slit:sxpar(*hdrs[j],'MASKNAME',COMMENT=cslit), $
        grat:sxpar(*hdrs[j],'GRORDER',COMMENT=cgrating), $
        grattilt:sxpar(*hdrs[j],'SFTILT',COMMENT=cgrattilt), $
        mjd:sxpar(*hdrs[j],'MJD-OBS',COMMENT=cmjd), $
        camfocus:sxpar(*hdrs[j],'FOCUSTEP',COMMENT=ccamfocus), $
        telfocus:sxpar(*hdrs[j],'P2FOCUS',COMMENT=ctelfocus), $
        focus:sxpar(*hdrs[j],'FOCUS',COMMENT=cfocus), $
        log_date:date, $
        weather:weather, $
        TO:sxpar(*hdrs[j],'SSA',COMMENT=cssa), $
        binary:binary[ii], $
        abcmp:companion_AB[j], $
        beam:strupcase(pattern_vec[j]), $
        scif:flamingos2_numtofile(date,fits_list_i[j]), $
        flatf:file_basename(flats_comb_name_i), $
        calf:file_basename(ar_cal_file[long(filter[ii] eq 'HK')]), $
        telf:telluric_frames[ii], $
        lampf:file_basename(lamps_comb_name_i), $
        DISP001:strtrim(abs(best_coeff[1L])), $
        O1MINPX:strjoin(strtrim(useful_region[0L],2),','), $
        O1MAXPX:strjoin(strtrim(useful_region[1L],2),','), $
        O1LC0:strtrim(best_coeff[0L],2), $
        O1LC1:strtrim(best_coeff[1L],2), $
        O1LC2:strtrim(best_coeff[2L],2), $
        CFLTSHP:strtrim(correct_spectral_structure_lamp,2), $
        doflat:strtrim(apply_flat_field,2), $
        NSGA:strtrim(nsig_cleanindtrace,2), $
        NSGAB:strtrim(nsig_cleandifftrace,2), $
        NVRT:strtrim(nvertical_cleantrace,2), $
        N2DFLT:strtrim(n2dfilt_cleantrace,2), $
        LCWID:strtrim(lamp_corr_wid,2), $
        DETSM:strtrim(extract_smooth,2), $
        NCOR:strtrim(n_cor_sec,2), $
        STNSIG:strtrim(straighten_nsig,2), $
        STNDEG:strtrim(straighten_ndeg,2), $
        STBCUT:strtrim(straighten_bordercut,2), $
        STNIT:strtrim(straighten_niter,2), $
        CMT:comments[ii]}
      coms = { $
        origin:'', $
        telescop:' Name of telescope', $
        instrume:cinstrume, $
        observer:cobserver, $
        object:' Name of object', $
        airmass:cairmass, $
        prgid:'Program ID', $
        obsid:cobsid, $
        ra:cra, $
        dec:cdec, $
        elevatio:celevatio, $
        azimuth:cazimuth, $
        crpa:ccrpa, $
        ha:cha, $
        ut:cut, $
        date:cdate, $
        exptime:cexptime, $
        coadds:ccoadds, $
        gain:cgain, $
        detector:cdetector, $
        filter:cfilter, $
        camera:ccamera, $
        pixscale:' Pixel scale on the detector', $
        slit:cslit, $
        grat:cgrating, $
        grattilt:cgrattilt, $
        mjd:cmjd, $
        camfocus:ccamfocus, $
        telfocus:ctelfocus, $
        focus:cfocus, $
        log_date:' Date in the observing log', $
        weather:'IQ|CC|WV|BG-percentile wheather', $
        TO:cssa, $
        binary:' Binary flag (1:binary,0:single)', $
        abcmp:' Companion file for AB subtraction', $
        beam:' Beam position (A or B)', $
        scif:' Raw science file', $
        flatf:' Raw flats files', $
        lampf:' Raw calibration lamp files', $
        calf:' File containing F2 calibrated lamp spectra', $
        telf:' File containing F2 calibrated telluric spectra', $
        DISP001:' Dispersion (um pixel-1) for order 01', $
        O1MINPX:' Pixel position for 2 min. pos. of Order 1', $
        O1MAXPX:' Pixel position for 2 max. pos. of Order 1', $
        O1LC0:' Pix-Lambda 0-order coefficient for Order 1', $
        O1LC1:' Pix-Lambda 1-order coefficient for Order 1', $
        O1LC2:' Pix-Lambda 2-order coefficient for Order 1', $
        CFLTSHP:' 0/1 flag for correction of flat spectral struc.', $
        DOFLAT:' 0/1 flag for application of flat fields', $
        NSGA:' N-sigma filt. to detect bad pix. in indiv. image', $
        NSGAB:' N-sigma fil. to detect bad pix. in A-B image', $
        NVRT:' Pix. wid. of vertical filt. to detect bad pix.', $
        N2DFLT:' Pix. wid. of 2D filter to detect bad pix.', $
        LCWID:' Pix. wid. of horizontal filter to correct flat struc.', $
        DETSM:' Pix. wid. of smooth factor to detect traces',$
        NCOR:' Num. of sec. on which to apply wv corr. on each order', $
        STNSIG:' N-Sig. detec. threshold for trace straigthening', $
        STNDEG:' Polynomial degree for trace straightening', $
        STBCUT:' Fraction of border cut for trace straightening', $
        STNIT:' Number of iterations for trace straightening', $
        CMT:' Comments from the log file'}
      hdr_info = {vals:vals,coms:coms}
      
      nap = 1L
      outname = flamingos2_numtofile(date,fits_list_i[j])
      strkill, outname, '.gz'
      strkill, outname, '.fits'
      
      pix_size = sxpar(*hdrs[j],'MASKNAME')
      strkill, pix_size, 'pix-slit'
      
      slitwarc = float(pix_size) * pixscale
      slitarc = 5.*60.
      slitwpix = long(pix_size)
      slitpix = 5.*60. / pixscale
      lamp_file = file_basename(lamps_comb_name_i,'.sav')
      flat_file = file_basename(flats_comb_name_i,'.sav')

      ;Reform all orders into a fits file
      nxprime = max(sx)
      
      data = fltarr(nxprime, 3L, norders)
      data[*,2L,*] = 1d3
      for l=0L, norders-1L do begin
        ss = sort(*lambda_ext[l])
        data[*,0L,l] = (*lambda_ext[l])[ss]
        data[0:sx[l]-1L,1L,l] = (*spectra_ext[j,l])[ss]
        data[0:sx[l]-1L,2L,l] = (*sky_ext[j,l])[ss]
        if sx[l] lt max(sx) then begin
          data[sx[l]:*,1,l] = weighted_median(*spectra_ext[j,l],medval=.1)
          data[sx[l]:*,2,l] = weighted_median(*sky_ext[j,l],medval=.8)*10
        endif
      endfor

      startpos = 0L
      stoppos = nxprime-1L
      writespec_ps,data,files[j],strupcase(pattern_vec[j]),lamp_file,$
        flat_file,nap,lindgen(norders)+1,startpos,stoppos,hdr_info,[ap_pos[j]],ap_wid[j],strtrim(vals.GRAT,2),$
        slitpix,slitarc,slitwpix,slitwarc,'um','DN / s','!7k!5 (!7l!5m)',$
        'f (!5DN s!u-1!N)','flamingos_extract.pro', /FITS
      
    endfor
    
    ;Use "pattern_vec" to to the right combination and telluric correction
    gAB = where(pattern_vec eq 'a' or pattern_vec eq 'b', ngAB)
;    gA = where(pattern_vec eq 'a', ngA, complement=gB, ncomplement=ngB)
    outname = flamingos2_numtofile(date,fits_list_i[0L])
    
    ;Do the fringing correction if it is needed
    do_fringing_i = 0
    if do_fringing_correction eq 1L and strlowcase(type[ii]) eq 'telluric' then do_fringing_i = 1
    if do_fringing_correction eq 2L and strlowcase(type[ii]) eq 'science' then do_fringing_i = 1
    if do_fringing_correction eq 3L then do_fringing_i = 1
    if do_fringing_i then begin
      
      ;Read the all data files of this target to perform a fringing correction
      files_this_target = outdir+file_basename(flamingos2_numtofile(date,fits_list_i),'.gz')
      
      if min(sxpar_mul(files_this_target,'FRCORR')) ne 1 or keyword_set(force) then begin
        nf = n_elements(files_this_target)
        if ~keyword_set(nx) then $
          nx = (size(readfits(files_this_target[0],/silent)))[1]
        sps = fltarr(nx,nf)+!values.f_nan
        esps = fltarr(nx,nf)+!values.f_nan
        for ll=0L, nf-1L do begin
          im = readfits(files_this_target[ll],hdrll,/silent)
          if ll eq 0L then begin
            lam = im[*,0]
            sps[*,ll] = im[*,1]
            esps[*,ll] = im[*,2]
          endif else begin
            sps[*,ll] = interpol2(im[*,1],im[*,0],lam,/repairnan,badval=!values.f_nan)
            esps[*,ll] = interpol2(im[*,2],im[*,0],lam,/repairnan,badval=!values.f_nan)
          endelse
        endfor
        
        ;residuals = sps/(mean(sps,dim=2)#make_array(nf,value=1.,/float))
        residuals = fltarr(nx,nf)+!values.f_nan
        for ll=0L, nf-1L do $
          residuals[*,ll] = sps[*,ll] / supersmooth(sps[*,ll],64L)
        
        fringing_mask_ii = fringing_mask[*,*,long(filter[ii] eq 'HK')]
        nmask = (size(fringing_mask_ii))[2]
        
        ;Mask telluric regions for fringing fit
        for lli=0L, nmask-1L do begin
          maski = where(lam ge fringing_mask_ii[0L,lli] and $
            lam le fringing_mask_ii[1L,lli], nmaski)
          if nmaski eq 0L then continue
          residuals[maski,*] = !values.f_nan
        endfor
        
        ;Fit fringing solutions
        reconstructed_fringing = dblarr(nx,nf)+1.
        fringing_amplitudes = dblarr(nf,2L)
        for ll=0L, nf-1L do begin
          for llr=0L, 1L do begin
            reg_llr = where(([1,-1])[llr]*lam lt ([1,-1])[llr]*fringing_regions[long(filter[ii] eq 'HK')], nreg_llr)
            yfit = !NULL
            weights = exp(-(findgen(n_elements(reg_llr))-double(n_elements(reg_llr))/2.)^2/(150.^2))
            if filter[ii] eq 'JH' and llr eq 0L then begin
              bad = where(lam[reg_llr] lt 0.94 or lam[reg_llr] gt 1.335, nbad)
              if nbad ne 0L then weights[bad] = 0.
            endif
            if filter[ii] eq 'JH' and llr eq 1L then begin
              bad = where(lam[reg_llr] lt 1.5 or lam[reg_llr] gt 1.75, nbad)
              if nbad ne 0L then weights[bad] = 0.
            endif
            if filter[ii] eq 'HK' and llr eq 0L then begin
              bad = where(lam[reg_llr] lt 1.45 or lam[reg_llr] gt 1.8, nbad)
              if nbad ne 0L then weights[bad] = 0.
            endif
            if filter[ii] eq 'HK' and llr eq 1L then begin
              bad = where(lam[reg_llr] lt 2.05 or lam[reg_llr] gt 2.36, nbad)
              if nbad ne 0L then weights[bad] = 0.
            endif
            estimated_periods = [.04,.045,.05,.055,.06]
            ;La vraie quantité physique c'est 2 * cette periode, en unites de microns
            nper = n_elements(estimated_periods)
            chis = fltarr(nper)+!values.d_nan
            for kk=0L, nper-1L do begin
              estimated_period = estimated_periods[kk]
              par_fringing_kk = fringing_fit_1d(lam[reg_llr], residuals[reg_llr,ll], ESTIMATED_PERIOD=estimated_period, $
                  YFIT=yfit, ERR=err, /SINUS,WEIGHTS=weights)
              chis[kk] = robust_mean((residuals[reg_llr,ll]/yfit-1.)^2)
            endfor
            void = min(chis,/nan,wm)
            estimated_period = estimated_periods[wm]
            par_fringing = fringing_fit_1d(lam[reg_llr], residuals[reg_llr,ll], ESTIMATED_PERIOD=estimated_period, $
              YFIT=yfit, ERR=err, /SINUS, weights=weights)
            ;if name[ii] eq 'HIP 13917' then stop
            if ~keyword_set(yfit) or n_elements(yfit) eq 1 or err ne '' then continue
            ;Only keep the fit if it decreased median-chi^2
            if robust_mean((residuals[reg_llr,ll]/yfit-1.)^2)/robust_mean((residuals[reg_llr,ll]-1.)^2) ge 1. then continue
            ;Only keep the fit if fringing period makes sense (within 25% of estimate)
            if abs((!dpi/par_fringing[1]))*2 lt .03 or abs((!dpi/par_fringing[1]))*2 gt .065 then continue
            reconstructed_fringing[reg_llr,ll] = yfit
            fringing_amplitudes[ll,llr] = abs(par_fringing[0L])
          endfor
        endfor
        if display_diagnostics then begin
          wset, 0
          device, decomposed=1
          del = 0.3
          for ll=0L, nf-1L do begin & $
            if ll eq 0L then $
              plot, residuals[*,ll] - median(residuals[*,ll]) + del*(ll+1), yrange=[0.,(nf+1)*del] else $
              oplot, residuals[*,ll] - median(residuals[*,ll]) + del*(ll+1) & $
            oplot, reconstructed_fringing[*,ll] - 1 + del*(ll+1), color=255 & $
            oplot, residuals[*,ll]/reconstructed_fringing[*,ll] - median(residuals[*,ll]/reconstructed_fringing[*,ll]) + del*(ll+1) + del/2., color=rvb_hex(100,255,100) & $
          endfor
        endif
        
        ;Smooth out fringing in the telluric regions
        lli=1L
        maski = where(lam ge fringing_mask_ii[0L,lli] and $
          lam le fringing_mask_ii[1L,lli], nmaski)
        if nmaski ne 0L then begin
          for ll=0L, nf-1L do begin
            if fringing_amplitudes[ll,0] lt 1d-5 then continue
            smoothening_vector = dindgen(nmaski)/double(nmaski-1L) * (fringing_amplitudes[ll,1]-fringing_amplitudes[ll,0]) + fringing_amplitudes[ll,0]
            smoothening_vector /= fringing_amplitudes[ll,0]
            reconstructed_fringing[maski,ll] = (reconstructed_fringing[maski,ll]-1.)*smoothening_vector+1.
          endfor
        endif
        
        ;Apply this correction to the spectra
        for ll=0L, nf-1L do begin
          im = readfits(files_this_target[ll],/silent,hdr)
          ;Verify that it was not already corrected
          if sxpar(hdr,'FRCORR') then continue
          sxaddpar, hdr, 'FRCORR', 1, 'Fringing was corrected with fringing_fit_1d.pro'
          ;Apply the fringing correction on spectrum and error
          fringing_interpol = interpol2(reconstructed_fringing[*,ll],lam,im[*,0],badval=1.)
          im[*,1] /= fringing_interpol
          im[*,2] /= fringing_interpol
          writefits, files_this_target[ll], im, hdr, /silent
        endfor
      endif
    endif
    
    if prereduc eq 1 then $ 
      continue
    
    ;Combine raw extractions
    if sep_ab_beams eq 0 then begin
      
      ;Combine all files (both the A and B beams)
      save_file_ab = name[ii]+'_BEAMS_AB_'+date+'_'+strjoin(strtrim(fits_list_i[gab],2),'_')+'_comb'
      string_replace, save_file_ab, ' ', '_'
      if ~file_test(outdir+save_file_ab+'.fits') or keyword_set(force) then begin
        ;If there is only one file, don't combine
        if ngab eq 1L then begin
          file_copy, outdir+flamingos2_numtofile(date,fits_list_i[gab[0L]],/nogz), outdir+save_file_ab[0L]+'.fits', /overwrite
        endif else begin
          xcombspec,'FLAMINGOS2.dat',PATH_INPUT=outdir,PREFIX_INPUT=strmid(outname,0,10),FILES_INPUT=strjoin(strtrim(fits_list_i[gab],2),','),$
            SAVE_INPUT=save_file_ab
          print, ' Enter .continue when done !'
          stop
        endelse
      endif
    
    endif else begin
      
      ;Combine files of the A beam
      save_file_a = name[ii]+'_BEAM_A_'+date+'_'+strjoin(strtrim(fits_list_i[ga],2),'_')+'_comb'
      string_replace, save_file_a, ' ', '_'
      if ~file_test(outdir+save_file_a+'.fits') or keyword_set(force) then begin
        ;If there is only one file, don't combine
        if nga eq 1L then begin
          file_copy, outdir+flamingos2_numtofile(date,fits_list_i[ga[0L]],/nogz), outdir+save_file_a[0L]+'.fits', /overwrite
        endif else begin
          xcombspec,'FLAMINGOS2.dat',PATH_INPUT=outdir,PREFIX_INPUT=strmid(outname,0,10),FILES_INPUT=strjoin(strtrim(fits_list_i[ga],2),','),$
            SAVE_INPUT=save_file_a
          print, ' Enter .continue when done !'
          stop
        endelse
      endif
      
      ;Combine files of the B beam
      save_file_b = name[ii]+'_BEAM_B_'+date+'_'+strjoin(strtrim(fits_list_i[gb],2),'_')+'_comb'
      string_replace, save_file_b, ' ', '_'
      if ~file_test(outdir+save_file_b+'.fits') or keyword_set(force) then begin
        ;If there is only one file, don't combine
        if ngb eq 1L then begin
          file_copy, outdir+flamingos2_numtofile(date,fits_list_i[gb[0L]],/nogz), outdir+save_file_b[0L]+'.fits', /overwrite
        endif else begin
          xcombspec,'FLAMINGOS2.dat',PATH_INPUT=outdir,PREFIX_INPUT=strmid(outname,0,10),FILES_INPUT=strjoin(strtrim(fits_list_i[gb],2),','),$
            SAVE_INPUT=save_file_b
          print, ' Enter .continue when done !'
          stop
        endelse
      endif
      
    endelse
    
    ;If it's a science file, apply telluric correction
    if strlowcase(type[ii]) eq 'science' then begin
      
      gtell = where(dates eq dates[ii] and fits eq telluric_frames[ii], ngtell)
      if ngtell ne 1L then $
        message, ' Ambiguous definition of telluric files to use !' 
      pattern_vec_tell = strtrim(transpose(byte(pattern[gtell])),2L)
      name_tell = name[gtell[0L]]
      string_replace, name_tell, ' ', '_'
      
      if sep_ab_beams eq 0 then begin
        gtellab = where(pattern_vec_tell eq 'a' or pattern_vec_tell eq 'b', ngtellab)
        if ngab ne 0L and ngtellab eq 0L then $
          message, ' There is no A-beam or B-beam telluric observation defined for this target !'
      endif else begin
        gtella = where(pattern_vec_tell eq 'a', ngtella)
        gtellb = where(pattern_vec_tell eq 'b', ngtellb)
        if nga ne 0L and ngtella eq 0L then $
          message, ' There is no A-beam telluric observation defined for this target !'
        if ngb ne 0L and ngtellb eq 0L then $
          message, ' There is no B-beam telluric observation defined for this target !'
      endelse
      
      ;Apply telluric corrections
      if sep_ab_beams eq 0 then begin
        if ngab ne 0L then begin
          tell_file_ab = name_tell+'_BEAMS_AB_'+date+'_'+strjoin(strtrim((nbrlist(telluric_frames[ii]))[gtellab],2),'_')+'_comb'
          if ~file_test(outdir+tell_file_ab+'.fits') then $
            message, ' Some combined telluric files are missing ! Always place telluric files *before* the associated science files in the reduction logs !'
          save_tellcor_ab = save_file_ab
          string_replace, save_tellcor_ab, '_comb', '_tc'
          
          if ~file_test(outdir+save_tellcor_ab+'.fits') or keyword_set(force) then begin
            hdr = headfits(outdir+save_file_ab+'.fits',/silent)
            if n_elements(hdr) le 20 then $
              message, ' The header for this image is suspiciously short. A problem probably happened in the previous steps !'
            
            ngood = 0L
            if strpos(strlowcase(name[gtell[0L]]),'hip') ne -1 then begin
              if ~keyword_set(hipparcos) then begin
                restore, resources_dir+'hipparcos.sav', /ver
                hipparcos = a
              endif
              good = where('hip '+strtrim(hipparcos.HIP,2) eq strlowcase(name[gtell[0]]), ngood)
            endif
            if ngood ne 0L then begin
              bmag = hipparcos[good[0]].BTMAG
              vmag = hipparcos[good[0]].VTMAG
            endif else begin
              sim = simbad_data(name[gtell[0L]],/US)
              if ~isa(sim,'struct') then $
                message, ' Could not download SIMBAD data to get B and V magnitudes for the telluric !'
              if ~finite(sim.bmag*sim.vmag) then $
                message, ' No SIMBAD information on B and V magnitudes for the telluric !'
              bmag = trim(sim.bmag)
              vmag = trim(sim.vmag)
            endelse
            ;Start by refining the wavelength solution of the science target, using the standard and comparing
            ; the relative positions of telluric absorption lines  
            im_ref = readfits(outdir+tell_file_ab+'.fits', hdr_ref,/silent)
            im_obs = readfits(outdir+save_file_ab+'.fits', hdr_obs,/silent)
            if filter[ii] eq 'JH' then $
              use_lines = [1.133,1.268,1.394,1.431,1.47] else $
              use_lines = [1.394,1.431,1.47,1.82125,1.9205,1.957,2.005,2.368,2.4315]
            lam_obs = refine_wavelength_solution_relative(im_obs[*,0], im_obs[*,1], im_ref[*,0], im_ref[*,1], POLYPAR=polypar,USE_LINES=use_lines,NODISPLAY=~display_wavelength_refinement,WID=14L,/REMOVECONT, /SLOPE)
            ;lam_obs = refine_wavelength_solution_relative(im_obs[*,0], im_obs[*,1]/savitzky_golay(im_obs[*,1],80), im_ref[*,0], im_ref[*,1]/savitzky_golay(im_ref[*,1],80), POLYPAR=polypar,USE_LINES=use_lines,NODISPLAY=~display_wavelength_refinement,WID=12L)
            im_obs[*,0] = lam_obs
            sxaddhist, 'Wavelength solution was refined with "refine_wavelength_solution_relative.pro"', hdr_obs
            sxaddhist, '  using the telluric standard "'+tell_file_ab+'.fits"', hdr_obs
            sxaddhist, '  and its telluric absorption lines', hdr_obs
            writefits, outdir+save_file_ab+'.fits',im_obs, hdr_obs 
            
            ;Fix wavelength solution with the HI lines and the telluric standard
            xtellcor_general,PATH_INPUT=outdir,SCIENCE_INPUT=save_file_ab+'.fits',TELL_INPUT=tell_file_ab+'.fits', B_INPUT=bmag, V_INPUT=vmag, $
              SAVE_INPUT=save_tellcor_ab, FWHM=float(sxpar(hdr,'SLTW_PIX'))*abs(float(sxpar(hdr,'O1LC1'))), SPT_INPUT='A0';strmid(numtospectral(round(spectraltonum(telluric_frames[gtell[0]]))),0,2)
            print, ' Enter .continue when done !'
            stop
            ;Put back the header, which is not preserved by "xtellcor_general.pro"
            imi = readfits(outdir+save_tellcor_ab+'.fits',hdri,/silent)
            sxaddhist, 'Telluric correction was performed with xtellcor_general.pro', hdr
            writefits, outdir+save_tellcor_ab+'.fits', imi, hdr
          endif
        endif
      endif else begin
        if nga ne 0L then begin
          tell_file_a = name_tell+'_BEAM_A_'+date+'_'+strjoin(strtrim((nbrlist(telluric_frames[ii]))[gtella],2),'_')+'_comb'
          if ~file_test(outdir+tell_file_a+'.fits') then $
            message, ' Some combined telluric files are missing ! Always place telluric files *before* the associated science files in the reduction logs !'
          save_tellcor_a = save_file_a
          string_replace, save_tellcor_a, '_comb', '_tc'
          
          if ~file_test(outdir+save_tellcor_a+'.fits') or keyword_set(force) then begin
            hdr = headfits(outdir+save_file_a+'.fits',/silent)
            if n_elements(hdr) le 20 then $
              message, ' The header for this image is suspiciously short. A problem probably happened in the previous steps !'
            xtellcor_general,PATH_INPUT=outdir,SCIENCE_INPUT=save_file_a+'.fits',TELL_INPUT=tell_file_a+'.fits', B_INPUT=bmag, V_INPUT=vmag, $
                SAVE_INPUT=save_tellcor_a, FWHM=float(sxpar(hdr,'SLTW_PIX'))*abs(float(sxpar(hdr,'O1LC1')))
            print, ' Enter .continue when done !'
            stop
            ;Put back the header, which is not preserved by "xtellcor_general.pro"
            imi = readfits(outdir+save_tellcor_a+'.fits',hdri,/silent)
            sxaddhist, 'Telluric correction was performed with xtellcor_general.pro', hdr
            writefits, outdir+save_tellcor_a+'.fits', imi, hdr
          endif
        endif
        
        if ngb ne 0L then begin
          tell_file_b = name_tell+'_BEAM_B_'+date+'_'+strjoin(strtrim((nbrlist(telluric_frames[ii]))[gtellb],2),'_')+'_comb'
          if ~file_test(outdir+tell_file_b+'.fits') then $
            message, ' Some combined telluric files are missing ! Always place telluric files *before* the associated science files in the reduction logs !'
          save_tellcor_b = save_file_b
          string_replace, save_tellcor_b, '_comb', '_tc'
          if ~file_test(outdir+save_tellcor_b+'.fits') or keyword_set(force) then begin
            hdr = headfits(outdir+save_file_b+'.fits',/silent)
            if n_elements(hdr) le 20 then $
              message, ' The header for this image is suspiciously short. A problem probably happened in the previous steps !'
            xtellcor_general,PATH_INPUT=outdir,SCIENCE_INPUT=save_file_b+'.fits',TELL_INPUT=tell_file_b+'.fits', B_INPUT=bmag, V_INPUT=vmag, $
              SAVE_INPUT=save_tellcor_b, FWHM=float(sxpar(hdr,'SLTW_PIX'))*abs(float(sxpar(hdr,'O1LC1')))
            print, ' Enter .continue when done !'
            stop
            ;Put back the header, which is not preserved by "xtellcor_general.pro"
            imi = readfits(outdir+save_tellcor_b+'.fits',hdri,/silent)
            sxaddhist, 'Telluric correction was performed with xtellcor_general.pro', hdr
            writefits, outdir+save_tellcor_b+'.fits', imi, hdr
          endif
        endif
      endelse
    endif
    
    ;If this is a science object and the last observation of this target, then combine all data logs together.
    if strlowcase(type[ii]) eq 'science' and max(where(name eq name[ii])) eq ii then begin
      
      science_to_reduce = where(name eq name[ii])
      names = name[science_to_reduce]
      filters = filter[science_to_reduce]
      id = names+'_'+filters
      uniq_id = id
      sortself, uniq_id, /uniq
      nu = n_elements(uniq_id)
      
      ;First, combine identical filters together
      uniq_names = []
      uniq_filters = []
      uniq_filenames = []
      for ji=0L, nu-1L do begin
        ggi = where(id eq uniq_id[ji], nggi)
        if nggi eq 0L then continue
        out_file_i = name[science_to_reduce[ggi[0L]]]+'_ALLCOMB_'+filter[science_to_reduce[ggi[0L]]]+'_tc'
        
        if sep_ab_beams eq 0 then begin
          filenames_ab = []
          for kki=0L, nggi-1L do begin
            pattern_vec = strtrim(transpose(byte(pattern[science_to_reduce[ggi[kki]]])),2L)
            gab = where(pattern_vec eq 'a' or pattern_vec eq 'b', ngab)
            if ngab ne 0L then begin
              fits_ab = strjoin(strtrim((nbrlist(fits[science_to_reduce[ggi[kki]]]))[gab],2L),'_')
              filenames_ab = [filenames_ab, name[science_to_reduce[ggi[kki]]]+'_BEAMS_AB_'+dates[science_to_reduce[ggi[kki]]]+'_'+fits_ab+'_tc.fits']
            endif
          endfor
          all_files_i = filenames_AB
        endif else begin
          filenames_a = []
          filenames_b = []
          for kki=0L, nggi-1L do begin
            pattern_vec = strtrim(transpose(byte(pattern[science_to_reduce[ggi[kki]]])),2L)
            ga = where(pattern_vec eq 'a', nga)
            if nga ne 0L then begin
              fits_a = strjoin(strtrim((nbrlist(fits[science_to_reduce[ggi[kki]]]))[ga],2L),'_')
              filenames_a = [filenames_a, name[science_to_reduce[ggi[kki]]]+'_BEAM_A_'+dates[science_to_reduce[ggi[kki]]]+'_'+fits_a+'_tc.fits']
            endif
            gb = where(pattern_vec eq 'b', ngb)
            if ngb ne 0L then begin
              fits_b = strjoin(strtrim((nbrlist(fits[science_to_reduce[ggi[kki]]]))[gb],2L),'_')
              filenames_b = [filenames_b, name[science_to_reduce[ggi[kki]]]+'_BEAM_B_'+dates[science_to_reduce[ggi[kki]]]+'_'+fits_b+'_tc.fits']
            endif
          endfor
          all_files_i = [filenames_A,filenames_B]
        endelse
        
        ;Sanity check on file headers
        for kik=0L, n_elements(all_files_i)-1L do begin
          file_kik = outdir+all_files_i[kik]
          if ~file_test(file_kik) then $
            message, ' File "'+file_kik+'" was not found !'
          if n_elements(headfits(file_kik,/silent)) le 20 then begin
            infile_hdr = file_basename(file_kik)
            string_replace, infile_hdr, '_tc', '_comb'
            hdr = headfits(outdir+infile_hdr,/silent)
            if n_elements(hdr) le 20 then $
              message, ' There seems to be a problem - the header of this file is suspiciously short !'
            imi = readfits(file_kik,hdri,/silent)
            sxaddhist, 'Telluric correction was performed with xtellcor_general.pro', hdr
            writefits, file_kik, imi, hdr
          endif
        endfor
        
        if ~file_test(outdir+out_file_i+'.fits') or keyword_set(force) then begin
          ;If there is only one file, don't combine
          nfi = n_elements(all_files_i)
          if nfi eq 1L then begin
            file_copy, outdir+all_files_i[0], outdir+out_file_i[0]+'.fits', /overwrite
          endif else begin
            xcombspec,'FLAMINGOS2.dat',PATH_INPUT=outdir,FILES_INPUT=strjoin(all_files_i,','),$
              SAVE_INPUT=out_file_i, /FULL_FILE_NAMES
            print, ' Enter .continue when done !'
            stop
          endelse
        endif
        
        ;ADD TO UNIQ_NAMES AND UNIQ_FILTERS
        uniq_names = [uniq_names, name[science_to_reduce[ggi[0L]]]]
        uniq_filters = [uniq_filters, filter[science_to_reduce[ggi[0L]]]]
        uniq_filenames = [uniq_filenames, out_file_i]
      endfor
      
      ;Then, combine JH & HK filters together
      id = uniq_names
      uniq_id = id
      nid = n_elements(id)
      sortself, uniq_id, /uniq
      nu = n_elements(uniq_id)

      ;First, change file headers so that xspextool understands that it now has 2 apertures
      ; that must be combined
      for ji=0L, nu-1L do begin
        ggi = where(id eq uniq_id[ji], nggi)
        for j=0L, nggi-1L do begin
          file_in = uniq_filenames[ggi[j]]+'.fits'
          file_out = file_in
          strkill, file_out, ['_JH','_HK']
          im = readfits(outdir+file_in,hdr,/silent)
          if n_elements(hdr) le 20 then $
            message, ' The header for this image is suspiciously short. A problem probably happened in the previous steps !'
          if j eq 0L then begin
            imcube = dblarr((size(im))[1],(size(im))[2],nggi)
            par = sxpar(hdr,'NORDERS',COMMENT=cmt)
            sxaddpar, hdr, 'NORDERS', nggi, cmt
            par = sxpar(hdr,'ORDERS',COMMENT=cmt)
            sxaddpar, hdr, 'ORDERS', strjoin(strtrim(lindgen(nggi)+1L,2L),','), cmt
            hdrcube = hdr
            apposs = fltarr(1L,nggi)
          endif else begin
            par = sxpar(hdr,'O1LC0',COMMENT=cmt)
            sxaddpar, hdrcube, 'O'+strtrim(j+1L,2L)+'LC0', par, cmt
            par = sxpar(hdr,'O1LC1',COMMENT=cmt)
            sxaddpar, hdrcube, 'O'+strtrim(j+1L,2L)+'LC1', par, cmt
            par = sxpar(hdr,'O1LC2',COMMENT=cmt)
            sxaddpar, hdrcube, 'O'+strtrim(j+1L,2L)+'LC2', par, cmt
            par = sxpar(hdr,'DISP001',COMMENT=cmt)
            sxaddpar, hdrcube, 'DISP00'+strtrim(j+1L,2L), par, cmt
            par = sxpar(hdr,'APPOSO01',COMMENT=cmt)
            apposs[0L,j] = float(par)
            sxaddpar, hdrcube, 'APPOSO0'+strtrim(j+1L,2L), par, cmt
          endelse
          imcube[*,*,j] = im
        endfor

        ;Apply correction factor for the brighstar leaks near array borders
        if nid gt 1 then begin
          if file_test(corr_slope_files[0]) then begin
            restore, corr_slope_files[0]
            corr = !NULL
            corr = interpol2(slope_corr_JH,lam_corr_JH,imcube[*,0,1],badval=1.)
            imcube[*,1,1] *= corr
            imcube[*,2,1] *= corr
          endif else message, 'The data files for slope corrections were not found ! The corrections will thus not be applied !', /continue
          if file_test(corr_slope_files[1]) then begin
            restore, corr_slope_files[1]
            corr = !NULL
            corr = interpol2(slope_corr_HK,lam_corr_HK,imcube[*,0,0],badval=1.)
            imcube[*,1,0] *= corr
            imcube[*,2,0] *= corr
          endif else message, 'The data files for slope corrections were not found ! The corrections will thus not be applied !', /continue
        endif else begin
          if min(file_test(corr_slope_files)) eq 1 then begin
            restore, corr_slope_files[long(filters[0] eq 'HK')]
            if filters[0] eq 'HK' then $
              corr = interpol2(slope_corr_HK,lam_corr_HK,imcube[*,0,0],badval=1.) else $
              corr = interpol2(slope_corr_JH,lam_corr_JH,imcube[*,0,0],badval=1.)
            imcube[*,1] *= corr
            imcube[*,2] *= corr
          endif else $
            message, 'The data files for slope corrections were not found ! The corrections will thus not be applied !', /continue
        endelse
        ;Save file for the combination
        writespec_ps,imcube,outdir+file_basename(file_out,'.fits'),'A',sxpar(hdrcube,'LAMPF'),$
          sxpar(hdrcube,'FLATF'),sxpar(hdrcube,'NAPS'),lindgen(long(sxpar(hdrcube,'NORDERS')))+1,sxpar(hdrcube,'START'),sxpar(hdrcube,'STOP'),$
          replicate(gethdrinfo(hdrcube,!NULL),1L),apposs,float(sxpar(hdrcube,'APRADIUS')),strtrim(sxpar(hdrcube,'GRAT'),2),$
          sxpar(hdrcube,'SLTH_PIX'),sxpar(hdrcube,'SLTH_ARC'),sxpar(hdrcube,'SLTW_PIX'),sxpar(hdrcube,'SLTW_ARC'),'um','DN / s','!7k!5 (!7l!5m)',$
          'f (!5DN s!u-1!N)','flamingos_extract.pro', /FITS
        
        file_merge = file_out
        string_replace, file_merge, '_tc.fits', '_merge'
        
        if nid eq 1 then begin
          file_merge = file_out
          strkill, file_merge, '.fits'
        endif
        
        if ~file_test(outdir+'MERGE/'+file_merge+'.fits') or keyword_set(force) then begin
          false_pos = 0L
          
          if nid gt 1 then begin
            xmergeorders, SCIENCE_INPUT=outdir+file_out, $
              SAVE_INPUT=file_merge
            stop
            print, ' Enter .continue when done !'
          endif
          
          ;Move the output file
          folder_check, outdir+'MERGE/'
          file_move, outdir+file_merge+'.fits', outdir+'MERGE/'+file_merge+'.fits', /OVERWRITE
          
          ;Do basic analysis
          if false_pos eq 0L then begin
;            spt_compare, outdir+'MERGE/'+file_merge+'.fits', spt_cen='M8', /big
;  
;            help, allers_index(outdir+'MERGE/'+file_merge+'.fits')
;            print,name[ii]
;            stop
;            ;spt_compare_allknown,outdir+'MERGE/'+file_merge+'.fits',spt='L1',graph=graph,nsmooth=10,ngolay=20,leg_xshift=.5,leg_yshift=-.08,mask=[2.4,2.5],tag='J2235-3844 (L1$\gamma$), J-K$_S$ = 2.15' & graph.yrange=[0.2,1.153]
;            ;graph.save,'J2235-3844_comp_all.png'
            
          endif else begin
            plot3, outdir+'MERGE/'+file_merge+'.fits'
            stop
          endelse
          
          print, ' Enter .continue when done !'
        endif
        
      endfor
    endif
  endfor
  
  print, 'All Done !'
  stop
  
End