Function detect_traces_AB, diff_in, SMOOTH=smooth, DISPLAY=display, BINARY=binary, $
  PEAK_VALUES=peak_values, FIT_PARAMS=fit_params, STOP=stop, SKY_IND=sky_ind, POST=post, $
  PEAK_WIDTHS=peak_widths, PEAK_INDICES=peak_indices, SUBIND=subind, UNRESOLVED=unresolved, $
  ONLY_AMPLITUDE=only_amplitude, XRANGE=xrange, FIXED_SEP=fixed_sep, TITLE=title, DIFFS=diffs
  ;On doit fournir une pose A-B en entree
  ;SMOOTH donne la taille du filtre pour adoucir la trace. 0 = pas de smooth.
  ;DISPLAY sert a afficher les etapes. Mettre un string pour s'en servir comme
  ; titre du plot.
  ;BINARY indique que chaque trace est binaire. Mettre = 1 pour
  
  forward_function remove, mpfitfun, rvb_hex, robust_sigma, AB_moffat_fixed_sep, tvread, saveimage
  
  ;On evite d'alterer la trace en entree
  diff = diff_in
  ny = (size(diff))[2]
  
  ;On fait une mediane verticale pour aller detecter la trace
  if keyword_set(diffs) then begin
    nd = (size(diffs))[3]
    dd = []
    for jj=0L, nd-1L do $
      dd = [dd,diffs[*,*,jj]]
    trace = median(dd,dim=1)
  endif else $
    trace = median(diff,dim=1)
  bad = where(~finite(trace), nbad)
  if nbad ne 0 then trace[bad] = 0.
  
  ;On adoucit la trace
  if keyword_set(smooth) then $
    smtrace = smooth(trace,smooth) else $
    smtrace = trace
  
  ;Remove sky
  diff -= median(smtrace)
  trace -= median(smtrace)
  smtrace -= median(smtrace)
  
  ;On determine la position des deux pires extrema ainsi que leur valeur
  nshift = 3L
  if keyword_set(unresolved) then nshift = 1L
  imax = where(smtrace[nshift:ny-nshift-1L] ge (shift(smtrace,nshift))[nshift:ny-nshift-1L] and smtrace[nshift:ny-nshift-1L] ge (shift(smtrace,-nshift))[nshift:ny-nshift-1L])+nshift
  imin = where(smtrace[nshift:ny-nshift-1L] le (shift(smtrace,nshift))[nshift:ny-nshift-1L] and smtrace[nshift:ny-nshift-1L] le (shift(smtrace,-nshift))[nshift:ny-nshift-1L])+nshift
  void = max(smtrace[imax],wmax)
  void = min(smtrace[imin],wmin)
  wmax = imax[wmax]
  wmin = imin[wmin]
  
  ;On estime la largeur a mi-hauteur des pics principaux
  lowsig_pos = where(smtrace lt smtrace[wmax]/2.)
  delpos = min(abs(lowsig_pos - wmax))
  lowsig_neg = where(-1.*smtrace lt -1.*smtrace[wmin]/2.)
  delneg = min(abs(lowsig_neg-wmin))
  
  ;On estime le niveau zero
  zerolevel = median(smtrace)
  
  ;Si on a une étoile binaire, on essaie de trouver les extrema secondaires
  xind = findgen(ny)
  if keyword_set(binary) then begin
    
    ;On fait fitter 4 fonctions de Moffat d'un coup dans la trace.
    sig_init = 3.
    wleft = (wmax < wmin)
    wright = (wmax > wmin)
    nan = !values.f_nan
    p0 = [zerolevel,wleft,nan,wright-wleft,1.,$
      trace[wleft],sig_init,nan,sig_init,trace[wright]/trace[wleft],sig_init,nan,sig_init]
    
    ;Parameter inputs
    parinfo = replicate({VALUE:!values.f_nan,FIXED:0B,LIMITED:[0B,0B],LIMITS:[0.,0.]},n_elements(p0))
    parinfo.VALUE = p0
    ;On limite les valeurs possibles des positions des centres de pics principaux a ± 1-sigma
    parinfo[[1L,3L]].LIMITED = bytarr(2,2)+1B
    ;On limite les positions des primaires a ± 1-sigma
    parinfo[[1L,3L]].LIMITS[1L] = p0[[1L,3L]] + sig_init
    ;Never reverse the left / right traces
    if parinfo[3L].limits[1L] lt 0. then parinfo[3L].limits[1L] = 0.
    ;On limite les largeurs entre 1.5 et 5 pixels
    parinfo[[6L,10L]].LIMITED = bytarr(2L,2L)+1B
    parinfo[[6L,10L]].LIMITS[0] = 1.5
    parinfo[[6L,10L]].LIMITS[1] = 15.5
    ;On limite les valeurs possibles pour l'indice de Moffat
    parinfo[4L].LIMITED = bytarr(2L)+1B
    parinfo[4L].LIMITS[0] = 0.51
    parinfo[4L].LIMITS[1] = 6.5
    
    nsep = 100
    seps_try = (findgen(nsep/2L)+1)/float(nsep/2L)*p0[3]/2.
    seps_try = [-reverse(seps_try),seps_try]
    prim = AB_moffat_fixed_sep(xind,p0*[1,1,1,1,1,1,1,0.,1,1,1,0.,1])
    chis = dblarr(nsep)+!values.d_nan
    pars = dblarr(n_elements(p0),nsep)+!values.d_nan
    g = where(finite(trace),ng)
    yfits = dblarr(ng,nsep)+!values.d_nan
    for isep=0L, nsep-1L do begin & $
      ;Start from a position on the grid 
      parinfo[2].value = seps_try[isep] & $
      parinfo[2].LIMITED = bytarr(2L) & $
      ;Estimate amplitude
      wleftseci = long(wleft+seps_try[isep]) & $
      wrightseci = long(wright+seps_try[isep]) & $
      amp_ratio_i = ((trace[wleftseci] - prim[wleftseci])/p0[5] > (trace[wrightseci] - prim[wrightseci])/(p0[5]*p0[9])) > .1 & $
      parinfo[[7L,11L]].LIMITED = bytarr(2L,2L) & $
      parinfo[7L].VALUE = amp_ratio_i & $
      parinfo[11L].VALUE = amp_ratio_i & $
      pari = mpfitfun('AB_moffat_fixed_sep',xind[g],trace[g],0., WEIGHTS=1D, YFIT=yfiti, /QUIET, PARINFO=parinfo,STATUS=status,ERRMSG=errmsg) & $
      ;plot, xind[g], trace[g] & $
      ;oplot, xind[g], yfiti, color=255 & $
      ;oplot, [wleftseci,wrightseci], interpol2(trace[g],xind[g],[wleftseci,wrightseci]), color=rvb_hex(100,255,100), psym=4, symsize=3, thick=2 & $
      ;oplot, xind[g], trace[g]-yfiti, color=rvb_hex(100,255,100) & $
      ;wait, .2 & $
      chis[isep] = total((trace[g]-yfiti)^2,/nan) & $
      pars[*,isep] = pari & $
      yfits[*,isep] = yfiti & $
    endfor
    ;set_plot,'X'
    ;plot, seps_try, 1./chis, /ps
    ;stop
    void = max(1./chis,/nan,wbest)
    par = pars[*,wbest]
    yfit = yfits[*,wbest]
    
    ;If the companion separation is negative (left) then change definition so that
    ; the A component is on the left. Unless the distinction is very clear
;Ya un probleme icitte, le code d'extraction plante
;    if par[2] lt 0. and par[7] gt .3 then begin; and (par[7] lt .5 or par[7] gt 2.)
;      par0 = par
;      par[1L] = par0[1L] + par0[2L]
;      par[2L] = -1 * par0[2L]
;      par[5L] = par0[5L] * par0[7L]
;      par[6L] = par0[8L]
;      par[7L] = 1. / par0[7L]
;      par[8L] = par0[6L]
;      par[9L] = par0[11] / par0[7] * par0[9]
;      par[10L] = par0[12L]
;      par[11L] = 1. / par0[11L]
;      par[12L] = par0[10L]
;    endif
    
  endif else begin
    
    ;Indices correspondant aux positions des pics
    subind = [2,6]
    
    ;On fait fitter les 2 fonctions de Moffat d'un coup dans la trace
    if ~keyword_set(only_amplitude) then begin
      p0 = [zerolevel,trace[wmax],wmax,delpos,1.,trace[wmin],wmin,delneg,1.]
      parinfo = replicate({VALUE:!values.f_nan,FIXED:0B,LIMITED:[0B,0B],LIMITS:[0.,0.]},n_elements(p0))
      parinfo.VALUE = p0
    endif else begin
      p0 = fit_params
      parinfo = replicate({VALUE:!values.f_nan,FIXED:0B,LIMITED:[0B,0B],LIMITS:[0.,0.]},n_elements(p0))
      parinfo.VALUE = p0
      parinfo.FIXED = 1B
      parinfo[subind-1L].FIXED = 0B
    endelse
    
    ;On empeche la largeur de devenir negative
    parinfo[subind+1L].LIMITED = [1B,1B]
    parinfo[subind+1L].LIMITS = [0.001,20.]
    parinfo[subind+2L].LIMITED = [1B,1B]
    parinfo[subind+2L].LIMITS = [0.51,6.]
    
    g = where(finite(trace),ng)
    par = mpfitfun('nmoffat_gmos',(findgen(ny))[g],trace[g],0., WEIGHTS=1D, YFIT=yfit, /QUIET, PARINFO=parinfo,STATUS=status,ERRMSG=errmsg)
    if status le 0 then message, 'MPFITFUN has failed ! Error message :'+errmsg
    
    ;On ordonne les traces de gauche a droite
    subind = subind[sort(par[subind])]
    
  endelse
  
  ;On affiche les resultats
  if keyword_set(display) then begin
    
    if display eq 2L then begin
      if !d.name ne 'Z' then set_plot,'Z'
      ;IND     0  1    2    3    4    5   6
      RED =   [0, 255, 255, 100, 100, 255,100]
      GREEN = [0, 255, 0,   255, 100, 100,155]
      BLUE =  [0, 255, 0,   100, 255, 100,100]
      ;Black,White,Red,Green,Blue,Pink
      TVLCT, RED, GREEN, BLUE
      cblack = 0
      cwhite = 1
      cred = 2
      cgreen = 3
      cblue = 4
      cpink = 5
      cdgreen = 6
      cviolet = 4
    endif else begin
      cblack = rvb_hex(0,0,0)
      cwhite = rvb_hex(255,255,255)
      cred = rvb_hex(255,0,0)
      cgreen = rvb_hex(100,255,100)
      cblue = rvb_hex(100,100,255)
      cviolet = rvb_hex(230,120,255)
      cpink = rvb_hex(255,100,100)
      cdgreen = rvb_hex(100,155,100)
    endelse
    
    if ~keyword_set(only_amplitude) then begin 
      xrangevec = [wmax+[-1,1]*delpos*10.,wmin+[-1,1]*delneg*10.]
      if ~keyword_set(xrange) then $
        xrange = [min(xrangevec,/nan),max(xrangevec,/nan)]
      if isa(display,'string') then title = display
    endif
    wset, 0
    plot, trace, TITLE=title, yrange=yrange, xrange=xrange
    if ~keyword_set(only_amplitude) then begin
      oplot, imin, trace[imin], color=cblue, /ps
      oplot, imax, trace[imax], color=cred, /ps
    endif
    oplot, xind, yfit, color=cgreen
    if ~keyword_set(only_amplitude) then begin
      if keyword_set(binary) then begin
        oplot, [par[1]],[trace[long(par[1])]], color=cred, psym=4,symsize=2, thick=2
        oplot, [par[1]+par[3]],[trace[long(par[1]+par[3])]], color=cred, psym=4,symsize=2, thick=2
        oplot, [par[1]+par[2]],[trace[long(par[1]+par[2])]], color=cviolet, psym=4,symsize=2, thick=2
        oplot, [par[1]+par[2]+par[3]],[trace[long(par[1]+par[2]+par[3])]], color=cviolet, psym=4,symsize=2, thick=2
      endif else begin
        oplot, [wmin], [trace[wmin]], color=cred, psym=4,symsize=3, thick=3
        oplot, [wmax], [trace[wmax]], color=cred, psym=4,symsize=3, thick=3
      endelse
    endif
    
    if display eq 2L then begin
      im = tvread(true=3)
      TVLCT, r, g, b, /Get
      im3d = BytArr(3, (size(im))[1], (size(im))[2])
      im3d[0,*,*] = r[im]
      im3d[1,*,*] = g[im]
      im3d[2,*,*] = b[im]
      date = curcompdate(/detail)
      string_replace, date, '@', '|'
      if ~keyword_set(title) then title = 'detect_traces_AB'
      dir = 'extract_traces_AB_logs/'
      file_mkdir,dir
      write_jpeg, dir+title+'_'+date+'.jpg', im3d, /true
    endif
    
  endif
  
  ;On choisit une region qui contient seulement du ciel
  ;Region "temoin" de 25 points contenant seulement du ciel
  stats = abs(trace/robust_sigma(trace,/nan))
  wsky = where(stats lt weighted_median(stats,medval=.3), nwsky)
  ;Max. 25 points de ciel
  nwsky <= 25L
  
  if ~keyword_set(binary) then begin
    trace_pos = par[subind]
  endif else begin
    trace_pos = [par[1L],par[1L]+par[2L],par[1L]+par[3L],par[1L]+par[2L]+par[3L]]
  endelse
  
  if nwsky ne 0 then begin
  
    sky_ind = list(LENGTH=(2L+2L*long(keyword_set(binary))))
    
    wsky = wsky[sort(abs(trace_pos[0L]-wsky))]
    sky_ind[0L] = wsky[0L:(nwsky-1L)]
    wsky = wsky[sort(abs(trace_pos[1L]-wsky))]
    sky_ind[1L] = wsky[0L:(nwsky-1L)]
    
    if keyword_set(binary) then begin
      wsky = wsky[sort(abs(trace_pos[2L]-wsky))]
      sky_ind[2L] = wsky[0L:(nwsky-1L)]
      wsky = wsky[sort(abs(trace_pos[3L]-wsky))]
      sky_ind[3L] = wsky[0L:(nwsky-1L)]
    endif
    
    if keyword_set(display) then $
      if display eq 1 then $
        for i=0L, n_elements(sky_ind)-1L do $
          oplot, xind[sky_ind[i]], trace[sky_ind[i]], color=rvb_hex(100,100,255), /ps
  endif
  
  if keyword_set(display) then begin
    if display eq 3 then begin
      date = curcompdate(/detail)
      string_replace, date, '@', '|'
      if ~keyword_set(title) then title = 'n'
      dir = 'extract_traces_AB_logs/'
      if ~file_test(dir) then file_mkdir,dir
      saveimage, dir+title+'_trace_detection_'+date+'.png', /PNG
      saveimage, dir+title+'_trace_detection.png', /PNG
    endif
  endif
  if keyword_set(stop) then stop
  
  ;On renvoie les positions des pics
  fit_params = par
  if ~keyword_set(binary) then begin
    peak_values = par[subind-1L]
    peak_widths = par[subind+1L]
    peak_indices = par[subind+2L]
  endif
  return, trace_pos
  
End