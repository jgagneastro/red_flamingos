Function simbad_data, name, US=us, DEBUG=debug, VERBOSE=verbose, SILENT=silent, COORD=coord, MAXDIST=maxdist, HTTP10=http10
;+
; NAME:
;       SIMBAD_DATA
;       
; PURPOSE:
;       RETRIEVE SIMBAD DATA FOR AN ARRAY OF OBJECT NAMES.
;       
; CALLING SEQUENCE:
;       structure = SIMBAD_DATA( name[, /US, /DEBUG, /VERBOSE, /SILENT] )
;
; INPUTS:
;       NAME = String array containing SIMBAD-resolvable names.
;       COORD = [nx2]-elements array containing object coordinates [not done yet].
;
; OPTIONAL INPUT KEYWORD:
;       /US      - Use cfa.harvard.edu server instead of u-strasbg.fr.
;       /DEBUG   - Stop when an error occurs.
;       /VERBOSE - Shows iteration number.
;       /SILENT  - Prevents WEBGET from displaying messages.
;       MAXDIST  - Maximum distance in arcseconds when giving COORDS. Default = 1"
;
; OUTPUTS:
;       STRUCTURE - An array of structures containing SIMBAD data for each object name provided.
;
; EXAMPLE:
;       Get data for 'HIP 6' and 'Barnard Star'.
;
;       IDL> structure = SIMBAD_DATA(['HIP 6','Barnard Star'])
;       IDL> help, structure, /structure
; 
; RESTRICTIONS:
;       (1) Names provided must be SIMBAD-resolvable.
;       (2) Coordinates search has not been implemented yet.
;       (3) Not all SIMBAD data is currently retrieved. However, adding "variables" for reading is easily done within the header of this code (see *ADD VARIABLES* below).
;       (4) Some symbols in object names might not properly convert to HTML. You can easily add some cases in the "HTML conversion" section just below the "ADD VARIABLES" section.
;
; PROCEDURES USED:
;       REPSTR(), ADD_TAGS, WEBGET(), VALID_NUM(), REMOVE, TEN()
;
; MODIFICATION HISTORY:
;       WRITTEN, Jonathan Gagne, August, 19 2012
;       Fixed cases where SIMBAD would give only part of the coordinates, like for NAME ORI LOOP
;-
  
  ;Subroutines
  forward_function repstr, add_tags, webget, valid_num, remove, ten
  
  ; **** ADD VARIABLES HERE ********************************************************************************************************
  
  ;2D measurements
  vars2D = ['ICRS','FK5','FK4','GAL','pm'];,'size'
  tags2D = ['ICRS,ep=J2000,eq=2000','FK5,ep=J2000,eq=2000','FK4,ep=B1950,eq=1950',$
            'Gal,ep=J2000,eq=2000','Proper motions:'];,'Angular size:'
  ;Whether or not to use ten for converting sexagecimal to decimal data. Sexagecimal data will be stored as a double-precision scalar.
  ;It is also assument that sexadecimal data will be in format (hh,mm,ss), (dd,mm,ss) like it's the case for (RA,DEC)
  useten = [1B,1B,1B,0B,0B];,0B
  ;Whether or not there is an Opt / IR tag for this data.
  ir_tag = [1B,1B,1B,1B,0B];,0B
  
  ;1D measurements
  vars1D = ['plx','vrad','reds','cz','UMAG','BMAG','VMAG','GMAG','RMAG','IMAG','ZMAG','YMAG','JMAG','HMAG','KMAG']
  tags1D = ['Parallax:','Radial Velocity:','Redshift:','cz:','Flux U :','Flux B :','Flux V :','Flux G :','Flux R :','Flux I :','Flux Z :','Flux Y :','Flux J :','Flux H :','Flux K :']
  
  ;1D string measurements
  vars1DS = ['spt']
  tags1DS = ['Spectral type:']
  
  ;Stack string measurements
  varsSS = ['ids','bibcodes']
  tagsSS = ['Identifiers (','Bibcodes ']
  
  ; ********************************************************************************************************************************
  
  ;HTML conversion
  if keyword_set(name) then begin
    object = repstr(name,'+','%2B')
    object = repstr(object,' ','%20')
    object = repstr(object,'[','%5B')
    object = repstr(object,']','%5D')
  endif
  
  ;Server
  ;us = 1
  web = 'u-strasbg.fr'
  if keyword_set(us) then web = 'cfa.harvard.edu'
  if keyword_set(coord) and ~keyword_set(maxdist) then maxdist = 1.
  
  if (size(coord))[0] ne 1L then $
    if (size(coord))[2] gt (size(coord))[1] and (size(coord))[2] ne 2L then $
      message, 'You have probably inverted the dimensions for COORD ! They must be N x 2 !' 
  
  if keyword_set(coord) then begin
    signs = strarr((size(coord))[1])+'-'
    bad = where(coord[*,1L] ge 0., nbad)
    if nbad ne 0 then signs[bad] = '++'
    URL = 'http://simbad.'+web+'/simbad/sim-coo?&Coord=%20'+strtrim(coord[*,0L],2)+'%20'+signs+strtrim(abs(coord[*,1L]),2)+'&output.format=ASCII'
  endif else $
    URL = 'http://simbad.'+web+'/simbad/sim-id?Ident='+object+'&output.format=ASCII'
  
  ;Variable names
  post_1DS = ['','_qual','_ref']
  pre_1D = ['','e','','']
  post_1D = ['','','_qual','_ref']
  pre_2D = ['','','e','e','','','']
  post_2D = ['_ra','_dec','_ra','_dec','_corr','_qual','_ref']
  
  ;Build a list of all variables
  nvars2D = n_elements(vars2D)
  nvars1D = n_elements(vars1D)
  nvars1DS = n_elements(vars1DS)
  nvarsSS = n_elements(varsSS)
  n_2D = n_elements(pre_2D)
  n_1D = n_elements(pre_1D)
  n_1DS = n_elements(post_1DS)
  vars = ['id','type']
  types = [7L,7L]
  vars = [vars, (transpose(pre_2D[transpose(lindgen(n_2D,nvars2D) mod n_2D)] + vars2D[lindgen(nvars2D,n_2D) mod nvars2D] + post_2D[transpose(lindgen(n_2D,nvars2D) mod n_2D)]))[*]]
  types = [types, (transpose(([4L,5L])[useten[lindgen(nvars2D,n_2D) mod nvars2D]]))[*]]
  vars = [vars, (transpose(pre_1D[transpose(lindgen(n_1D,nvars1D) mod n_1D)] + vars1D[lindgen(nvars1D,n_1D) mod nvars1D] + post_1D[transpose(lindgen(n_1D,nvars1D) mod n_1D)]))[*]]
  types = [types, lonarr(n_1D*nvars1D)+4L]
  vars = [vars, (transpose(vars1DS[lindgen(nvars1DS,n_1DS) mod nvars1DS] + post_1DS[transpose(lindgen(n_1DS,nvars1DS) mod n_1DS)]))[*]]
  types = [types, lonarr(n_1DS*nvars1DS)+7L]
  vars = [vars, varsSS, 'notes']
  types = [types, lonarr(nvarsSS+1L)+7L]
  nvars = n_elements(vars)
  bad = where(strpos(vars,'_qual') ne -1 or strpos(vars,'_ref') ne -1, nbad)
  if nbad ne 0 then types[bad] = 7L
  
  ;Create the final structure that will contain everything
  nans = ['','0B','-1','-1L','!values.f_nan','!values.d_nan','complex(1,1)*!values.f_nan','''NaN''',$
          '','dcomplex(1,1)*!values.d_nan','','','0U','0UL','-1LL','0ULL']
  add_tags, {ID:'NaN'}, vars[1:*], nans[types[1:*]], struct
  if keyword_set(name) then $
    ns = n_elements(name)
  if keyword_set(coord) then $
    ns = (size(coord))[1]
  struct = replicate(struct,ns) 
  
  ;Loop over possible array
  for i=0L, ns-1L do begin
    if keyword_set(name) then $
      text_name = name[i] else $
      text_name = '"'+strjoin(strtrim(reform(coord[i,*]),2),',')+'"'
    if keyword_set(verbose) then print, ' Parsing ['+strtrim(i+1,2)+'/'+strtrim(ns,2)+'] : '+text_name+' ...'
    
    s = webget(URL[i], SILENT=silent, HTTP10=http10)
    ;Avoid using webget, which sucks
    ;ourl = obj_new('IDLnetURL')
    ;text = ourl->get(/string_array,url=URL[i])
    ;s = {text:text}
    if ~isa(s,'struct') then begin
      if ~keyword_set(silent) then message, ' Object '+strtrim(i,2)+' : No data downloaded !', CONTINUE=~keyword_set(debug)
      continue
    endif
    if n_elements(s.text) le 4 then begin
      if ~keyword_set(silent) then message, ' Object '+strtrim(i,2)+' : No data downloaded !', CONTINUE=~keyword_set(debug)
      continue
    endif
    
    if keyword_set(coord) then begin
      gline = (where(strpos(s.text,'1|') ne -1L,ngl))[0]
      if ngl eq 0 then begin
        message, ' Object '+strtrim(i,2)+' : No data downloaded !', CONTINUE=~keyword_set(debug)
        continue
      endif
      data = strsplit(s.text[gline],'|',/extract)
      if float(data[1]) gt maxdist then begin
        message, ' Object '+strtrim(i,2)+' : No data within radius !', CONTINUE=~keyword_set(debug)
        continue
      endif
      name_search = strtrim(data[2],2)
      structi = simbad_data(name_search, US=us, DEBUG=debug, VERBOSE=verbose, SILENT=silent)
      if ~isa(structi,'struct') then continue
      struct[i] = structi
      continue
    endif
    
    ;Remove first 4 lines
    text = s.text[4:*]
    
    ;Object
    id = 'NaN'
    type = 'NaN'
    gobj = (where(strpos(text,'Object') ne -1, ngobj))[0]
    if ngobj ne 0 then begin
      id = strtrim(strmid(text[gobj],strpos(text[gobj],'Object ')+7L,strpos(text[gobj],' --- ')-(strpos(text[gobj],'Object ')+7L)),2)
      void = strtrim(strmid(text[gobj],strpos(text[gobj],' --- ')+5L),2)
      type = strtrim(strmid(void,0,strpos(void,' --- ')),2)
    endif
    
    for j=0L, nvars2D-1L do begin
      
      ;The order is messed up in SIMBAD for data with (Opt) or (IR) tag -- Only when using ASCII output (this is strange).
      ; See example ; http://simbad.u-strasbg.fr/simbad/sim-id?Ident=HIP%206&output.format=ASCII
      ; The order issue might get fixed sooner and later, then we'll need to update this code.
      if ir_tag[j] then begin
        pre_2Dj = [pre_2D[0:1],'',pre_2D[5],pre_2D[2:4],pre_2D[6]]
        post_2Dj = [post_2D[0:1],'_type',post_2D[5],post_2D[2:4],post_2D[6]]
        ;Initialize NaN values
        void = execute(strjoin(pre_2Dj+vars2D[j]+post_2Dj+' = '+[strarr(2)+'!values.'+(['d','f'])[useten[j]]+'_nan',strarr(2)+'''NaN''',strarr(3)+'!values.'+(['d','f'])[useten[j]]+'_nan','''NaN'''],' & '))
      endif else begin
        pre_2Dj = pre_2D
        post_2Dj = post_2D
        ;Initialize NaN values
        void = execute(strjoin(pre_2Dj+vars2D[j]+post_2Dj+' = '+[strarr(5)+'!values.'+(['d','f'])[useten[j]]+'_nan',strarr(2)+'''NaN'''],' & '))
      endelse
      
      void *= execute(strjoin(['a','b','p']+' = '+[strarr(3)+'!values.'+(['d','f'])[useten[j]]+'_nan'],' & '))
      if void eq 0 then message, ' A code execution has failed !'
      
      ;Search for data
      gi = (where(strpos(text,tags2D[j]) ne -1, ngi))[0]
      if ngi eq 0 then continue
      
      ;Cut the text string
      text[gi] = repstr(text[gi],' )','')
      text[gi] = repstr(text[gi],')','')
      text[gi] = repstr(text[gi],'(','')
;      tt = text[gi]
;      strkill, tt, 'eq=2000:'
;      text[gi] = tt
      tsplit = strtrim(strsplit(strmid(text[gi],strpos(text[gi],tags2D[j])+strlen(tags2D[j])),' ',/extract),2)
      tsplit = repstr(tsplit,'[','')
      tsplit = repstr(tsplit,']','')
      if strtrim(tsplit[0],2) eq ':' then $
        remove, 0, tsplit
      nsp = n_elements(tsplit)
      
      ;Read first part of data
      void = 1B
      if useten[j] then begin
        ;Two special cases where SIMBAD gives incomplete coordinates
        if nsp ge 8L then begin
          ncoord = nsp - 6L
          symbol_pos = (where(strpos(tsplit,'-') ne -1 or strpos(tsplit,'+')  ne -1, nf))[0]
          if nf ne 0 and symbol_pos le ncoord then begin
            ravec = fltarr(3)
            decvec = fltarr(3)
            for kk=0L, symbol_pos-1L do $
              ravec[kk] = tsplit[kk]
            for kk=symbol_pos, ncoord-1L do $
              decvec[kk-symbol_pos] = tsplit[kk]
            void = 1B
            if valid_num(tsplit[0]) and valid_num(tsplit[1]) then $
              void *= execute(pre_2Dj[0]+vars2d[j]+post_2Dj[0]+' = ten(ravec[0],ravec[1],ravec[2])*15d0')
            if valid_num(tsplit[2]) then $
              void *= execute(pre_2Dj[1]+vars2d[j]+post_2Dj[1]+' = ten(decvec[0],decvec[1],decvec[2])')
            if void eq 0B then message, ' A code execution has failed !'
          endif
          k0 = ncoord
        endif else begin
          if nsp ge 3L then $
            if valid_num(tsplit[0]) and valid_num(tsplit[1]) and valid_num(tsplit[2]) then $
              void *= execute(pre_2Dj[0]+vars2d[j]+post_2Dj[0]+' = ten(tsplit[0],tsplit[1],tsplit[2])*15d0')
          if nsp ge 6L then $
            if valid_num(tsplit[3]) and valid_num(tsplit[4]) and valid_num(tsplit[5]) then $
              void *= execute(pre_2Dj[1]+vars2d[j]+post_2Dj[1]+' = ten(tsplit[3],tsplit[4],tsplit[5])')
          if void eq 0B then message, ' A code execution has failed !'
          k0 = 6L
        endelse
      endif else begin
        if valid_num(tsplit[0]) then $
          void *= execute(pre_2Dj[0]+vars2d[j]+post_2Dj[0]+' = tsplit[0]')
        if valid_num(tsplit[1]) then $
          void *= execute(pre_2Dj[1]+vars2d[j]+post_2Dj[1]+' = tsplit[1]')
        k0 = 2L
      endelse
      
      ; Only the first coordinates line is complete in the ASCII ouput, compared with HTML.
      ; See for comparison ; http://simbad.u-strasbg.fr/simbad/sim-id?Ident=HIP%206
      ; and http://simbad.u-strasbg.fr/simbad/sim-id?Ident=HIP%206&output.format=ASCII
      ; That's why we skip the rest of the loop if the "problem" happens.
      if nsp le k0 then continue
      
      ;Read second part of data (OPT/IR, quality, reference)
      void = execute(vars2D[j]+'_type = tsplit['+strtrim((where(post_2Dj eq '_type'))[0]+k0-2L,2)+']')
      void *= execute(vars2D[j]+'_qual = tsplit['+strtrim((where(post_2Dj eq '_qual'))[0]+k0-2L,2)+']')
      void *= execute(vars2D[j]+'_ref = tsplit['+strtrim((where(post_2Dj eq '_ref'))[0]+k0-2L,2)+']')
      void *= execute(strjoin('if '+vars2D[j]+['_type','_qual','_ref']+' eq ''~'' then '+vars2D[j]+['_type','_qual','_ref']+' = ''NaN''',' & '))
      if void eq 0B then message, ' A code execution has failed !'
      
      ;Read the last part (error ellipse)
      ind = [where(post_2Dj eq '_ra' and pre_2Dj eq 'e'),where(post_2Dj eq '_dec' and pre_2Dj eq 'e'),where(post_2Dj eq '_corr' and pre_2Dj eq '')]+k0-2L
      void = execute(strjoin('if valid_num(tsplit['+strtrim(ind,2)+']) then '+['a','b','p']+' = '+(['float','double'])[useten[j]]+'(tsplit['+strtrim(ind,2)+'])',' & '))
      if void eq 0B then message, ' A code execution has failed !'
      
      ;If the error ellipse is provided, convert to the actual error and correlation coefficient
      if finite(a*b*p) then begin
        era = sqrt(sin(p*!dtor)^2*a^2 + cos(p*!dtor)^2*b^2)
        edec = sqrt(cos(p*!dtor)^2*a^2 + sin(p*!dtor)^2*b^2)
        corr = (a^2-b^2)*cos(p*!dtor)*sin(p*!dtor)/(era*edec)
        ind -= (k0-2L)
        void = execute(strjoin(pre_2Dj[ind]+vars2d[j]+post_2Dj[ind]+' = '+['era','edec','corr'],' & '))
        if void eq 0B then message, ' A code execution has failed !'
      endif
    endfor
    
    ;Loop over all variables to read the data
    for j=0L, nvars1D-1L do begin
      
      ;Initialize NaN values
      void = execute(strjoin(pre_1D+vars1D[j]+post_1D+' = '+[strarr(2)+'!values.f_nan', strarr(2)+'''NaN'''],' & '))
      if void eq 0B then message, ' A code execution has failed !'
      
      ;Search for data
      gi = (where(strpos(text,tags1D[j]) ne -1, ngi))[0]
      if ngi eq 0 then continue
      
      ;Cut the text string
      tsplit = strtrim(strsplit(strmid(text[gi],strpos(text[gi],tags1D[j])+strlen(tags1D[j])),' ',/extract),2)
      tsplit = repstr(tsplit,'[','')
      tsplit = repstr(tsplit,']','')
      if strtrim(tsplit[0],2) eq ':' then $
        remove, 0, tsplit
      
      ;Read data
      void = execute(strjoin(['if valid_num(tsplit['+['0','1']+']) then ',strarr(2)]+pre_1D+vars1D[j]+post_1D+' = '+[strarr(2)+'float(',strarr(2)]+'tsplit['+strtrim(indgen(4),2)+']'+[strarr(2)+')',strarr(2)],' & '))
      void *= execute(strjoin('if '+vars1D[j]+['_qual','_ref']+' eq ''~'' then '+vars1D[j]+['_qual','_ref']+' = ''NaN''',' & '))
      if void eq 0B then message, ' A code execution has failed !'
    endfor
    
    ;Loop over all variables to read the data
    for j=0L, nvars1DS-1L do begin
      
      ;Initialize NaN values
      void = execute(strjoin(vars1DS[j]+post_1DS+' = '+[strarr(3)+'''NaN'''],' & '))
      if void eq 0B then message, ' A code execution has failed !'
      
      ;Search for data
      gi = (where(strpos(text,tags1DS[j]) ne -1, ngi))[0]
      if ngi eq 0 then continue
      
      ;Cut the text string
      tsplit = strtrim(strsplit(strmid(text[gi],strpos(text[gi],tags1DS[j])+strlen(tags1DS[j])),' ',/extract),2)
      if strtrim(tsplit[0],2) eq ':' then $
        remove, 0, tsplit
      
      ;Read data
      void = execute(strjoin(vars1DS[j]+post_1DS+' = tsplit['+strtrim(indgen(3),2)+']',' & '))
      void *= execute(strjoin('if '+vars1DS[j]+post_1DS+' eq ''~'' then '+vars1DS[j]+post_1DS+' = ''NaN''',' & '))
      if void eq 0B then message, ' A code execution has failed !'
    endfor
    
    ;Loop over all variables to read the data
    for j=0L, nvarsSS-1L do begin
      ;Initialize NaN values
      void = execute(varsSS[j]+' = ''NaN''')
      if void eq 0B then message, ' A code execution has failed !'
      
      ;Search for data
      gi = (where(strpos(text,tagsSS[j]) ne -1, ngi))[0]
      if ngi eq 0 then continue
      
      ;Check number of entries
      nen = fix(repstr(repstr(repstr((strsplit(text[gi],'(',/extract))[-1],')',''),' ',''),':',''))
      gi += 1L
      
      ;Read data
      vec = []
      l = 0L
      for k=0L, 9999L do begin
        if gi+k ge n_elements(text)-1L then break
        data = strtrim(strsplit(repstr(text[gi+k],'  ','|'),'|',/extract),2)
        bad = where(strlen(data) le 3, nbad)
        if nbad eq n_elements(data) then continue
        if nbad ne 0 then remove, bad, data
        vec = [vec, data]
        l += n_elements(data)
        if l ge nen then break
      endfor
      vec = strjoin(vec,'|')
      
      void = execute(varsSS[j]+' = vec')
      if void eq 0B then message, ' A code execution has failed !'
    endfor
    
    ;Notes
    notes = 'NaN'
    gi = (where(strpos(text,'Notes (') ne -1, ngi))[0]
    
    if ngi ne 0 then begin
      
      ;Check number of entries
      nen = fix(repstr(repstr(repstr((strsplit(text[gi],'(',/extract))[-1],')',''),' ',''),':',''))
      if nen ne 0 then begin
        
        ;Read the notes
        gi += 1L
        vec = []
        l = 0L
        for k=0L, 9999L do begin
          if gi+k ge n_elements(text)-1L then break
          if strmid(strtrim(text[gi+k],2),0,1) eq '(' then continue
          vec = [vec, text[gi+k]]
          l += 1L
          if l ge nen then break
        endfor
        
        notes = strjoin(strtrim(vec,2),'|')
      
      endif
    endif
    
    void = 1B
    for j=0L, nvars-1L do begin
      void = execute('struct[i].'+vars[j]+' = '+vars[j])
      if void eq 0B then message, ' A code execution has failed !'
    endfor
    if void eq 0B then message, ' A code execution has failed !'
    
  endfor
  return, struct
  
End