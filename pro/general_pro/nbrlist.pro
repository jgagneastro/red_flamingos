Function nbrlist, tlist, RNS=rns  
;  Retourne un array de type integer contenant tous les nombres
;  specifies par une liste du format '1-5,10,11,18-20'
;  
;  tlist: string, specification de la liste de nombres
;                 "A-B": de A jusqu'a B inclusivement  
;                 ",": pour separer les series de nombre consecutifs
;                 ex. '1-10,12,15-18'
;  
;  keyword optionel:
;   rns=[r,s]: "Read'N Skip", pour repeter un patron dans une sequence
;              a partir du debut de la liste, inclut les r premiers
;              nombres et saute les s suivants, et recommence
;              jusqu'a la fin de la liste
  
  on_error, 2
  
  if tlist eq '' then return, !NULL
  
  txtlist=tlist
  nchar=strlen(txtlist)
  
  if (keyword_set(rns)) then period=rns[0]+rns[1]
  
  while nchar gt 0L do begin
    comma = strpos(txtlist,',')
    if comma eq -1L then comma=nchar+1L
    dash=strpos(txtlist,'-')
    if dash gt comma then dash=-1L
    ;if there is a dash in this part of the string
    if dash ne -1 then begin
        ni = long(strmid(txtlist,0L,dash))
        nf = long(strmid(txtlist,dash+1L,comma-dash-1))

        if (keyword_set(rns)) then begin
            for n=ni,nf do begin
              if ((n-ni) mod period lt rns[0L]) then $
                  if (n_elements(list) eq 0L) then list = n else list = [list, n]
            endfor
        endif else begin
            for n=ni,nf do $
                if (n_elements(list) eq 0L) then list = n else list = [list, n]
        endelse
    ;if no dash in this part of the string
    endif else begin
        n=long(strmid(txtlist,0L,comma))
        if (n_elements(list) eq 0L) then list = n else list = [list, n]
    endelse

    ;number of char not read yet
    nchar = nchar - comma - 1L
    nchar = (nchar > 0L)
    
    ;trim txtlist to remove characters read
    if nchar gt 0 then txtlist = strmid(txtlist, comma+1L, nchar)
  endwhile
  
  return,list
End
