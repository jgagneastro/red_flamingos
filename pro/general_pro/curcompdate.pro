Function curcompdate, KEEPNIGHT=keepnight, DETAIL=detail, FORMAT=format, PATH=path, UTC=utc, RANDOM=random
;+
; NAME:
;       CURCOMPDATE
;       
; PURPOSE:
;       Returns the current computer date in the format AAMMDD.
;       
; CALLING SEQUENCE:
;       string = CURCOMPDATE( [, FORMAT=format, /KEEPNIGHT, /DETAIL, /PATH] )
;
; INPUTS:
;       None.
;
; OPTIONAL INPUT:
;       FORMAT - If you choose format=2, the the output format will be AA-MM-DD. Format=1
;                is the default one.
;                    
; OPTIONAL INPUT KEYWORD:
;       /KEEPNIGHT - If this keyword is set then 0.5 days are removed from the computer date.
;                    This is useful during an night of observation for example. 
;       /DETAIL    - Adds the computer time in the format @hh:mm:ss
;       /PATH      - Replaces the : by _ if run on a Windows machine (: are not allowed in paths).
; 
; OUTPUTS:
;       STRING - A scalar string consisting of the computer date in the format AAMMDD. 
;
; RESTRICTIONS:
;       (1) Your system time as returned by systime() must be in english.
;
; PROCEDURES USED:
;       DATEWRITE, SYSJUL(), ADDZERO(), STRREPLACE()
;
; MODIFICATION HISTORY:
;       WRITTEN, Jonathan Gagne, September, 3 2011
;-
  
  compile_opt hidden ;Evite de toujours afficher 'Compiled module'
  forward_function datewrite, sysjul, addzero, string_replace, getseed
  on_error, 2
  
  if ~keyword_set(format) then format = 1
  if keyword_set(path) and n_elements(detail) eq 0 then detail = 1
  
  if keyword_set(keepnight) then begin
    datewrite, sysjul()-0.5d0, date, /julian
    return, date
  endif
  
  months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
  ss = systime(UTC=utc)
  mon = where(months eq strmid(ss,4,3))
  if mon[0] eq -1 then return, '-1'
  case format of
    1 : finstr = strmid(ss,22,2)+addzero(fix(strcompress(string(mon[0]+1),/remove_all)))+addzero(fix(strmid(ss,8,2))) + ( keyword_set(detail) ? '@'+strmid(ss,11,8) : '' )
    2 : finstr = strmid(ss,20,4)+'-'+addzero(fix(strcompress(string(mon[0]+1),/remove_all)))+'-'+addzero(fix(strmid(ss,8,2))) + ( keyword_set(detail) ? '@'+strmid(ss,11,8) : '' )
    else : return, 'ERROR'
  endcase
  
  if keyword_set(path) then string_replace, finstr, ':', '_';and !version.os_family eq 'Windows' then
  
  if keyword_set(random) then finstr += strtrim(long(getseed()*1d5),2)
  
  return, finstr
  
End
