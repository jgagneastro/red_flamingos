Function path_library, marker, HOSTNAME=hostname, SOFT=soft
;  This function takes a marker in entry, then returns the path corresponding to this marker, 
;    depending on which machine IDL is now running.
;  To get the exact name of the current machine, just type "spawn, 'hostname', hostname & print, hostname"
;    in the IDL console. Then you can start adding markers with their paths.
  
  forward_function gpathf
  compile_opt hidden
  on_error, 2
   
  ;Check on which machine this code is running
  if ~keyword_set(hostname) then $
    spawn, 'hostname', hostname
  hostname = strlowcase(hostname)
  if hostname eq 'antares' then hostname = 'arcturus'
  
  ;Find the path corresponding to the marker
  case strlowcase(marker) of
    'redflamingos_basedir' : $
      begin
      case hostname of
        strlowcase('antares') : p = ['/Users/gagne/Documents/IDL/IDL_library/Public/red_flamingos/']
        else : goto, badhost
      endcase
    end
    'idl_sav' : $
      begin
      case hostname of
        strlowcase('antares') : p = ['/Users/gagne/Documents/IDL/IDL_sav/']
        else : goto, badhost
      endcase
    end
    'idl_csv': $
      begin
      case hostname of
        strlowcase('antares') : p = ['/Users/gagne/Documents/IDL/IDL_csv/']
        else : goto, badhost
      endcase
    end
    'reduced_data': $
      begin
      case hostname of
        strlowcase('antares') : p = ['/Users/gagne/Dropbox/Data_Repository/REDUCED/']
        else : goto, badhost
      endcase
    end
    else : if keyword_set(soft) then return, '' else message, 'The marker "'+strtrim(marker,2)+'" could not be recognized.'
  endcase
  
  if n_elements(p) eq 1 then return, p[0]
  return, p
  
  BADHOST :
  if keyword_set(soft) then return, ''
  message, 'The machine "'+strtrim(hostname,2)+'" is not associated with a "'+strtrim(marker,2)+'" marker.'
End