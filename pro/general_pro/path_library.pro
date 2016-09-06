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
  
  ;Find the path corresponding to the marker
  case strlowcase(marker) of
    'redflamingos_basedir' : $
      begin
      case hostname of
        strlowcase('antares') : p = ['/Users/gagne/Documents/IDL/IDL_library/Public/red_flamingos/']
        else : goto, badhost
      endcase
    end
    'redflamingos_sav' : $
      begin
      case hostname of
        strlowcase('antares') : p = [gpath('redflamingos_basedir')+'idl_sav/']
        else : goto, badhost
      endcase
    end
    'redflamingos_logdir': $
      begin
      case hostname of
        strlowcase('antares') : p = [gpath('redflamingos_basedir')+'flamingos_logfiles/']
        else : goto, badhost
      endcase
    end
    'redflamingos_raw_data': $
      begin
      case hostname of
        strlowcase('antares') : p = [gpath('redflamingos_basedir')+'sample_data/raw/']
        else : goto, badhost
      endcase
    end
    'redflamingos_reduced_data': $
      begin
      case hostname of
        strlowcase('antares') : p = [gpath('redflamingos_basedir')+'sample_data/reduced/']
        else : goto, badhost
      endcase
    end
    'redflamingos_screenshots': $
      begin
      case hostname of
        strlowcase('antares') : p = [gpath('redflamingos_basedir')+'screenshots/']
        else : goto, badhost
      endcase
    end
    'redflamingos_darks': $
      begin
      case hostname of
        strlowcase('antares') : p = [gpath('redflamingos_raw_data')+'DARKS/']
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