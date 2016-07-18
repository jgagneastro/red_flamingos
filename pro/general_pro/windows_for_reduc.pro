Pro windows_for_reduc, MORE=more, LONGTOP=longtop
  
  ;Use !version.os_family
  spawn,'system_profiler SPDisplaysDataType',output
  thunderbolt_connected = max(strpos(output,'Thunderbolt Display:') ne -1,/nan)
  
  factor = 1.
  if keyword_set(longtop) then factor = 2.8
  if thunderbolt_connected then begin
    window,2,xsize=900.,xpos=1285.,ysize=400,ypos=140
    if keyword_set(more) then begin
      window,3,xsize=300.,xpos=2190.,ysize=1280,ypos=140;Not adjusted
      window,4,xsize=300.,xpos=2495.,ysize=1280,ypos=140;Not adjusted
    endif
    window,1,xsize=900.*factor,xpos=1285.,ysize=400,ypos=575
    window,0,xsize=900.*factor,xpos=1285.,ysize=400,ypos=1100
  endif else begin
    if keyword_set(more) then begin
      window,2,xsize=450.,xpos=10.,ysize=200,ypos=70
      window,3,xsize=200.,xpos=465.,ysize=710,ypos=600;Not adjusted
      window,4,xsize=200.,xpos=670.,ysize=710,ypos=600;Not adjusted
      window,1,xsize=450.*factor,xpos=10.,ysize=200,ypos=275
      window,0,xsize=450.*factor,xpos=10.,ysize=300,ypos=500
    endif else begin
      window,2,xsize=600.,xpos=10.,ysize=200,ypos=70
      window,1,xsize=600.*factor,xpos=10.,ysize=200,ypos=275
      window,0,xsize=600.*factor,xpos=10.,ysize=300,ypos=500
    endelse
  endelse
End