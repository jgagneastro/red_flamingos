Function read_abba_pattern, pattern_in, repeated_or_not
  
  pattern_vec = strtrim(transpose(byte(strlowcase(pattern_in))),2L)
  
  ;First, locate all the "a" positions
  pos_a = where(pattern_vec eq 'a', npos_a, complement=pos_b, ncomplement=npos_b)
  if npos_b eq 0 or npos_a eq 0 then message, ' There must be at least one A and one B positions in the ABBA patterns !'
  ;We'll also keep trace of the a or b positions that were already used.
  used_b = intarr(npos_b)
  ncouples = (npos_a > npos_b)
  couples = lonarr(2L,ncouples)-1L
  repeated_or_not = intarr(2L,ncouples)-1
  for i=0L, npos_a-1L do begin
    ;Find the nearest "b" position that was not already used
    gg = where(used_b eq 0, ngg)
    if ngg ne 0 then begin
      void = min(abs(pos_b[gg]-pos_a[i]),wmin)
      chosen_b = gg[wmin]
      couples[*,i] = [pos_a[i],pos_b[chosen_b]]
      repeated_or_not[*,i] = [0,0]
      used_b[chosen_b] = 1
    endif else begin
      ;If all b positions were already used, then just choose the nearest one.
      void = min(abs(pos_b-pos_a[i]),chosen_b)
      couples[*,i] = [pos_a[i],pos_b[chosen_b]]
      repeated_or_not[*,i] = [0,1]
    endelse
  endfor
  
  if n_elements(couples) eq 2L then begin
    couples = reform(couples,2L,1L)
    repeated_or_not = reform(repeated_or_not,2L,1L)
  endif
  
  ;If all b positions were already used here, then return
  if min(used_b) eq 1 then return, couples
  
  ;If not all b positions were already used, then use them now
  gg = where(~used_b, ngg)
  for i=0L, ngg-1L do begin
    void = min(abs(pos_a-pos_b[gg[i]]),chosen_a)
    couples[*,npos_a+i] = [pos_a[chosen_a],pos_b[gg[i]]]
    repeated_or_not[*,npos_a+i] = [1,0]
  endfor
  
  return, couples
  
End