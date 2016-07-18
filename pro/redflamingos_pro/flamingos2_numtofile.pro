Function flamingos2_numtofile, date, file_number, NOGZ=nogz
  forward_function addzero
  name = 'S20'+date+'S'+addzero(file_number,4L)+'.fits'
  if ~keyword_set(nogz) then name += '.gz'
  return, name
End