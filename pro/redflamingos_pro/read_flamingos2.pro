Function read_flamingos2, file, useful_region, HEADER=header, REMOVEBRIGHTLINES=removebrightlines
  forward_function mrd_hread2
  im = double(readfits(file,/silent,ext=1L))
  header = headfits(file)
  if ~keyword_set(useful_region) then return, transpose(im)
  ;return, transpose(im[0:useful_region[0L],useful_region[1L]:*])
  return, transpose(im[useful_region[0L]:useful_region[1L],useful_region[2L]:useful_region[3L]])
End