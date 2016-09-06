function clipnsmooth, xvals, datav, varv, specv, eval, coeffv

; CHECK INPUTS
nx = n_elements(xvals)
coeffv = [mean(datav/specv)]

; EVALUATE
gv = where(specv ne 0,numgood, complement = bv) ; good and bad pixel locations
if (gv[0] eq -1 or n_elements(gv) le 3) then return, dblarr(nx) ; a row of zeros of right length
goodx = xvals[gv]

savgolfilt = SAVGOL(15, 15, 0, 4, /DOUBLE)
if n_elements(goodx) lt n_elements(savgolfilt) then $  
  savgolfilt = savgol(round(n_elements(goodx)/2)-1,round(n_elements(goodx)/2)-1, 0, 4, /double)

estg = CONVOL(datav[gv]/specv[gv], savgolFilt, /edge_truncate)

fiteval = dblarr(nx)
if (bv[0] ne -1) then begin      ; interpolate the bad from good pixels
  badx = xvals[bv]
  estb = interpol(estg, goodx, badx, /quadratic)
  fiteval[badx ] = estb
endif
fiteval[goodx] = estg

return, fiteval

end
