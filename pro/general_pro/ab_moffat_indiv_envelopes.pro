Function AB_moffat_indiv_envelopes, xind, fit_params_in, NORM=norm
  
  forward_function AB_moffat_fixed_sep
  
  fit_params = fit_params_in
  fit_params_left = fit_params*[0.,1.,1.,1.,1.,1.,1.,0.,1.,0.,1.,0.,1.]
  primary_envelope_left = AB_moffat_fixed_sep(xind,fit_params_left)
  fit_params_right = fit_params_left
  fit_params_right[1L] = fit_params[1L] + fit_params[3L]
  fit_params_right[5L] = fit_params[5L] * fit_params[9L]
  fit_params_right[6L] = fit_params[10L]
  primary_envelope_right = AB_moffat_fixed_sep(xind,fit_params_right)
  sfit_params_left = fit_params_left
  sfit_params_left[1L] = fit_params[1L] + fit_params[2L]
  sfit_params_left[5L] = fit_params[5L] * fit_params[7L]
  sfit_params_left[6L] = fit_params[8L]
  sec_envelope_left = AB_moffat_fixed_sep(xind,sfit_params_left)
  sfit_params_right = fit_params_right
  sfit_params_right[1L] = fit_params[1L] + fit_params[2L] + fit_params[3L]
  sfit_params_right[5L] = fit_params[5L] * fit_params[9L] * fit_params[11L]
  sfit_params_right[6L] = fit_params[12L]
  sec_envelope_right = AB_moffat_fixed_sep(xind,sfit_params_right)
  
  envelopes = [[primary_envelope_left],[primary_envelope_right],[sec_envelope_left],[sec_envelope_right]]
  
  if keyword_set(norm) then $
    envelopes /= (make_array(n_elements(xind),value=1d,/double)#total(envelopes,1,/double,/nan))
  
  return, envelopes
End