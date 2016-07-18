Pro sortmultiple, order, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, $
              v13, v14, v15, v16, v17, v18, v19, v20, v21, v22, v23, $
              v24, v25, v26, v27, v28, v29, v30, v31, v32, v33, v34, $
              v35, v36, v37, v38, v39, v40, v41, v42, v43, v44, v45, $
              v46, v47, v48, v49, v50, v51, v52, v53, v54, v55, $
              v56, v57, v58, v59, v60, v61, v62, v63, v64, v65, $
              v66, v67, v68, v69, v70, v71, v72, v73, v74, v75, $
              v76, v77, v78, v79, v80, SORT=sort
;  Sert a classer plusieurs vecteurs selon un ordre donne
  compile_opt hidden
  on_error, 2
  nparams = n_params()-1
  for i=1, nparams do $
    R = execute('v'+strtrim(i,2)+' = v'+strtrim(i,2)+'[order]')
  if keyword_set(sort) then begin
    cmd = 'sortmultiple, sort(order), '
    cmdv = 'v'+strtrim(indgen(nparams)+1,2)
    cmd += strjoin(cmdv,', ')
    void = execute(cmd)
  endif
End