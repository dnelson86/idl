pro get_siminfo, Base, SnapBase, Num, NumGroupFiles, NumSnapFiles, Ngroups, Nsubgroups, NumPartType, Ti, BoxSize

  if num ge 1000 then begin
      exts='0000'
      exts=exts+strcompress(string(Num),/remove_all)
      exts=strmid(exts,strlen(exts)-4,4)
   endif else begin
      exts='000'
      exts=exts+strcompress(string(Num),/remove_all)
      exts=strmid(exts,strlen(exts)-3,3)
   endelse


  fnr = 0

  f = Base + "/groups_" + exts +"/fof_subhalo_tab_"+exts +"."+strcompress(string(fnr),/remove_all) + ".hdf5"
  print, f

  openr, 1, f, error = err
  if err ne 0 then begin
     f = Base + "/fof_subhalo_tab_"+exts + ".hdf5"
  endif else begin
     close, 1
  endelse
  
  
  print, "reading file:", f
   
  file_id = h5f_open(f)
  
  grp_id = h5g_open(file_id, 'Header')
  
  NumGroupFiles = read_attribute(grp_id, 'NumFiles')
  
  Ngroups  = read_attribute(grp_id, 'Ngroups_Total')
  
  Nsubgroups  = read_attribute(grp_id, 'Nsubgroups_Total')
  
  h5g_close, grp_id
  h5f_close, file_id
  

  
  f = Base + "/snapdir_" + exts +"/"+snapbase+"_"+exts +"."+strcompress(string(fnr),/remove_all) + ".hdf5"
  print, f
  
         
  openr, 1, f, error = err
  if err ne 0 then begin
     f = Base +"/"+snapbase+"_"+exts + ".hdf5"
  endif else begin
     close, 1
  endelse
  

  print, "reading file:", f
  
  file_id = h5f_open(f)
  
  grp_id = h5g_open(file_id, 'Header')

  NumSnapFiles = read_attribute(grp_id, 'NumFilesPerSnapshot')
  NumPartType  = read_attribute(grp_id, 'NumPart_Total')

  NumPart_Total_HighWord  = read_attribute(grp_id, 'NumPart_Total_HighWord')
  
  NumPartType = long64(NumPartType) + long64(NumPart_Total_HighWord)*(long64(2)^32)


  Ti = read_attribute(grp_id, 'Time')

  BoxSize = read_attribute(grp_id, 'BoxSize')


  h5g_close, grp_id
  h5f_close, file_id



end

