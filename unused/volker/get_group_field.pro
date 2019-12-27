function get_group_field, FieldName, Base, Num, NumGroupFiles, Ngroups

      
  if num ge 1000 then begin
     exts='0000'
     exts=exts+strcompress(string(Num),/remove_all)
     exts=strmid(exts,strlen(exts)-4,4)
  endif else begin
     exts='000'
     exts=exts+strcompress(string(Num),/remove_all)
     exts=strmid(exts,strlen(exts)-3,3)
  endelse
  
   skip = 0L
   skipids=0L
   skipsub=0L
   fnr = 0L

   repeat begin

      f = Base + "/groups_" + exts +"/fof_subhalo_tab_"+exts +"."+strcompress(string(fnr),/remove_all) + ".hdf5"

      if fnr eq 0 then begin
                                ; check whether file exists, otherwise
                                ; assume the data is in a single file                                                                        
         openr, 1, f, error = err
         if err ne 0 then begin
            f = Base + "/fof_subhalo_tab_"+exts + ".hdf5"
         endif else begin
            close, 1
         endelse
      endif

     ; print, "reading file:", f

      file_id = h5f_open(f)

      grp_id = h5g_open(file_id, 'Header')
      nloc = read_attribute(grp_id, 'Ngroups_ThisFile')


  
      if fnr eq 0 then begin
  
         dataset_id= h5d_open(file_id, 'Group/' + FieldName)

         dataspace_id = H5D_GET_SPACE(dataset_id)
         datatype_id = H5D_GET_TYPE(dataset_id)

         idltype = H5T_IDLTYPE(Datatype_id)
         dimensions = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)

         
         ;print, "IdlType=", idltype, "  Dims=", dimensions
         
         h5d_close, dataset_id

         if n_elements(dimensions) eq 2 then begin
            firstdim = dimensions(0)
            if (idltype eq 3) or (idltype eq 13) then begin
               data = lonarr(firstdim, Ngroups) 
            endif 
            if idltype eq 4 then begin
               data = fltarr(firstdim, Ngroups) 
            endif 
            if idltype eq 5 then begin
               data = dblarr(firstdim, Ngroups) 
            endif 
      endif else begin
            firstdim = 1
            if (idltype eq 3) or (idltype eq 13) then begin
               data = lonarr(Ngroups) 
            endif 
            if idltype eq 4 then begin
               data = fltarr(Ngroups) 
            endif 
            if idltype eq 5 then begin
               data = dblarr(Ngroups) 
            endif 
         endelse
      endif


      if nloc gt 0 then begin

         if firstdim gt 1 then begin
            data(*, skip:skip + nloc -1) = read_dataset(file_id, 'Group/' + FieldName)
         endif else begin
            data(skip:skip + nloc -1) = read_dataset(file_id, 'Group/' + FieldName)
         endelse
      endif

      h5g_close, grp_id
      h5f_close, file_id
      
      fnr++
      skip += nloc

   endrep until fnr eq NumGroupFiles



   
   return, data
end

