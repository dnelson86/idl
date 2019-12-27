function get_snap_field, Type, FieldName, Base, SnapBase, Num, NumSnapFiles, NumPartType

  if num ge 1000 then begin
     exts='0000'
     exts=exts+strcompress(string(Num),/remove_all)
     exts=strmid(exts,strlen(exts)-4,4)
  endif else begin
     exts='000'
     exts=exts+strcompress(string(Num),/remove_all)
     exts=strmid(exts,strlen(exts)-3,3)
  endelse
  
   skip = 0LL
   skipids=0LL
   skipsub=0LL
   fnr = 0L


   Len = NumPartType[Type]
   Off = 0LL

   ts = string(Type, format ='(I1)')

   nleft = Len

   if Len gt 0 then begin
                                ; go through the snapshot file and
                                ; filter out the particles we need
      fnr = 0L
      flag = 0

      repeat begin
         
         f = Base + "/snapdir_" + exts +"/"+snapbase+"_"+exts +"."+strcompress(string(fnr),/remove_all) + ".hdf5"

         if fnr eq 0 then begin
                                ; check whether file exists, otherwise
                                ; assume the data is in a single file
            openr, 1, f, error = err
            if err ne 0 then begin
               f = Base +"/"+snapbase+"_"+exts + ".hdf5"
            endif else begin
               close, 1
            endelse
         endif
            

            
         file_id = h5f_open(f)

         grp_id = h5g_open(file_id, 'Header')
            
         nlocType = read_attribute(grp_id, 'NumPart_ThisFile')
         
         nloc = nlocType[type]

         if (flag eq 0) and (nloc gt 0) then begin
  
            flag = 1

            dataset_id= h5d_open(file_id, 'PartType'+ts+'/' + FieldName)

            dataspace_id = H5D_GET_SPACE(dataset_id)
            datatype_id = H5D_GET_TYPE(dataset_id)
            
            idltype = H5T_IDLTYPE(Datatype_id)
            dimensions = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace_id)

          ;  print, "IdlType=", idltype, "  Dims=", dimensions
            
            h5d_close, dataset_id
            
            if n_elements(dimensions) eq 2 then begin
               firstdim = dimensions(0)
               if idltype eq 3 then begin
                  data = lonarr(firstdim, Len) 
               endif 
               if idltype eq 15 then begin
                  data = lon64arr(firstdim, Len) 
               endif 
               if idltype eq 4 then begin
                  data = fltarr(firstdim, Len) 
               endif 
               if idltype eq 5 then begin
                  data = dblarr(firstdim, Len) 
               endif 
            endif else begin
               firstdim = 1
               if idltype eq 3 then begin
                  data = lonarr(Len) 
               endif 
               if idltype eq 15 then begin
                  data = lon64arr(Len) 
               endif 
               if idltype eq 4 then begin
                  data = fltarr(Len) 
               endif 
               if idltype eq 5 then begin
                  data = dblarr(Len) 
               endif 
            endelse
         endif


         if nloc gt off then begin ; we may have something in this file
            nstart = off
            if nloc - off gt nleft then begin ; there are more particles in the file then we need
               ncount = nleft  
            endif else begin
               ncount = nloc - off
            endelse


            if firstdim gt 1 then begin

               dims  = [firstdim, nloc]
               start = [0, nstart]
               count = [firstdim, ncount]

               Data(*, skip:skip + ncount -1)  = read_dataslab(file_id, 'PartType'+ts+'/'+FieldName, dims, start, count)
            endif else begin

               Data(skip:skip + ncount -1)    = read_dataslab(file_id, 'PartType'+ts+'/'+FieldName, [nloc], [nstart], [ncount])

            endelse


            nleft -= ncount
            skip += ncount
            off += ncount
         endif
            
         off -= nloc
         
         h5g_close, grp_id
         h5f_close, file_id
         
         fnr++
            
      endrep until (fnr eq NumSnapFiles) || (nleft eq 0)
      
   endif
   
   return, data
end

