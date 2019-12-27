pro get_group_len_and_offset, Base, Num, NumGroupFiles, Ngroups, GroupLenType, GroupOffsetType

  GroupLenType = get_group_field('GroupLenType', Base, Num, NumGroupFiles, Ngroups)
 
  ntypes = n_elements(GroupLenType(*,0))
 
   ;; create the offset-table                                                                                                                                             
  GroupOffsetType = lon64arr(ntypes, Ngroups)

  for i=1L, Ngroups-1 do begin
     GroupOffsetType[*,i] =  GroupOffsetType[*,i-1] + GroupLenType[*,i-1]
  endfor


end

