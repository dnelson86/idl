
BaseDir =   "/virgo/simulations/Millennium/"    ; The output-directory of the simulation
SnapBase=   "snap_millennium"                   ; The basename of the snapshot files

Num = 63                          ; The number of the snapshot we look at

group = 1234567 ; group number

;;; The object-file of the compile C-library for accessing the group catalogue

ObjectFile = "/u/dnelson/idl/IdlGroupLib/idlgrouplib.so"

Num=long(Num)
exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)
Outputdir  = Basedir + "/snapdir_"+exts+"/"

;;;;;;;;;; First, we get the number of groups

Ngroups = CALL_EXTERNAL(ObjectFile, $
                       'get_total_number_of_groups', /UL_VALUE, $
                        OutputDir, $
                        Num)

print, "Number of groups in the catalogue: ", Ngroups

;;;;;;;;;; Now we load the group catalogue

GroupLen = lonarr(Ngroups)
GroupFileNr = lonarr(Ngroups)
GroupNr = lonarr(Ngroups)

Ngroups = CALL_EXTERNAL(ObjectFile, $
                       'get_group_catalogue', /UL_VALUE, $
                        OutputDir, $
                        Num, $
                        GroupLen, $
                        GroupFileNr, $
                        GroupNr)

dataset = {_NAME:'GroupLen', _TYPE:'Dataset', _DATA:GroupLen}
h5_create, 'test.hdf5', dataset

;;;;;;;;;; Now we get the size of the hashtable

NFiles=0L

HashTabSize = CALL_EXTERNAL(ObjectFile, $
                       'get_hash_table_size', /UL_VALUE, $
                        OutputDir, $
                        Num, $
                        SnapBase, $
                        NFiles)

;;;;;;;;; Now we load the hash-table

HashTable= lonarr(HashTabSize)
FileTable= lonarr(HashTabSize)
LastHashCell = lonarr(NFiles)
NInFiles = lonarr(NFiles)

result = CALL_EXTERNAL(ObjectFile, $
                       'get_hash_table', /UL_VALUE, $
                        OutputDir, $
                        Num, $
                        SnapBase, $
                        HashTable, $
                        FileTable, $
                        LastHashCell, $
                        NInFiles)


;;;;;; Now we are all set to read out individual groups (repeatedly if desired)

finr=  GroupFileNr(group) ; determines in which file this group is stored
grnr = GroupNr(group)     ; gives the group number within this file
Len = GroupLen(group)     ; gives the group length

Pos = Fltarr(3,Len)
Sx=0.0
Sy=0.0
Sz=0.0

result = CALL_EXTERNAL(ObjectFile, $
                     'get_group_coordinates', /UL_VALUE, $
                      OutputDir, $
                      Num, $
                      SnapBase, $
                      HashTable, $
                      FileTable, $
                      HashTabSize, $
                      LastHashCell, $
                      NInFiles, $ 
                      grnr, $
                      finr, $
                      Len, $
                      Pos, $
                      Sx, Sy, Sz)


 ;Plot, Pos(0,*), Pos(1,*), psym=3
 mean_x = mean(Pos[0,*])
 mean_y = mean(Pos[1,*])
 mean_z = mean(Pos[2,*])

 print, "GroupNr=", group, " length=", Len, "   center of mass="
 print, format='(f12)', sx, sy, sz
 print, "Mean x, y, z = ", mean_x, mean_y, mean_z

