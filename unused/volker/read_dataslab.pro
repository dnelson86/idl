function read_dataslab, file_id, dataset_name, dims, str, cnt

    dataset_id= h5d_open(file_id, dataset_name)

    dataspace_id = H5S_CREATE_SIMPLE(dims)
    h5s_select_hyperslab, dataspace_id, str, cnt ,/reset
    memoryspace_id = H5S_CREATE_SIMPLE(cnt)

    data  = h5d_read(dataset_id, FILE_SPACE=dataspace_id, MEMORY_SPACE=memoryspace_id)

    h5s_close, memoryspace_id
    h5s_close, dataspace_id
    h5d_close, dataset_id
   
    return, data
end

