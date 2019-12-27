function read_dataset, file_id, dataset_name

    dataset_id= h5d_open(file_id, dataset_name)
    data  = h5d_read(dataset_id)
    h5d_close, dataset_id
   
    return, data
end

