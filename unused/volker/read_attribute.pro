function read_attribute, grp_id, attr_name

    attr_id = h5a_open_name(grp_id, attr_name)
    attr = h5a_read(attr_id)
    h5a_close, attr_id

    return, attr
end

