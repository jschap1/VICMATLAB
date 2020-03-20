function dat = hdfreader_for_datastore(fname)

    
    dat = hdfread(fname, 'EOSGRID', 'Fields', 'T2M');
    
end