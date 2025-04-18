function col_idx_struct = col_idx_to_name(feature_names)
    col_idx_struct = struct() ;
    for i = 1 : length(feature_names)
        col = char(feature_names(1,i));
        col_idx_struct.(col) = i;
    end
end