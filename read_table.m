function [matrix, col_idx_struct] = read_table(table_file)
    input_struct = importdata(table_file);
    matrix = input_struct.data;
    cols = input_struct.colheaders;
    col_idx_struct = col_idx_to_name(cols);
end