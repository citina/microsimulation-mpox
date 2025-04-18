% Get the row corresponding to the demographic group
% The common_col_struct provides all the fields we want to match
function demog_mat_indices = find_demog_rows(target_table, target_table_col_struct, demog_group_def, demog_col_struct, common_col_struct)

    % Account for the two types of age paradigms (age range with min/max,
    % and just age as an integer)
    
    demog_mat_indices = 1:size(target_table,1);
    
    if isfield(target_table_col_struct, 'age')
        for field = fieldnames(common_col_struct).'
            field_string = string(field);
            if strcmp(field_string,'age_min')
                demog_group_age_min = demog_group_def(demog_col_struct.age_min);
                demog_mat_indices = find_indices(target_table, demog_mat_indices, target_table_col_struct.age, '>=', demog_group_age_min); 
            elseif strcmp(field_string, 'age_max')
                demog_group_age_max = demog_group_def(demog_col_struct.age_max);
                demog_mat_indices = find_indices(target_table, demog_mat_indices, target_table_col_struct.age, '<=', demog_group_age_max); 
            else
                demog_field_value = demog_group_def(demog_col_struct.(field_string));
                demog_mat_indices = find_indices(target_table, demog_mat_indices, target_table_col_struct.(field_string), '=', demog_field_value); 
            end
        end
    else
        for field = fieldnames(common_col_struct).'
            field_string = string(field);
            demog_field_value = demog_group_def(demog_col_struct.(field_string));
            demog_mat_indices = find_indices(target_table, demog_mat_indices, target_table_col_struct.(field_string), '=', demog_field_value); 
        end
    end
end


% in target_table_col_struct => age
% in demog_col_struct => age_min, age_max