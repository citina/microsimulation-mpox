% This function takes in a state matrix and outputs all row indices that
% 1) are in the eligible_people array
% 2) meet the criteria that its column variable is in the valid_values
% array
%
% Input:
% - state_matrix : matrix of (people, features)
% - eligible_people : vector of eligible row indices in state_matrix
% - column_idx : a column index specifying what variable to restrict on
% - comparison_type : defines what comparison is applied
% - valid_values : vector of valid values for the column specified by
% column_idx
%
% Output:
% - indices: a vector of row indices that are within eligible_people and
% meet the condition specified by column_idx and valid_values
%
function mat_indices = find_indices(state_matrix, eligible_people, column_idx, comparison_type, valid_values)
    
    % Set minimum and max in case of a between comparison_type
    if ismember(comparison_type, ['[]', '[)', '(]', '()'])
        if length(valid_values) ~= 2
            msg = 'Please enter a minimum and maximum value for a "between" comparison_type formatted as [min,max]';
            error(msg);
        end
        min = valid_values(1);
        max = valid_values(2);
    end
    switch comparison_type
        case '='
            relative_indices = find(ismember(state_matrix(eligible_people,column_idx), valid_values));
        case '[]'
            relative_indices = find(state_matrix(eligible_people,column_idx) >= min & state_matrix(:,column_idx) <= max);
        case '[)'
            relative_indices = find(state_matrix(eligible_people,column_idx) >= min & state_matrix(:,column_idx) < max);
        case '(]'
            relative_indices = find(state_matrix(eligible_people,column_idx) > min & state_matrix(:,column_idx) <= max);
        case '()'
            relative_indices = find(state_matrix(eligible_people,column_idx) > min & state_matrix(:,column_idx) < max);   
        case '<='
            relative_indices = find(state_matrix(eligible_people,column_idx) <= valid_values(1));
        case '<'
            relative_indices = find(state_matrix(eligible_people,column_idx) < valid_values(1));   
        case '>='
            relative_indices = find(state_matrix(eligible_people,column_idx) >= valid_values(1));
        case '>'
            relative_indices = find(state_matrix(eligible_people,column_idx) > valid_values(1));   
        otherwise
            msg = 'Please enter a valid comparison_type.';
            error(msg)
    end
    
    mat_indices = eligible_people(relative_indices);

end