function [updated_state_matrix, transitionTally] = transition5(transition, transition_path,...
                                    state_matrix, StateMatCols, state_mat_demog_group_idx,...
                                    future_state_param_names, calibrationConstant)
    
    % make a copy of the state matrix such that a tally can be made of the number of changes that occur
    state_matrix_copy = state_matrix;

    % Read in transition matrix
    [state_trans_mat, StateTransCols] = read_table(transition_path);
    
    % Find people in the starting state within the demographic group
    eligible_rows_demog = state_mat_demog_group_idx;
    
    % iterate through the rows of the transition matrix
    for initialState = 1:size(state_trans_mat,1)
        
        % reset possible eligible rows for each iteration of a new inital state
        eligible_rows = eligible_rows_demog;
          
        % iterate through the field names in the transition
        % use this to find the final set of eligible rows
        for field = fieldnames(StateTransCols)'
            field_string = string(field);
            if(~contains(field_string,'state') && ~contains(field_string,'prob'))
                start_state_value = state_trans_mat(initialState, StateTransCols.(field_string));
                % check if fieldname contains age_min or age_max
                if(contains(field_string,'age_min'))
                    eligible_rows = find_indices(state_matrix, eligible_rows, StateMatCols.age, '>=', start_state_value);
                elseif (contains(field_string,'age_max'))
                    eligible_rows = find_indices(state_matrix, eligible_rows, StateMatCols.age, '<=', start_state_value);
                else
                    eligible_rows = find_indices(state_matrix, eligible_rows, StateMatCols.(field_string), '=', start_state_value); 
                end
            end
        end


        % Get all new_state and new_prob index
        state_idx = contains(fieldnames(StateTransCols), future_state_param_names(1));
        prob_idx = contains(fieldnames(StateTransCols), future_state_param_names(2));

        % Assign the parameter 
        drawComp = state_trans_mat(initialState, prob_idx);
        newState = state_trans_mat(initialState, state_idx);   

        drawComp = drawComp *calibrationConstant;

        % create a random draw value for all eligible rows
        unif_prob = rand([1,length(eligible_rows)]);
        
        % iterate through each comparison of the cumulative frequencies to
        % determine 
        to_transition = drawComp >= unif_prob;

        % identify which transitions should occur
        state_mat_transition_rows = eligible_rows(to_transition);

        % transition people to appropriate stae
        state_matrix(state_mat_transition_rows, StateMatCols.(transition)) = newState;
            
    end

    A = state_matrix_copy(:, StateMatCols.(transition)) ~= state_matrix(:, StateMatCols.(transition));


    transitionTally = sum(A);
    
    updated_state_matrix = state_matrix;
end
      
