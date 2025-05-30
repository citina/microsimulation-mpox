function [demog_table, DemogTblCols] = create_demog_groups(demog_var_def_file)
    
    % Read in csv file that defines demographic variable categories
    % age min and max MUST be listed first in excel sheet
    [demog_var_def, DemogTblCols] = read_table(demog_var_def_file);
        
    % Create possible age_min, age_max pairs if age is a factor
    if isfield(DemogTblCols,'age_min') && isfield(DemogTblCols,'age_max')
        age_buckets = [demog_var_def(:,DemogTblCols.age_min), demog_var_def(:,DemogTblCols.age_max)];
        
        % Template of variables we've already created combinations for
        template = age_buckets;
    else
        % template is blank if we dont have age buckets
        template = [];
    end

    % Placeholder matrix where we can build on top of the template 
    placeholder = [];
    
    % For every variable that isn't age
    for variable = fieldnames(DemogTblCols)'
        variable_string = string(variable);

        if ~strcmp(variable_string, 'age_min') && ~strcmp(variable_string, 'age_max')

            column = demog_var_def(:,DemogTblCols.(variable_string));
            column = column(~isnan(column));
            
            % Number of times to repeat the category is based on number of
            % buckets already created
            num_repeats = size(template, 1);
            
            % if number of repeats is zero, the number of repeats should
            % be reset to have a value of 1 (value of zero just means no
            % buckets created yet)
            if num_repeats == 0
                num_repeats = 1;
            end
            
            % For every row in each column
            for categ = column'
                to_append = [template repmat(categ, num_repeats, 1)];
                placeholder = [placeholder; to_append]; %#ok<AGROW>
            end 
            
            template = placeholder;
            placeholder = [];
        end
    end
    
    demog_table = template;

end
