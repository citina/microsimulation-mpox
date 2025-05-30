% This function will take in the current state matrix, and column defs,
% and the initial age and race proportion definitions  and proportions
% that are used to define the inital population, and the inflow proportion

function [updated_state_matrix, birthTally] = birth2(state_matrix, StateMatCols, age_def, race_def, race_prop, inflow)

	age_buckets = age_def;

	% check to see if inflow is proportion of total population alive or fixed number
	if inflow < 1
		numPeople = round(inflow * sum(state_matrix(:,StateMatCols.alive)));
	else
		numPeople = inflow;
	end

	% generate age attribute
	age = [];
    
	% index
	i = 1;
    
	% iterate through age buckets
	while i < length(age_buckets)

		% define min and max range
		age_min = age_buckets(i);
		age_max = age_buckets(i+1);

		popSize = round(numPeople);
		age = [age; randi([age_min, age_max], popSize,1)];

		i = i + 2;
	end

	% generate race attribute
	race = [];
    
	% iterate through all the race_index
	for i = 0:race_def(end) 
		% generate proper proportions for each race
		propSize = round(race_prop(i+1)*numPeople);
		race = [race; zeros(propSize,1)+i];
	end

	% because proportions may not be the same, find the smaller of the 2 and use that as the number 
	% of births
	numPeople = min(length(age), length(race));
	age = age(1:numPeople,1);
	race = race(1:numPeople,1);

	% treatment, aware, prep, and status will all be at the null (0)
	% alive will be 1
    [age_wk ,hiv_status, hiv_aware,...
     pox_status, pox_aware, vaccinated, vax_wk, vax2_wk,...
     treatment, isolation, added, ve, infectious_wks] = deal(zeros(numPeople,1));
    hiv_age_bucket = ones(numPeople,1);
    age_bucket = ones(numPeople,1);
    alive = ones(numPeople,1);
    
	birthMatrix = [age, age_wk, race, hiv_status, hiv_aware, pox_status, pox_aware, vaccinated, vax_wk, vax2_wk, treatment, isolation, alive, age_bucket, hiv_age_bucket, added, ve, infectious_wks];

	updated_state_matrix = [state_matrix ; birthMatrix];
	birthTally = numPeople;

% 	figure
% 	histogram(age,age_buckets)

end