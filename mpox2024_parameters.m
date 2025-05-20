
% inputWorkbook = 'Inputs_monkeypox.xlsx';

% Number of weeks to run the simulation
T = xlsread(SIM_INPUT_FILE, 'SimDuration');

% inflow proportion (weekly)
inflow = xlsread(SIM_INPUT_FILE, 'InflowProportion'); 

age_def = xlsread(SIM_INPUT_FILE, 'AgeDef')';

% race def and prop
race_def = xlsread(SIM_INPUT_FILE, 'RaceDef')';
race_prop = xlsread(SIM_INPUT_FILE, 'RaceProp')';

% Infection Calibration Constant
infectCalib_1 = xlsread(SIM_INPUT_FILE, 'InfectCalib_1'); %week 1
infectCalib_2 = xlsread(SIM_INPUT_FILE, 'InfectCalib_2'); %week 2-x
infectCalib_3 = xlsread(SIM_INPUT_FILE, 'InfectCalib_3'); %week x+1 and after
hivInfect = xlsread(SIM_INPUT_FILE, 'HIVInfect'); % relative risk for HIV+

% Multiplier for increased force of infection for young
youngInfect = xlsread(SIM_INPUT_FILE, 'YoungInfect');
midInfect = xlsread(SIM_INPUT_FILE, 'MidInfect');
oldInfect = xlsread(SIM_INPUT_FILE, 'OldInfect');

% Multiplier for increased force of infection based on race
raceInfect = xlsread(SIM_INPUT_FILE, 'RaceInfect');

% qurantine adherence
isolation_adherence = xlsread(SIM_INPUT_FILE, 'isolationAdh');

% vaccination efficiency
vac1_plwh = 0.51;
vac2_plwh = 0.702;
vac1_normal = 0.721;
vac2_normal = 0.878;

[~, filePaths] = xlsread(SIM_INPUT_FILE, 'FilePaths', 'B1:B19');

% Initial population state matrix input file
init_pop_file = filePaths{1};

% File path for table that defines all demographic variable categories
demog_var_def_file = filePaths{2};

% File path for table that defines all the demographic variables for mixing matrix
mixing_mat_def_file = filePaths{3};

% File path for mixing matrix
mixing_mat_file = filePaths{4};

% File path for number of partners weekly based on mixing stratification;
num_partners_file = filePaths{5};

% Paths defined for natural death transition
deathNatural_transition_path = filePaths{7};

% Paths defined for vaccination transition
% vac1_wk1 = filePaths{8};
% vac1_wk2 = filePaths{9};
% vac1_wk3 = filePaths{10};
vac1_transition_path = filePaths{11};
vac2_transition_path = filePaths{12};

% Paths defined for awareness transition
aware1_transition_path = filePaths{13}; % for asymptomatic 
aware2_transition_path = filePaths{14}; % for symptomatic 

% Paths defined for mpox status transition
asym2recover_transition_path = filePaths{15};
asym2sym_transition_path = filePaths{16}; % for asymptomatic 
sym2recover_transition_path = filePaths{17}; % for symptomatic 

% Paths defined for get on treatment transition
ontrt_transition_path = filePaths{18}; 

% Define future state parameter names
future_state_param_names = {'new_state', 'new_prob'};










