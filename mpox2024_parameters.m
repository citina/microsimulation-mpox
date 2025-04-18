
% inputWorkbook = 'Inputs_monkeypox.xlsx';

% Number of weeks to run the simulation
T = xlsread(sim_inputFile, 'SimDuration');

% inflow proportion (weekly)
inflow = xlsread(sim_inputFile, 'InflowProportion'); 

age_def = xlsread(sim_inputFile, 'AgeDef')';

% race def and prop
race_def = xlsread(sim_inputFile, 'RaceDef')';
race_prop = xlsread(sim_inputFile, 'RaceProp')';

% Infection Calibration Constant
infectCalib_1 = xlsread(sim_inputFile, 'InfectCalib_1'); %week 1
infectCalib_2 = xlsread(sim_inputFile, 'InfectCalib_2'); %week 2-x
infectCalib_3 = xlsread(sim_inputFile, 'InfectCalib_3'); %week x+1 and after
hivInfect = xlsread(sim_inputFile, 'HIVInfect'); % relative risk for HIV+

% Multiplier for increased force of infection for young
youngInfect = xlsread(sim_inputFile, 'YoungInfect');
midInfect = xlsread(sim_inputFile, 'MidInfect');
oldInfect = xlsread(sim_inputFile, 'OldInfect');

% Multiplier for increased force of infection based on race
raceInfect = xlsread(sim_inputFile, 'RaceInfect');

% qurantine adherence
isolation_adherence = xlsread(sim_inputFile, 'isolationAdh');

% vaccination efficiency
vac1_plwh = 0.51;
vac2_plwh = 0.702;
vac1_normal = 0.721;
vac2_normal = 0.878;

[~, filePaths] = xlsread(sim_inputFile, 'FilePaths', 'B1:B19');

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

% Paths defined for isolation transition
iso_transition_path = '../input/isolationOn_0.2.csv';%filePaths{19}; 

% % Paths for policies
% vax_05x = filePaths{20}; 
% vax_2xblack = filePaths{21}; 
% vax_2xhis = filePaths{22}; 
% vax_2xplwh = filePaths{23}; 
% vax_novax = filePaths{31};
% vac1_wk1_plwh = filePaths{32};
% vac1_wk2_plwh = filePaths{33};
% vac1_wk3_plwh = filePaths{34};
% vac1_wk9_plwh = filePaths{35};
% vac1_wk1_hiplwh = filePaths{36};
% vac1_wk2_hiplwh = filePaths{37};
% vac1_wk9_hiplwh = filePaths{38};
% iso_05 = filePaths{24}; 
% vac1_wk1_2xb = filePaths{39};
% vac1_wk2_2xb = filePaths{40};
% vac1_wk3_2xb = filePaths{41};
% vac1_wk9_2xb = filePaths{42};
% vac1_wk1_2xh = filePaths{43};
% vac1_wk2_2xh = filePaths{44};
% vac1_wk3_2xh = filePaths{45};
% vac1_wk9_2xh = filePaths{46};
% vac1_wk1_2xw = filePaths{47};
% vac1_wk2_2xw = filePaths{48};
% vac1_wk3_2xw = filePaths{49};
% vac1_wk9_2xw = filePaths{50};
% 
% vac1_wk1_double = filePaths{51};
% vac1_wk2_double = filePaths{52};
% vac1_wk3_double = filePaths{53};
% vac1_wk9_double = filePaths{54};


% % Paths for sensitivity analysis
% asym2sym_transition_path_3 = filePaths{26};
% asym2sym_transition_path_5 = filePaths{27};
% asym2sym_transition_path_9 = filePaths{28};

% Define future state parameter names
future_state_param_names = {'new_state', 'new_prob'};










