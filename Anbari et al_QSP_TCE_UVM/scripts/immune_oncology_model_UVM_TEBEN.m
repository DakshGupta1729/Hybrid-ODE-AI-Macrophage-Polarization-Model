%% Immune Oncology Model
% Script for setting up and running the immune oncology model in simbiology
clear
close all
sbioreset

%% Add use-defined units into both Simbiology and Symbolic libraries
% Add 'cell' unit to SimBiology and Symbolic Toolboxes
if (isempty(sbioshowunits('cell')))
    cell_unit = sbiounit('cell','molecule');
    sbioaddtolibrary(cell_unit);
end
% Symbolic Unit
u = symunit;
try u.cell;
catch
    newUnit('cell',u.molecule);
end

% Add 'mU' unit to SimBiology and Symbolic Toolboxes
mU_unit = sbiounit('mU','mole/liter');
sbioaddtolibrary(mU_unit);
% Symbolic Unit
u = symunit;
try u.mU;
catch
    newUnit('mU',u.molarity);
end

%% Setup Parameters
% Setup Parameters
params_in     = parameters_UVM_TEBEN;
params_out    = load_parameters(params_in);

%% Create the SimBiology Model
% Model Settings
model_name = 'Immune Oncology Model';
start_time = 0.0; % [days]
time_step = 1; % [days] 0.01 days ~ 15 mins
end_time = 400; % [days]
absolute_tolerance = 1e-9;
relative_tolerance = 1e-6;
% solver = 'ode15s';
solver = 'sundials';
% Model Object
time = start_time:time_step:end_time;
model = simbio_init(model_name,time,solver,absolute_tolerance,relative_tolerance,params_out);
% Maximal simulation time
config = getconfigset(model);
set(config, 'MaximumWallClock', 120)
% set(config, 'MaximumNumberOfLogs', 1)
set(config.SolverOptions, 'AbsoluteToleranceScaling', false)

%% Add Modules to the Model
% Cancer Modules
model = cancer_module(model,'C1',params_out); 
% T cell Modules
model = Treg_module(model,params_out);
model = Teff_module(model,'1',params_out,{'C1'});
% APC Module
model = APC_module(model,params_out);
% Antigen Modules
antigenCP = create_antigen({'C1'},5.4e-13,'antigenID',0);
model = antigen_module(model,'0',params_out,antigenCP);
antigen   = create_antigen({'C1'},5.4e-13,'antigenID',1);
model = antigen_module(model,'1',params_out,antigen);

% PK Modules
params_aPD1    = pk_parameters('pembrolizumab');
params_aPDL1   = pk_parameters('durvalumab');
params_aCTLA4  = pk_parameters('tremelimumab');
params_teben = pk_parameters('tebentafusp');
model = pk_module(model,'aPD1',params_aPD1);
model = pk_module(model,'aPDL1',params_aPDL1,'n');
model = pk_module(model,'aCTLA4' ,params_aCTLA4);
model = pk_module(model,'teben' ,params_teben);
% Checkpoint Modules
model = checkpoint_module(model,params_out,'T','C1');
model = checkpoint_module(model,params_out,'T','APC');


% QSPIO-TNBC Modules
model = Th_module(model,params_out);
model = MDSC_module(model,params_out,{'C1'},'inostat',0,'drugName','entinostat');
model = macrophage_module(model,params_out,{'C1'},'aCD47',0); % PK module in development for aCD47

params_TCE_T1 = TCEs_parameters(params_out,'T1','C1');
params_TCE_T0 = TCEs_parameters(params_out,'T0','C1');
model = TCEs_module_treg(model,params_TCE_T1, 'T1','C1' ,params_teben);
model = TCEs_module_treg(model,params_TCE_T0, 'T0','C1' ,params_teben);
model = TCEs_module(model,params_out,'T1',{'C1'});

% Setup Dosing
% dose_schedule = [];
dose_schedule = schedule_dosing({'tebentafusp20','tebentafusp30','tebentafusp68'});
dbstop if warning

% Initialize and Run the Model (should run with realistic baseline parameters)
% (should be commented out when conducting in silico virtual clinical trial)
tic
[model,success,simDataInit] = initial_conditions(model);
% [model,success] = initial_conditions(model);
toc

% Generate a list of parameters and species for debug
modelComp = listModelComp(model);
%%
% Run Simulation
if (success)
    tic
    simData = sbiosimulate(model,[],[],dose_schedule);
    toc
else
    simData = simDataInit;
    disp('Tumour did not reach specified initial tumour diameter with current parameters');
end
% 
% % % Plots
% % % Plot diagnostics
% % if (success)
%    diagnostic_plot(simData,model);
% %    diagnostic_plot_H(simData);
% %    diagnostic_plot_KPR(simData,model);
% % % end
