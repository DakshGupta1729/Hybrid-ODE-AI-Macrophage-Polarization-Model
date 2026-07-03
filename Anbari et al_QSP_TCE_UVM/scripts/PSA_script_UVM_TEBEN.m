%% In silico Virtual Clinical Trial (Paramster Sensitivity Analysis PSA)
% Script for setting up and running in silico clinical trial
clear
close all
sbioreset

%% Create the model
immune_oncology_model_UVM_TEBEN

%% Define dosing
% dose_schedule = [];
dose_schedule = schedule_dosing({'tebentafusp20','tebentafusp30','tebentafusp68'});

%% Define and Prepare the Input and Output Parameters

% % Generate New Parameter Sets (New Virtual Patient Cohort)
n_PSA = 1500;

% % Set distributions of the selected parameters for random sampling
params_in  = PSA_param_in_UVM_TEBEN;
% Set boundaries for model species (only include those added in the present model)
params_out = PSA_param_out(model);
% Add selected model outputs to sensitivity analysis
params_in  = PSA_param_obs(params_in);
% Randomly generate parameter sets using Latin-Hypercube Sa5-mpling
params_in  = PSA_setup(model,params_in,n_PSA);

% Save Virtual Patient Cohort
save('VP1500.mat', 'params_in', 'params_out')

% Run Batch Simulations
warning('off','all')
dbstop if warning
dbclear all

sbioaccelerate(model, dose_schedule)
tic
[simDataPSA, params_out] = simbio_PSA(model,params_in,params_out,dose_schedule);
toc

%% Postprocess
% Postprocess Data -> Calculate Clonality, Percentages and ...
simDataPSApost = PSA_post(simDataPSA,params_in,params_out);

% Add pre-treatment observables to the params_in
params_in = PSA_preObs(simDataPSA,simDataPSApost,params_in,params_out);

% Prepare the data for the rest of the analysis
params_out = PSA_prep(simDataPSA,simDataPSApost,params_out);

% Save and print data of interest (by assigning a unique code name for the trial)
% sprint_data(simDataPSA, simDataPSApost, params_in, params_out, 'teben_1500')

%% Perform and plot different types of analysis

% Partial Rank Correlation Coefficients
% PSA_PRCC(params_in,params_out,'plausible')
% PSA_PRCC(params_in,params_out)

% t-SNE Analysis
% PSA_tSNE(params_in,params_out,'plausible')
% PSA_tSNE(params_in,params_out,'patient')

% eFAST

% Principle Component Analysis

%% Plot Results
% close all
% Plot percent change in size and RECIST
PSA_plot_RECIST(simDataPSA,simDataPSApost,params_out)

% Tumor Size
% PSA_plot_TumSize(simDataPSA,simDataPSApost,params_out)

% Kaplan-Meier Progression-free survival (PFS)
% PSA_plot_KaplanMeier(simDataPSA,simDataPSApost,params_out)

%% Waterfall plots
% Plot percent change in size using waterfall plots for a parameter.
% Waterfall plots can be plotted using either end tumor sizes or best overall responses.
%
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'n_T1_clones')
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'k_P1_d1')
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'initial_tumour_diameter')
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'k_C1_growth')
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'k_Treg_TCE')
% PSA_plot_Waterfall(simDataPSApost,model,params_in,params_out,'k_T1')



% PSA_plot_Waterfall_color(simDataPSApost,model,params_in,params_out,[.2 .4 .7]) % blue
% PSA_plot_Waterfall_color(simDataPSApost,model,params_in,params_out,[.2 .6 .4]) % green
% PSA_plot_Waterfall_color(simDataPSApost,model,params_in,params_out,[.8 .6 .2]) % yellow/orange
