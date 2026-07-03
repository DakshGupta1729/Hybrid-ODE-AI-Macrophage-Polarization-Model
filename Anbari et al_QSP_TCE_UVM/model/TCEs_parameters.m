% Function to generate bites module parameters from physical parameters
%
% Inputs: params_in  -- object containing the default parameters 
%         Tname      -- name of the T cell forming the checkpoint synapse 
%         Cname      -- name of the cancer or APC cell forming the 
%                       checkpoint synapse 
%
% Output: params_out -- object containing parameters to be used for a
%                       checkpoint module
%
% Created: April 17th, 2019 (Huilin Ma)
% Last Modified: MAy 20, 2019 (Huilin)

function params_out = TCEs_parameters(params_in,Tname,Cname)

% PD1/PDL1/PDL2 rates
% kon Values
params_out.kon_gp100_teben  = params_in.kon_gp100_teben ;
params_out.kon_gp100_teben2  = params_in.kon_gp100_teben / params_in.dT_syn;
params_out.kon_CD3_teben  = params_in.kon_CD3_teben ;
params_out.kon_CD3_teben2  = params_in.kon_CD3_teben /params_in.dT_syn;
params_out.kon_gp100_teben.Notes  = ['kon for gp100-teben ' params_in.kd_gp100_teben.Notes];
params_out.kon_CD3_teben.Notes  = ['kon for CD3-teben ' params_in.kd_CD3_teben.Notes];
params_out.kint  = params_in.kint ;
% koff Values
params_out.koff_gp100_teben  = params_in.kon_gp100_teben  * params_in.kd_gp100_teben; 
params_out.koff_CD3_teben  = params_in.kon_CD3_teben  * params_in.kd_CD3_teben;
params_out.koff_gp100_teben.Notes  = ['calculated based on the measured kd and kon ' params_in.kd_gp100_teben.Notes];
params_out.koff_CD3_teben.Notes  = ['calculated based on the measured kd and kon ' params_in.kd_CD3_teben.Notes];

% Synapse size
params_out.T_syn   = params_in.T_syn;
params_out.dT_syn   = params_in.dT_syn;

% Surface Area of Cells
params_out.A_Tcell = 4*pi*(params_in.D_Tcell/2)^2;
params_out.A_Tcell.Notes = ['calculated based on the average T cell diameter ' params_in.D_Tcell.Notes];
% Surface Area of Cancer Cells
params_out.A_cell = 4*pi*(params_in.D_cell/2)^2;
params_out.A_cell.Notes = ['calculated based on the average Cancer cell diameter ' params_in.D_cell.Notes];
% Surface Area of APC Cells
params_out.A_APC = params_in.A_s;

% Checkpoint Expression
% Select the expression of the right cell type based on cancer or APC
if Tname(2)=='1'  
    params_out.([Tname,'_CD3']) = params_in.T8_CD3;
elseif Tname(2)=='0'
    params_out.([Tname,'_CD3']) = params_in.T4_CD3;
end
if Cname(1)=='C'  
    params_out.([Cname,'_gp100']) = params_in.T8_gp100;
elseif Cname(1)=='A'
    params_out.([Cname,'_gp100']) = params_in.APC_gp100;
end

% Hill Function Parameters
params_out.gp100_CD3_50 = params_in.gp100_CD3_50;
params_out.n_gp100_CD3 = params_in.n_gp100_CD3;
