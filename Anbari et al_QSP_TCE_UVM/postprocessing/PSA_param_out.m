% Function to generate object with model parameters to include as outputs 
% of the parameter sensitivity analysis
%
% Output: params -- object containing output parameters
%                   -> for each parameter: 
%                       - Adds the name of the parameter to the list
%                       - defines upper and lower bounds of acceptable
%                       physiological range
%                       - specifies if the output is comaprtment,
%                       parameter, species or postprocessed parameters
%                           - for species: specifies compartment
%         Examples:      
%                 % k1
%                 params.names = [params.names; 'k1'];
%                 params.k1.UpperBound = 1; 
%                 params.k1.LowerBound = 0;
%                 params.k1.Units      = 'liter';
%                 params.k1.Type       = 'compartment';
%                 params.k1.ScreenName = 'k1 binding rate';
% 
%                 % k2
%                 params.names = [params.names; 'k2'];
%                 params.k2.UpperBound = 1; 
%                 params.k2.LowerBound = 0;
%                 params.k2.Units      = 'liter';
%                 params.k2.Type       = 'parameter';
%                 params.k2.ScreenName = 'k2 binding rate';
%
%                 % k3 in central
%                 params.names = [params.names; 'k3_C'];
%                 params.k3_C.UpperBound  = 1e12;
%                 params.k3_C.LowerBound  = 0;
%                 params.k3_C.Units       = 'cell';
%                 params.k3_C.Type        = 'species';
%                 params.k3_C.Name        = 'k2';
%                 params.k3_C.Compartment = 'V_C';
%                 params.k3.ScreenName    = 'k3 binding rate';
%
%                 % k4 post processed
%                 params.names = [params.names; 'k4'];
%                 params.k4.UpperBound  = 1;
%                 params.k4.LowerBound  = 0;
%                 params.k4.Units       = 'dimensionless';
%                 params.k4.Type        = 'post';
%                 params.k4.ScreenName  = 'k4 binding rate';
% NOTE: for postprocessed parameters ensure the name here is the same as 
%       the name assigned in the postprocess functions       

% Created: Jan 20, 2019 (Mohammad Jafarnejad)
% Last Modified: Jan 23, 2019 (MJ)

function params = PSA_param_out(varargin)

params.names = {};

% Tumor Volume
params.names = [params.names; 'V_T'];
params.V_T.UpperBound = 1e6;
params.V_T.LowerBound = 0.00;
params.V_T.Units      = 'microliter';
params.V_T.Type       = 'compartment';
params.V_T.ScreenName = 'Tumor Volume';


% Number of CD8 T cells in the tumor
params.names = [params.names; 'CD8_density'];
params.CD8_density.UpperBound  = 5e5;
params.CD8_density.LowerBound  = 0;
params.CD8_density.Units       = 'cell/milliliter';
params.CD8_density.Type        = 'pre';
params.CD8_density.ScreenName  = 'CD8 Density in Tumor';

% Number of Treg cells in the tumor
params.names = [params.names; 'Treg_density'];
params.Treg_density.UpperBound  = 5e5;
params.Treg_density.LowerBound  = 0;
params.Treg_density.Units       = 'cell/milliliter';
params.Treg_density.Type        = 'pre';
params.Treg_density.ScreenName  = 'Treg Density in Tumor';


% CD8 to Treg ratio
params.names = [params.names; 'CD8FoxP3ratio_T'];
params.CD8FoxP3ratio_T.UpperBound  = 10;
params.CD8FoxP3ratio_T.LowerBound  = 0.3;
params.CD8FoxP3ratio_T.Units      = 'dimensionless';
params.CD8FoxP3ratio_T.Type       = 'pre';
params.CD8FoxP3ratio_T.ScreenName = 'CD8 to Treg Ratio in Tumor';

% CD4 to Treg ratio
params.names = [params.names; 'CD4FoxP3ratio_T'];
params.CD4FoxP3ratio_T.UpperBound  = 20;
params.CD4FoxP3ratio_T.LowerBound  = 1;
params.CD4FoxP3ratio_T.Units      = 'dimensionless';
params.CD4FoxP3ratio_T.Type       = 'pre';
params.CD4FoxP3ratio_T.ScreenName = 'CD4 to Treg Ratio in Tumor';


% all Hill functions
if (nargin == 1)
    model = varargin{1};
    for i = 1:length(model.Parameters)
        if (model.Parameters(i).Name(1)=='H')
            params.names = [params.names; model.Parameters(i).Name];
            params.(model.Parameters(i).Name).UpperBound = 1;
            params.(model.Parameters(i).Name).LowerBound = 0;
            params.(model.Parameters(i).Name).Units      = model.Parameters(i).ValueUnits;
            params.(model.Parameters(i).Name).Type       = 'parameter';
            params.(model.Parameters(i).Name).ScreenName = strrep(model.Parameters(i).Name,'_','\_');
        end
    end
end
