% TCEs Module
%
% Inputs: model        -- SimBiology model object with four compartments
%         params       -- object containing model parameter Values, Units, and Notes:
%
% Outputs: model -- SimBiology model object with new TCEs module

function model = TCEs_module(model,params,Tname,cancer_types)

% Add Parameters
k_C_BTcell = addparameter(model,'k_C_BTcell',params.k_C_BTcell.Value,'ValueUnits',params.k_C_BTcell.Units,'ConstantValue',false);
    set(k_C_BTcell,'Notes',['Rate of cancer cell death by TCE activated T cells ' params.k_C_BTcell.Notes]);
k_C_BTreg = addparameter(model,'k_C_BTreg',params.k_C_BTreg.Value,'ValueUnits',params.k_C_BTreg.Units,'ConstantValue',false);
    set(k_C_BTreg,'Notes',['Rate of cancer cell death by TCE activated Tregs ' params.k_C_BTreg.Notes]);
k_Treg_TCE = addparameter(model,'k_Treg_TCE',params.k_Treg_TCE.Value,'ValueUnits',params.k_Treg_TCE.Units,'ConstantValue',false);
    set(k_Treg_TCE,'Notes',['Rate of T cell death by TCE activated Tregs ' params.k_Treg_TCE.Notes]);

% Add Reactions
reaction = addreaction(model,['V_T.' Tname ' -> V_T.T1_exh']);
    set(reaction,'ReactionRate',['k_Treg_TCE*V_T.' Tname '*Tregs_/(V_T.' Tname '+Tregs_+cell)*H_gp100_C1_T0']);
    set(reaction,'Notes','TCEs mediated T cell death from Tregs');

% Get Model Rules for Updating
model_rules = get(model,'Rules');

% Update TCEs Mediated Cancer Killing
for i = 1:length(cancer_types)
    reaction = addreaction(model,['V_T.' cancer_types{i} ' -> V_T.C_x']);
        set(reaction,'ReactionRate',['k_C_BTcell*V_T.' Tname '*V_T.' cancer_types{i} '/(K_T_C*V_T.' cancer_types{i} '+V_T.' Tname '+cell)*V_T.' Tname '/(V_T.' Tname '+K_T_Treg*Tregs_+cell)*H_gp100_C1_T1'...
            ' + k_C_BTreg*V_T.' Tname '*V_T.' cancer_types{i} '/(K_T_C*V_T.' cancer_types{i} '+V_T.' cancer_types{i} '+cell)*V_T.' cancer_types{i} '/(V_T.' cancer_types{i} '+K_T_Treg*Tregs_+cell)*H_gp100_C1_T0']);
        set(reaction,'Notes','TCEs Mediated Cancer cell killing by T cells');
    rule = get(model_rules(5),'Rule');
        set(model_rules(5),'Rule',[rule,' + k_C_BTcell*V_T.' Tname '*V_T.' cancer_types{i} '/(K_T_C*V_T.' cancer_types{i} '+V_T.' Tname '+cell)*V_T.' Tname '/(V_T.' Tname '+K_T_Treg*Tregs_+cell)*H_gp100_C1_T1'...
            ' + k_C_BTreg*V_T.' Tname '*V_T.' cancer_types{i} '/(K_T_C*V_T.' cancer_types{i} '+V_T.' cancer_types{i} '+cell)*V_T.' cancer_types{i} '/(V_T.' cancer_types{i} '+K_T_Treg*Tregs_+cell)*H_gp100_C1_T0']);

    for j = 1:length(model_rules)
        if ~isempty(strfind(model_rules(j).Rule, ['k_' cancer_types{i} '_therapy'])) && isempty(strfind(model_rules(j).Rule, 'k_C_BTcell'))
            model_rules(j).Rule = [model_rules(j).Rule ' + k_C_BTcell*V_T.' Tname '/(K_T_C*V_T.' cancer_types{i} '+V_T.' Tname '+cell)*V_T.' Tname '/(V_T.' Tname '+K_T_Treg*Tregs_+cell)*H_gp100_C1_T1'...
            ' + k_C_BTreg*V_T.' Tname '/(K_T_C*V_T.' cancer_types{i} '+V_T.' cancer_types{i} '+cell)*V_T.' cancer_types{i} '/(V_T.' cancer_types{i} '+K_T_Treg*Tregs_+cell)*H_gp100_C1_T0'];
        end
    end
end
