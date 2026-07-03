% TCEs Module
%
% Models Bites Interactions
%
% Inputs: model       -- simbio model object with four compartments
%         params      -- object containing the default parameters
%         Tname       -- name of the T cell forming the checkpoint synapse
%         Cname       -- name of the cancer forming the synapse
%         teben_params -- object containing the teben PK parameters
% Outputs: model -- simbio model object with new TCEs module
%
% Created: April 22, 2019 (Huilin Ma)
% Last Modified: May 20, 2019 (Huilin Ma)

function model = TCEs_module_treg(model,params,Tname,Cname,teben_params)

% select the right compartment based on cancer or APC
if Cname(1)=='C'
    compDrug = model.Compartment(3);
    gamma = 'gamma_T';
elseif Cname(1)=='A'
    compDrug = model.Compartment(4);
    gamma = 'gamma_LN';
end

% Add the synapse compartment
compt = addcompartment(model,['TCEsyn_',Tname,'_',Cname],params.T_syn.Value,'CapacityUnits',params.T_syn.Units);
    set(compt,'Notes',['synapse comparment between ',Tname,' and ',Cname,' ', params.T_syn.Notes]);

% Determine if first call
first_call = true;
try % see if synapse exist
    p = addparameter(model,'T_syn' ,params.T_syn.Value ,'ValueUnits',params.T_syn.Units);
    set(p,'Notes',['Surface area of the synapse ' params.T_syn.Notes]);
catch
    first_call = false;
end

if first_call
% Add Pharmacokinetics
% model = pk_module(model,'teben',teben_params);

% Add Synapse size
p = addparameter(model,'dT_syn' ,params.dT_syn.Value ,'ValueUnits',params.dT_syn.Units);
    set(p,'Notes',['Distance between two cells in the synapse ' params.dT_syn.Notes]);

% Add surface areas
p = addparameter(model,'SA_Tcell' ,params.A_Tcell.Value ,'ValueUnits',params.A_Tcell.Units);
    set(p,'Notes',['Surface area of the T cell ' params.A_Tcell.Notes]);
p = addparameter(model,'SA_cell' ,params.A_cell.Value ,'ValueUnits',params.A_cell.Units);
    set(p,'Notes',['Surface area of the Cancer cell ' params.A_cell.Notes]);
p = addparameter(model,'SA_APC' ,params.A_APC.Value ,'ValueUnits',params.A_APC.Units);
    set(p,'Notes',['Surface area of the APC ' params.A_APC.Notes]);

% Add kon Values
kon = addparameter(model,'kon_CD3_teben',params.kon_CD3_teben.Value,'ValueUnits',params.kon_CD3_teben.Units);
    set(kon,'Notes',['kon of CD3-teben binding ' params.kon_CD3_teben.Notes]);
kon = addparameter(model,'kon_gp100_teben',params.kon_gp100_teben.Value,'ValueUnits',params.kon_gp100_teben.Units);
    set(kon,'Notes',['kon of gp100-teben binding ' params.kon_CD3_teben.Notes]);
kon = addparameter(model,'kon_CD3_teben2',params.kon_CD3_teben2.Value,'ValueUnits',params.kon_CD3_teben2.Units);
    set(kon,'Notes',['kon of CD3-teben binding ' params.kon_CD3_teben.Notes]);
kon = addparameter(model,'kon_gp100_teben2',params.kon_gp100_teben2.Value,'ValueUnits',params.kon_gp100_teben2.Units);
    set(kon,'Notes',['kon of gp100-teben binding ' params.kon_CD3_teben.Notes]);

% Add koff Values
koff = addparameter(model,'koff_CD3_teben' ,params.koff_CD3_teben.Value ,'ValueUnits',params.koff_CD3_teben.Units);
    set(koff,'Notes',['koff of CD3-teben binding ' params.koff_CD3_teben.Notes]);
koff = addparameter(model,'koff_gp100_teben',params.koff_gp100_teben.Value,'ValueUnits',params.koff_gp100_teben.Units);
    set(koff,'Notes',['koff of gp100-teben binding ' params.koff_gp100_teben.Notes]);

% Bites Hill parameters
p = addparameter(model,'gp100_CD3_50',params.gp100_CD3_50.Value,'ValueUnits',params.gp100_CD3_50.Units);
    set(p,'Notes',['gp100/CD3 concentration for half-maximal T cell inactivation ' params.gp100_CD3_50.Notes]);
p = addparameter(model,'n_gp100_CD3',params.n_gp100_CD3.Value,'ValueUnits',params.n_gp100_CD3.Units);
    set(p,'Notes',['Hill coefficient for gp100/CD3 half-maximal T cell inactivation ' params.n_gp100_CD3.Notes]);

end

% Checkpoint Expressions
% Check if T cell has defined before
first_Tcell_call = true;
try
    p = addparameter(model,[Tname,'_CD3_total'],params.([Tname,'_CD3']).Value,'ValueUnits',params.([Tname,'_CD3']).Units,'ConstantValue',false);
catch
    first_Tcell_call = false;
end
if first_Tcell_call
        set(p,'Notes',['concentration of CD3 on ',Tname,' cells ' params.([Tname,'_CD3']).Notes]);
end
% Check if Cancer or APC has defined before
first_Ccell_call = true;
try
   p = addparameter(model,[Cname,'_gp100_total'],params.([Cname,'_gp100']).Value,'ValueUnits',params.([Cname,'_gp100']).Units,'ConstantValue',false);
catch
    first_Ccell_call = false;
end
if first_Ccell_call
        set(p,'Notes',['number of gp100 molecules per ',Cname,' cell ' params.([Cname,'_gp100']).Notes]);
end

% Add Species
x = addspecies(compt,'gp100_teben',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of gp100-teben-gp100 complex');
x = addspecies(compt,'gp100_teben_CD3',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of gp100-teben-CD3 complex');
x = addspecies(compt,'CD3_teben',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of CD3-teben complex');
x = addspecies(compt,'gp100',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of gp100 in synapse');
x = addspecies(compt,'CD3',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of CD3 in synapse');



% Update Input Parameters
addrule(model,[compt.Name,'.CD3',' = '  ,Tname,'_CD3_total /SA_Tcell'   ] ,'initialAssignment');
if Cname(1)=='C'
       addrule(model,[compt.Name,'.gp100 = ' ,Cname,'_gp100_total /SA_cell'] ,'initialAssignment');
elseif Cname(1)=='A'
       addrule(model,[compt.Name,'.gp100 = ' ,Cname,'_gp100_total /SA_APC'] ,'initialAssignment');
end

%bivalent reactions
  R = addreaction(model,[compt.Name,'.gp100 <-> ',compt.Name,'.gp100_teben']);
    set (R, 'ReactionRate', ['kon_gp100_teben*(',compt.Name,'.gp100 * V_T.teben) - koff_gp100_teben*',compt.Name,'.gp100_teben']);
    set (R, 'Notes'       , 'binding and unbinding of gp100 to teben on surface in synapse');

 R = addreaction(model,[compt.Name,'.gp100_teben + ',compt.Name,'.CD3 <-> ',compt.Name,'.gp100_teben_CD3']);
    set (R, 'ReactionRate', ['1e-4*kon_CD3_teben2*(',compt.Name,'.CD3 * ',compt.Name,'.gp100_teben) -  koff_CD3_teben*',compt.Name,'.gp100_teben_CD3']);
    set (R, 'Notes'       , 'binding and unbinding of CD3 to gp100_teben on in synapse');
 R = addreaction(model,[compt.Name,'.CD3 <-> ',compt.Name,'.CD3_teben']);
    set (R, 'ReactionRate', ['kon_CD3_teben*(',compt.Name,'.CD3 * V_T.teben) -  koff_CD3_teben*',compt.Name,'.CD3_teben']);
    set (R, 'Notes'       , 'binding and unbinding of CD3 to teben on surface in synapse');
 R = addreaction(model,[compt.Name,'.CD3_teben + ',compt.Name,'.gp100 <-> ',compt.Name,'.gp100_teben_CD3']);
    set (R, 'ReactionRate', ['1e-4*kon_gp100_teben2*(',compt.Name,'.CD3_teben)*(',compt.Name,'.gp100)-  koff_gp100_teben*',compt.Name,'.gp100_teben_CD3']);
    set (R, 'Notes'       , 'binding and unbinding of gp100 to CD3_teben on in synapse');

    % Set PD1 Hill Function
p = addparameter(model,['H_gp100_',Cname,'_',Tname],0,'ValueUnits','dimensionless','ConstantValue',false);
    set(p,'Notes','Hill function of gp100 cancer cells in tumor');
addrule(model,['H_gp100_',Cname,'_',Tname,' = ((',compt.Name,'.gp100_teben_CD3)/gp100_CD3_50)^n_gp100_CD3/(((',compt.Name,'.gp100_teben_CD3)/gp100_CD3_50)^n_gp100_CD3 + 1)'],'repeatedAssignment');
