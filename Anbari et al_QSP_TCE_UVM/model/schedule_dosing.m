% Dose Schedule
%
% Generates Donsing Schedule
%
% Inputs: drugName -- drug name or character array of drug names
%         varargin -- Name-Value Pairs
%                     - drugName_dose (mg/kg)
%                     - drugName_schedule [start,interval,repeat]
%                     - patientWeight (kg)
%
%       vaild drugName values: nivolumab, atezolizumab, ipilimumab, pembrolizumab, entinostat, nab-paclitaxel, tremelimumab
%
% Outputs: dosing -- SimBiology model object with new antigen module


function dose_schedule = schedule_dosing(drugName,varargin)

% Check if drugName is cell array
if (iscell(drugName))
    N = length(drugName);
else
    N = 1;
    drugName = {drugName};
end

% Optional Inputs
in = inputParser;
% Pembrolizumab
addParameter(in,'pembrolizumab_dose',2.5); % 200 mg every three weeks
addParameter(in,'pembrolizumab_schedule',[0,21,30]);
% Nivolumab
addParameter(in,'nivolumab_dose',20); % 20 mg/kg every four weeks
addParameter(in,'nivolumab_schedule',[0,28,30]);
% Atezolizumab
addParameter(in,'atezolizumab_dose',15); % 15 mg/kg every three weeks
addParameter(in,'atezolizumab_schedule',[0,21,30]);
% Durvalumab
addParameter(in,'durvalumab_dose',10); % 10 g/kg every two weeks
addParameter(in,'durvalumab_schedule',[0,14,30]);
% Ipilimumab
addParameter(in,'ipilimumab_dose',1); % 1 mg/kg every three weeks
addParameter(in,'ipilimumab_schedule',[0,21,30]);
% Tremelimumab
addParameter(in,'tremelimumab_dose',1); % 1 mg/kg every four weeks
addParameter(in,'tremelimumab_schedule',[0,28,30]);
% Entinostat
addParameter(in,'entinostat_dose',3); % 3 mg every week
addParameter(in,'entinostat_schedule',[0,7,55]);
% NabPaclitaxel
addParameter(in,'nabPaclitaxel_dose',100); % 100 mg/m2 Q3/4W
addParameter(in,'nabPaclitaxel_schedule',[0,28,15]);
% Cibisatamab
addParameter(in,'cibisatamab_dose',1.43); % 1.43 0.86 0.57 2.29 mg/kg every 7 days
addParameter(in,'cibisatamab_schedule',[0,7,70]);
% aCD47
addParameter(in,'aCD47_dose',10); % 10 mg/kg in combo; 15 in mono
addParameter(in,'aCD47_schedule',[0,7,57]);
% Tebentafusp
addParameter(in,'tebentafusp20_dose',0.00027); 
addParameter(in,'tebentafusp20_schedule',[0,7,0])
addParameter(in,'tebentafusp30_dose',0.0004); 
addParameter(in,'tebentafusp30_schedule',[7,7,0]); 
addParameter(in,'tebentafusp68_dose',0.00091); 
addParameter(in,'tebentafusp68_schedule',[14,7,70]); 

% Patient Weight
addParameter(in,'patientWeight',75);
% Patient Body Surface area
addParameter(in,'patientBSA',1.9);

% Parse Inputs
parse(in,varargin{:});

% Pembrolizumab
dose_pembro = in.Results.pembrolizumab_dose;
schedule_pembro = in.Results.pembrolizumab_schedule;
% Nivolumab
dose_nivo = in.Results.nivolumab_dose;
schedule_nivo = in.Results.nivolumab_schedule;
% Atezolizumab
dose_atezo = in.Results.atezolizumab_dose;
schedule_atezo = in.Results.atezolizumab_schedule;
% Durvalumab
dose_durva = in.Results.durvalumab_dose;
schedule_durva = in.Results.durvalumab_schedule;
% Ipilimumab
dose_ipil = in.Results.ipilimumab_dose;
schedule_ipil = in.Results.ipilimumab_schedule;
% Tremelimumab
dose_treme = in.Results.tremelimumab_dose;
schedule_treme = in.Results.tremelimumab_schedule;
% Entinostat
dose_ENT = in.Results.entinostat_dose;
schedule_ENT = in.Results.entinostat_schedule;
% NabPaclitaxel
dose_nabp = in.Results.nabPaclitaxel_dose;
schedule_nabp = in.Results.nabPaclitaxel_schedule;
% Cibisatamab
dose_cibis = in.Results.cibisatamab_dose;
schedule_cibis = in.Results.cibisatamab_schedule;
% aCD47
dose_aCD47 = in.Results.aCD47_dose;
schedule_aCD47 = in.Results.aCD47_schedule;
% Tebentafusp
dose_teben20 = in.Results.tebentafusp20_dose;
schedule_teben20 = in.Results.tebentafusp20_schedule;
dose_teben30 = in.Results.tebentafusp30_dose;
schedule_teben30 = in.Results.tebentafusp30_schedule;
dose_teben68 = in.Results.tebentafusp68_dose;
schedule_teben68 = in.Results.tebentafusp68_schedule;


% Patient Weight
patient_weight = in.Results.patientWeight;
% Patient BSA
patient_BSA = in.Results.patientBSA;

% Pembrolizumab
MW_pembro = 1.49E8;
doseObj_pembro = sbiodose('pembro','Amount',patient_weight*dose_pembro/MW_pembro,'AmountUnits','mole','TargetName','V_C.aPD1');
doseObj_pembro.StartTime = schedule_pembro(1);
doseObj_pembro.Interval = schedule_pembro(2);
doseObj_pembro.TimeUnits = 'day';
doseObj_pembro.RepeatCount = schedule_pembro(3);
doseObj_pembro.Active = true;

% Nivolumab
MW_nivo = 1.436E8; % milligrams per mole
doseObj_nivo = sbiodose('nivo','Amount',patient_weight*dose_nivo/MW_nivo,'AmountUnits','mole','TargetName','V_C.aPD1');
doseObj_nivo.StartTime = schedule_nivo(1);
doseObj_nivo.Interval = schedule_nivo(2);
doseObj_nivo.TimeUnits = 'day';
doseObj_nivo.RepeatCount = schedule_nivo(3);
doseObj_nivo.Active = true;

% Entinostat
MW_ENT = 3.764085e5; % mg/mole
Bio = 0.18; % Fraction of dose through buccal absorption

doseObj_ENT_1 = sbiodose('ENT_1','Amount',Bio*dose_ENT/MW_ENT,'AmountUnits','mole','TargetName','V_C.inostat_Buccal');
doseObj_ENT_1.StartTime = schedule_ENT(1);
doseObj_ENT_1.Interval = schedule_ENT(2);
doseObj_ENT_1.TimeUnits = 'day';
doseObj_ENT_1.RepeatCount = schedule_ENT(3);
doseObj_ENT_1.Active = true;
doseObj_ENT_1.DurationParameterName = 'durP_inostat';

doseObj_ENT_2 = sbiodose('ENT_2','Amount',(1-Bio)*dose_ENT/MW_ENT,'AmountUnits','mole','TargetName','V_C.Dose_GI');
doseObj_ENT_2.StartTime = schedule_ENT(1);
doseObj_ENT_2.Interval = schedule_ENT(2);
doseObj_ENT_2.TimeUnits = 'day';
doseObj_ENT_2.RepeatCount = schedule_ENT(3);
doseObj_ENT_2.Active = true;
doseObj_ENT_2.LagParameterName = 'lagP_inostat';

% nabPaclitaxel
% MW = 853.9; % g/mol
doseObj_nabp_1 = sbiodose('nabp_1','Amount',dose_nabp*patient_BSA,'AmountUnits','milligram','TargetName','V_1.NabP');
doseObj_nabp_1.Rate = doseObj_nabp_1.Amount/30;
doseObj_nabp_1.RateUnits = 'milligram/minute';
doseObj_nabp_1.StartTime = schedule_nabp(1);
doseObj_nabp_1.Active = true;
doseObj_nabp_1.Interval = schedule_nabp(2);
doseObj_nabp_1.TimeUnits = 'day';
doseObj_nabp_1.RepeatCount = schedule_nabp(3);

doseObj_nabp_2 = sbiodose('nabp_2','Amount',dose_nabp*patient_BSA,'AmountUnits','milligram','TargetName','V_1.NabP');
doseObj_nabp_2.Rate = doseObj_nabp_2.Amount/30;
doseObj_nabp_2.RateUnits = 'milligram/minute';
doseObj_nabp_2.StartTime = schedule_nabp(1)+7;
doseObj_nabp_2.Active = true;
doseObj_nabp_2.Interval = schedule_nabp(2);
doseObj_nabp_2.TimeUnits = 'day';
doseObj_nabp_2.RepeatCount = schedule_nabp(3);

doseObj_nabp_3 = sbiodose('nabp_3','Amount',dose_nabp*patient_BSA,'AmountUnits','milligram','TargetName','V_1.NabP');
doseObj_nabp_3.Rate = doseObj_nabp_3.Amount/30;
doseObj_nabp_3.RateUnits = 'milligram/minute';
doseObj_nabp_3.StartTime = schedule_nabp(1)+14;
doseObj_nabp_3.Active = true;
doseObj_nabp_3.Interval = schedule_nabp(2);
doseObj_nabp_3.TimeUnits = 'day';
doseObj_nabp_3.RepeatCount = schedule_nabp(3);

% Atezolizumab
MW_atezo = 1.436E8; % milligrams per mole
doseObj_atezo = sbiodose('atezo','Amount',patient_weight*dose_atezo/MW_atezo,'AmountUnits','mole','TargetName','V_C.aPDL1');
doseObj_atezo.StartTime = schedule_atezo(1);
doseObj_atezo.Interval = schedule_atezo(2);
doseObj_atezo.TimeUnits = 'day';
doseObj_atezo.RepeatCount = schedule_atezo(3);
doseObj_atezo.Active = true;

% Durvalumab
MW_durva = 1.49E8; % milligrams per mole
doseObj_durva = sbiodose('durva','Amount',patient_weight*dose_durva/MW_durva,'AmountUnits','mole','TargetName','V_C.aPDL1');
doseObj_durva.StartTime = schedule_durva(1);
doseObj_durva.Interval = schedule_durva(2);
doseObj_durva.TimeUnits = 'day';
doseObj_durva.RepeatCount = schedule_durva(3);
doseObj_durva.Active = true;

% Cibisatamab
MW_cibis = 1.75E8; % milligrams per mole
doseObj_cibis = sbiodose('cibis','Amount',patient_weight*dose_cibis/MW_cibis,'AmountUnits','mole','TargetName','V_C.cibis');
doseObj_cibis.StartTime = schedule_cibis(1);
doseObj_cibis.Interval = schedule_cibis(2);
doseObj_cibis.TimeUnits = 'day';
doseObj_cibis.RepeatCount = schedule_cibis(3);
doseObj_cibis.Active = true;

% aCD47
MW_aCD47 = 7.8E7; % milligrams per mole
doseObj_aCD47 = sbiodose('aCD47','Amount',patient_weight*dose_aCD47/MW_aCD47,'AmountUnits','mole','TargetName','V_C.aCD47');
doseObj_aCD47.StartTime = schedule_aCD47(1);
doseObj_aCD47.Rate = doseObj_aCD47.Amount/60;
doseObj_aCD47.RateUnits = 'mole/minute';
doseObj_aCD47.Interval = schedule_aCD47(2);
doseObj_aCD47.TimeUnits = 'day';
doseObj_aCD47.RepeatCount = schedule_aCD47(3);
doseObj_aCD47.Active = true;

% Ipilimumab
MW_ipil = 1.486349E8; % milligrams per mole
doseObj_ipi = sbiodose('ipi','Amount',patient_weight*dose_ipil/MW_ipil,'AmountUnits','mole','TargetName','V_C.aCTLA4');
doseObj_ipi.StartTime = schedule_ipil(1);
doseObj_ipi.Interval = schedule_ipil(2);
doseObj_ipi.TimeUnits = 'day';
doseObj_ipi.RepeatCount = schedule_ipil(3);
doseObj_ipi.Active = true;

% Tremelimumab
MW_treme = 1.4638E8; % milligrams per mole
doseObj_treme = sbiodose('treme','Amount',patient_weight*dose_treme/MW_treme,'AmountUnits','mole','TargetName','V_C.aCTLA4');
doseObj_treme.StartTime = schedule_treme(1);
doseObj_treme.Interval = schedule_treme(2);
doseObj_treme.TimeUnits = 'day';
doseObj_treme.RepeatCount = schedule_treme(3);
doseObj_treme.Active = true;

% Tebentafusp
MW_teben = 7.7E7; % milligrams per mole
doseObj_teben20 = sbiodose('teben20','Amount',patient_weight*dose_teben20/MW_teben,'AmountUnits','mole','TargetName','V_C.teben');
doseObj_teben20.StartTime = schedule_teben20(1);
doseObj_teben20.Interval = schedule_teben20(2);
doseObj_teben20.TimeUnits = 'day';
doseObj_teben20.RepeatCount = schedule_teben20(3);
doseObj_teben20.Active = true;

doseObj_teben30 = sbiodose('teben30','Amount',patient_weight*dose_teben30/MW_teben,'AmountUnits','mole','TargetName','V_C.teben');
doseObj_teben30.StartTime = schedule_teben30(1);
doseObj_teben30.Interval = schedule_teben30(2);
doseObj_teben30.TimeUnits = 'day';
doseObj_teben30.RepeatCount = schedule_teben30(3);
doseObj_teben30.Active = true;

doseObj_teben68 = sbiodose('teben68','Amount',patient_weight*dose_teben68/MW_teben,'AmountUnits','mole','TargetName','V_C.teben');
doseObj_teben68.StartTime = schedule_teben68(1);
doseObj_teben68.Interval = schedule_teben68(2);
doseObj_teben68.TimeUnits = 'day';
doseObj_teben68.RepeatCount = schedule_teben68(3);
doseObj_teben68.Active = true;



% Dose Schedule Array
dose_schedule(N) = sbiodose('empty'); % preallocate array
for i = 1:N
    switch drugName{i}
        case 'nivolumab'
            dose_schedule(i) = doseObj_nivo;
        case 'pembrolizumab'
            dose_schedule(i) = doseObj_pembro;
        case 'atezolizumab'
            dose_schedule(i) = doseObj_atezo;
        case 'durvalumab'
            dose_schedule(i) = doseObj_durva;
        case 'ipilimumab'
            dose_schedule(i) = doseObj_ipi;
        case 'tremelimumab'
            dose_schedule(i) = doseObj_treme;
        case 'nabPaclitaxel'
            dose_schedule(i) = doseObj_nabp_1;
        case 'entinostat'
            dose_schedule(i) = doseObj_ENT_1;
        case 'aCD47'
            dose_schedule(i) = doseObj_aCD47;
        case 'cibisatamab'
            dose_schedule(i) = doseObj_cibis;
       case 'tebentafusp20'
            dose_schedule(i) = doseObj_teben20;
       case 'tebentafusp30'
            dose_schedule(i) = doseObj_teben30;
       case 'tebentafusp68'
            dose_schedule(i) = doseObj_teben68;
        otherwise
            error('No match for drug name');
    end
end

for i = 1:N
    switch drugName{i}
        case 'nabPaclitaxel'
            dose_schedule(end+1) = doseObj_nabp_2;
            dose_schedule(end+1) = doseObj_nabp_3;
        case 'entinostat'
            dose_schedule(end+1) = doseObj_ENT_2;
    end
end
