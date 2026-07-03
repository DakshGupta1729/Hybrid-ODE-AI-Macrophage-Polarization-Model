% SimBiology Diagnostic Plot Generator
%
% Plots SimBiology Data
%
% Inputs: simData - SimBiology output data
%         model   - SimBiology model object


function diagnostic_plot(simData, model)

numClones = howManyClones(simData);
numMHC    = howManyMHC(simData);

f = figure;
% % Set the window in a specific size
% set(f0100,'Position', [50 50 1300 800]);
% Maximized window
set(f,'units','normalized','outerposition',[0 0 1 1])

%% Tumor growth and antigens
subplot(3,5,1); hold on; box on;
simbio_plot(simData,'C1');
 xlabel('Time (days)'); ylabel('Number of Cancer Cells'); 
%  set(gca, 'YScale', 'log'); %ylim([1 1e12]);

subplot(3,5,2); hold on; box on;
yyaxis left
simbio_plot(simData,'V_T');
xlabel('Time (days)'); ylabel('Tumour Volume'); 
% set(gca, 'YScale', 'log'); %ylim([1 1e12]);
yyaxis right
[~,V_T,~] = selectbyname(simData, 'V_T');
D_T = 2*(3/(4*pi)*V_T).^(1/3);
plot(simData.time, D_T)
ylabel('Tumour Size (cm)');



%% Teff in LN, Central, Tumor, and Peripheral
subplot(3,5,6); hold on; box on;
for i =1:numClones
simbio_plot(simData,['T',num2str(i)],'CompartmentName','V_LN','LegendEntry',['$T_{',num2str(i),',eff,LN}$']);
end
xlabel('Time (days)'); ylabel('Number of Cells'); legend; set(gca, 'YScale', 'log'); %ylim([1 1e12]);

subplot(3,5,7); hold on; box on;
for i =1:numClones
simbio_plot(simData,['T',num2str(i)],'CompartmentName','V_C','LegendEntry',['$T_{',num2str(i),',eff,C}$']);
end
xlabel('Time (days)'); ylabel('Number of Cells'); legend; set(gca, 'YScale', 'log'); %ylim([1 1e12]);

subplot(3,5,8); hold on; box on;
for i =1:numClones
simbio_plot(simData,['T',num2str(i)],'CompartmentName','V_T','LegendEntry',['$T_{',num2str(i),',eff,T}$']);
end
xlabel('Time (days)'); ylabel('Number of Cells'); legend; set(gca, 'YScale', 'log'); %ylim([1 1e12]);

subplot(3,5,9); hold on; box on;
for i =1:numClones
simbio_plot(simData,['T',num2str(i)],'CompartmentName','V_P','LegendEntry',['$T_{',num2str(i),',eff,P}$']);
end
xlabel('Time (days)'); ylabel('Number of Cells'); legend; set(gca, 'YScale', 'log'); %ylim([1 1e12]);

%% Treg in LN, Central, Tumor, and Peripheral
subplot(3,5,10); hold on; box on;
simbio_plot(simData,'T0','CompartmentName','V_T' ,'LegendEntry','$T_{reg,T}$' );
simbio_plot(simData,'T0','CompartmentName','V_C' ,'LegendEntry','$T_{reg,C}$' );
simbio_plot(simData,'T0','CompartmentName','V_LN','LegendEntry','$T_{reg,LN}$');
simbio_plot(simData,'T0','CompartmentName','V_P','LegendEntry','$T_{reg,P}$');
xlabel('Time (days)'); ylabel('Number of Cells'); legend; set(gca, 'YScale', 'log'); %ylim([1 1e12]);


%% Checkpoints and drugs
subplot(3,5,11); hold on; box on;
simbio_plot(simData,'gp100_teben_CD3','compartment','TCEsyn_T1_C1','legend','gp100_teben_CD3');
xlabel('Time (days)'); ylabel('molecule/micrometer^2');legend;
subplot(3,5,12); hold on; box on;
simbio_plot(simData,'gp100_teben','compartment','TCEsyn_T1_C1','legend','gp100-teben');
xlabel('Time (days)'); ylabel('molecule/micrometer^2');legend;
subplot(3,5,13); hold on; box on;
simbio_plot(simData,'CD3_teben','compartment','TCEsyn_T1_C1','legend','CD3-teben');
xlabel('Time (days)'); ylabel('molecule/micrometer^2');legend;
%% Checkpoints and drugs
subplot(3,5,14); hold on; box on;
 simbio_plot(simData,'Mac_M1','CompartmentName','V_T' ,'LegendEntry','$MAC_M1$');
xlabel('Time (days)'); ylabel('Cell');legend;
subplot(3,5,15); hold on; box on;
 simbio_plot(simData,'Mac_M2','CompartmentName','V_T' ,'LegendEntry','$MAC_M2$')
xlabel('Time (days)'); ylabel('cell');legend;
