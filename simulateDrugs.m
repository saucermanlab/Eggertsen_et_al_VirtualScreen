%% If any inputs have changed, run the following section: 
close all; clear all; clc;

model='Hypertrophy_Model.xlsx';

% load('finalDrugOutputNetworkTargets.mat'); % Result from the webscraper NEED TO CLEAN UP WEBSCRAPER OUTPUT AND THEN CHANGE THIS TO BE READTABLE LIKE THE OTHER VARIABLES
warning off;
formattedReactions = table;
% Species/Reaction information from toy_model.xlsx or network, 'species'/'reactions' tab 
networkReactions = readtable(model, 'Sheet', 'reactions');
if ismember('=',cell2mat(networkReactions{1,3}))==0
    networkReactions(1,:)=[];  
end
% Formats network reactions to only show the product/output. 
for i = 1:height(networkReactions)
    reaction = string(networkReactions{i,3});
    nodeOfReaction = extractAfter(reaction, '=>'); 
    formattedReactions{i,1} = strtrim(nodeOfReaction);
end
formattedReactions.Properties.VariableNames(1) = {'ReactionOutputNode'};
save('formattedReactions.mat', 'formattedReactions');

%% Generate the drug modified ODE and parameter files from NETFLUX
[status, result] = generateDrugODEs(model);

%% Create DrugsToSimulate.csv
% % % [status, result] = generateDrugsToSimulate(model);

%% Inputs for simulations
opts = detectImportOptions('DrugsToSimulate.csv');
opts = setvartype(opts,{'AgonistTargetIndex','AntagonistTargetIndex'},'char');
drugsToSimulate = readtable('DrugsToSimulate.csv');% 'Z:\Taylor\Code\PPI\PPI_Edges\DrugsToSimulate2.csv',opts);% 

formattedReactions = load('formattedReactions.mat');
formattedReactions = formattedReactions.formattedReactions; % Extract from struct
% formattedReactions(1,:) = []; % remove empty (first) row

%% Set drug dose or doses (as a vector) & Sensitivity analysis

alteration_antag = [0:0.02:1] ;%%% [0:0.02:1]
% % % [padding,list_numDrugs, list_drugAction, list_drugType, list_drugTarget, rowLabels_drugs] = createRowLabels(drugsToSimulate);

inputs = ["AngII","ANPi","BNPi","CT1","EGF","ET1","FGF","IGF1","IL6","ISO","LIF","NE","NRG1","PE","Stretch","TGFB","TNFa"];
inp = 16; %input to network 
drugofinterest = 146;  %146, 33, 52

% % % for inp = 1:length(inputs)%%% Deleteme %%%

% Node parameters
[params,y0] = tempDrugODE_params;
[rpar,tau,ymax,speciesNames]=params{:};

% Steady-state Control pre-Simulation
tspan = [0 50]; options = [];
[pt,py] = ode15s(@tempDrugODE,tspan,y0,options,params);

w = rpar(1,:);
%New input value
% % % % ymax(61) = 0.01; %map3k23
% % % % ymax(70) = 0.01; %mekk1
% % % % ymax(90) = 0.01; %rac1
% % % % ymax(9) = 0.01; %ATF2
% % % % ymax(22) = 0.01; %cJun
% % % % y0(46) = ymax(46); tau(46) = 10^9; %GSK3b
% % % % y0(39) = ymax(39); tau(39) = 10^9; %Foxo
w(inp) = 0.1;

n = rpar(2,:);
EC50 = rpar(3,:);
dose = rpar(4,:);
drugBinding = rpar(5,:); 
drugAgonism = rpar(6,:); 
rpar = [w;n;EC50;dose;drugBinding;drugAgonism];
inputNodeW = num2cell(1:length(speciesNames)); % Nodes to test drug against
params=[rpar,tau,ymax,speciesNames];
% Steady-state Control Simulation
tspan = [0 50]; options = [];
[t,y] = ode15s(@tempDrugODE,tspan,y0,options,params);
yEnd = y(end,:)';

%%%deleteme%%%
figure; hold on
leg = [];
for k = [19 93]
    plot(t,y(end,k))
    leg = [leg,speciesNames(k)];
end
legend(leg)
%%%deleteme%%%

% Reset initial y values
y0 = real(yEnd); 

pred = zeros(length(alteration_antag),2);
sens = zeros(height(drugsToSimulate),length(inputNodeW));
for i = drugofinterest-1 %input drug of interest
    drugsToSimulate.Drug{i}
    drugsToSimulate.AgonistTarget{i}
    drugsToSimulate.AntagonistTarget{i}
    for j = 1:length(alteration_antag) 
        [t2,y2] = drugtreatment(dose,drugBinding,drugAgonism,drugsToSimulate,formattedReactions,alteration_antag,w,n,EC50,tau,ymax,speciesNames,y0,i,j);
        ySimEnd = y2(end,:)';

        count = 1; leg = [];
        specset = [19 11 12 18 74];
        for k = specset 
             responseset(count,j) = y2(end,k);
             leg = [leg, speciesNames(k)];
             count = count + 1;
        end

    end

    figure; hold on
    for kk=1:length(specset)
        plot(alteration_antag,responseset(kk,:)); ylim([0 1]);
    end
    xlabel('Drug Dose'); ylabel('Fractional Species Activation')
    legend(leg)
    title(drugsToSimulate.Drug{i})
    set(gca,'Xscale','log')

    for ii = 1:length(speciesNames)
        targeteffect(ii) = y2(end,ii)-y2(1,ii);
    end
    figure
    Xt = categorical(speciesNames);
    Xt = reordercats(Xt,speciesNames);
    bar(Xt,targeteffect)
    title(strcat(drugsToSimulate.Drug{i},';',inputs(inp)))
    ylabel('Change in Activation')

    %Write data to text file
    X = categorical(speciesNames');
    X = reordercats(X,string(X));
    T_Fin = y2(end,:)'; T0 = y2(1,:)';
    T = table(X,T0,T_Fin,T_Fin-T0,'VariableNames',{'Species_name','T0','T_Fin','T_Diff'});
    filename = strcat(char(drugsToSimulate.Drug{i}),'_',char(inputs(inp)),'.txt');
    writetable(T, filename)
end
    
% % %     ajax(inp) = responseset(1,46)-responseset(1,1);%%% Deleteme %%%
% % % end %%% Deleteme %%%

%% Graphing the predicted outcomes

nn=categorical({'Control','NE','NE + Celecoxib','NE + GSK3b OE','NE + Foxo OE','NE + GSK3b/Foxo OE','NE + Celecoxib + GSK3b KD','NE + Celecoxib + Foxo KD'});
nn=reordercats(nn,{'Control','NE','NE + Celecoxib','NE + GSK3b OE','NE + Foxo OE','NE + GSK3b/Foxo OE','NE + Celecoxib + GSK3b KD','NE + Celecoxib + Foxo KD'});
A=[0.1593,0.3418,0.2174,0.2561,0.2961,0.2044,1,1]; %{'Control','NE','NE + Celecoxib','NE + GSK3b KD','NE + Foxo KD','NE + GSK3b/Foxo KD'}
% % % A=[0.1593,0.7077,0.1808,0.1232,0.6049]; %AngII/Baricitinib/Ras KD/PI3K KD
% % % A=[0.1593,0.7077,0.1894,0.6292,0.6235,0.5225]; %{'Control','AngII','AngII + Baricitinib','AngII + ATF2 KD','AngII + cJun KD','AngII + ATF2/cJun KD'}
% % % A=[0.1593,0.7077,0.1808,0.5485,0.5903,0.5409,0.2713]; %{'Control','AngII','AngII + Baricitinib','AngII + Mekk1 KD','AngII + Rac1 KD','AngII + Map3k23 KD','AngII + Mekk1/Rac1/Map3k23 KD'}
% % % A=[0.1593,0.6604,0.1954,0.08988,0.652]; %TGFb/Midostaurin/Ras KD/PKD KD
% % % A=[0.1593,0.6604,0.1954,0.4956,0.4871,0.3057,1,1]; %{'Control','TGFb','TGFb + Midostaurin','TGFb + Mekk1 KD','TGFb + Map3k23 KD','TGFb + Mekk1/Map3k23 KD'} 
% % % A=[0.1593,0.7623,0.0625,0.08411,0.6886]; %Nrg1/Zanubrutinib/EGFR KD/ERBB KD/JAK KD
% % % A=[0.1593,0.7623,0.0625,0.5868,0.6253,0.5776,0.2299]; %Nrg1/Zanubrutinib/Mekk1 KD/Rac1 KD/Map3k23 KD/Mekk1,Rac1,Map3k23 KD
% % % A=[0.1593,0.7077,0.6845,0.1808,0.6152,0.3418,0.3157,0.2383,0.2174,0.6604,0.1954,0.5371,0.6231]; %control/AngII/AngII+Midostaurin/AngII+Baricitinib/AngII+Celecoxib/NE/NE+Midostaurin/NE+Baricitinib/NE+Celecoxib/TGFb/TGFb+Midostaurin/TGFb+Baricitinib/TGFb+Celecoxib

figure
b=bar(nn,A);
b.FaceColor = '#3182bd';
ylabel('Predicted Cell Area (a.u.)')
xtickangle(0)

nn=categorical({'Control','TGFb','TGFb + Midostaurin','TGFb + Mekk1 KD','TGFb + Map3k23 KD','TGFb + Mekk1/Map3k23 KD','TGFb + Mekk1 OE','TGFb + Map3k23 OE'});
nn=reordercats(nn,{'Control','TGFb','TGFb + Midostaurin','TGFb + Mekk1 KD','TGFb + Map3k23 KD','TGFb + Mekk1/Map3k23 KD','TGFb + Mekk1 OE','TGFb + Map3k23 OE'});
A=[0.1593,0.6604,0.1954,0.4956,0.4871,0.3057,1,1]; %{'Control','TGFb','TGFb + Midostaurin','TGFb + Mekk1 KD','TGFb + Map3k23 KD','TGFb + Mekk1/Map3k23 KD'}
figure
b=bar(nn,A);
b.FaceColor = '#3182bd';
ylabel('Predicted Cell Area (a.u.)')
xtickangle(0)

namestimuli = char(inputs(inp));
namestimuli_drug = strcat(namestimuli,'+',char(drugsToSimulate.Drug{i}));
outcontrol_CA =py(end,19); outcontrol_BNP=0.134; outcontrol_MHC=0.1357;
outstimuli_CA=responseset(1,1); outstimuli_MHC=responseset(2,1); outstimuli_BNP=responseset(3,1);
outstimuli_drug_CA=responseset(1,41);outstimuli_drug_MHC=responseset(2,41); outstimuli_drug_BNP=responseset(3,41);

x=categorical({'Control', namestimuli, namestimuli_drug});
x= reordercats(x,{'Control', namestimuli, namestimuli_drug});
y=[outcontrol_CA, outstimuli_CA, outstimuli_drug_CA];
% % % % % y=[outcontrol_MHC, outstimuli_MHC, outstimuli_drug_MHC];
% % % % % y=[outcontrol_BNP, outstimuli_BNP, outstimuli_drug_BNP];

%Input experimental data %%
exp_control= 100.8 ;
exp_control_std= 13.6 ;
exp_stim= 132 ;
exp_stim_std= 11.2 ;
exp_stim_drug= 115.6 ;
exp_stim_drug_std= 8.8 ;
ye=[exp_control, exp_stim, exp_stim_drug];
yerr = [exp_control_std, exp_stim_std, exp_stim_drug_std];

figure 
% % % subplot(2,1,2)
% % % b2=bar(x,ye); hold on
% % % er = errorbar(x,ye,yerr,yerr); er.Color=[0 0 0]; er.LineStyle='none';
% % % title('Experiment')
% % % ylabel('Cell Area (\mum^2)')
% % % b2.FaceColor = '#bdbdbd';
% % % 
% % % subplot(2,1,1)
b1=bar(x,y);
title('Model')
ylabel('Predicted Cell Area (a.u.)')
b1.FaceColor = '#3182bd';

CA = responseset(1,1);
CA_01Dr = responseset(1,11);
CA_04Dr = responseset(1,21);
CA_16Dr = responseset(1,41);
Dr = 0.1593;

xx=categorical({namestimuli, strcat(namestimuli,'+low Midostaurin'),strcat(namestimuli,'+medium Midostaurin'),strcat(namestimuli,'+high Midostaurin'),'Midostaurin'});
xx= reordercats(xx,{namestimuli, strcat(namestimuli,'+low Midostaurin'),strcat(namestimuli,'+medium Midostaurin'),strcat(namestimuli,'+high Midostaurin'),'Midostaurin'});
yy=[CA,CA_01Dr,CA_04Dr,CA_16Dr,Dr];

figure
bb=bar(xx,yy);
bb.FaceColor = '#92c5de';
ylabel('Predicted Cell Area')