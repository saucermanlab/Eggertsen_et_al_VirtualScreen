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

%% Generate the drug modified  ODE and parameter files from NETFLUX
[status, result] = generateDrugODEs(model);

%% Create DrugsToSimulate.csv
% % % [status, result] = generateDrugsToSimulate(model);

%% Inputs for simulations
opts = detectImportOptions('DrugsToSimulate.csv');
opts = setvartype(opts,{'AgonistTargetIndex','AntagonistTargetIndex'},'char');
drugsToSimulate = readtable('DrugsToSimulate_Polypharm.csv',opts);

formattedReactions = load('formattedReactions.mat');
formattedReactions = formattedReactions.formattedReactions; % Extract from struct
% formattedReactions(1,:) = []; % remove empty (first) row

%% Set drug dose or doses (as a vector) & Sensitivity analysis

alteration_antag = [0:0.02:1] ;%%% [0,0.5,0,0.5,0,0.5]; % variable 0.5
[padding,list_numDrugs, list_drugAction, list_drugType, list_drugTarget, rowLabels_drugs] = createRowLabels(drugsToSimulate);
poly = []; poly_label = [];

inputs = ["AngII","ANPi","BNPi","CT1","EGF","ET1","FGF","IGF1","IL6","ISO","LIF","NE","NRG1","PE","Stretch","TGFB","TNFa"];
for inp = 1:17 %input to network .
    poly_data = []; poly_name_ag = {}; poly_name_antag = {};
    inputs(inp)
    for drugofinterest = 22:28

        % Node parameters
        [params,y0] = tempDrugODE_params;
        [rpar,tau,ymax,speciesNames]=params{:};
        w = rpar(1,:);
        %New input value
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

        % Reset initial y values
        y0 = real(yEnd); 

        sens = zeros(height(drugsToSimulate),length(inputNodeW));
        for i = drugofinterest-1 %input drug of interest
            for j = 1:length(alteration_antag) 
                
                [t2,y2] = drugtreatment(dose,drugBinding,drugAgonism,drugsToSimulate,formattedReactions,alteration_antag,w,n,EC50,tau,ymax,speciesNames,y0,i,j);
                ySimEnd = y2(end,:)';

                count = 1; leg = [];
                specset = [19 11 12];
                for k = specset 
                     responseset(count,j) = y2(end,k);
                     leg = [leg, speciesNames(k)];
                     count = count + 1;
                end

            end
        end
        % Calculates the change in cell area (where greater inhibition is
        % more positive)
        poly_data = [poly_data, (y2(1,19)-y2(end,19))];
        poly_name_ag{end+1} = drugsToSimulate.AgonistTarget{i};
        poly_name_antag{end+1} = drugsToSimulate.AntagonistTarget{i};
    end
    
    poly = vertcat(poly,poly_data);
    poly_label = vertcat(poly_label,inputs(inp));
end

figure;
set(gca, 'Visible', 'on');
colormap default;
% % %maxVal = 2;
% % caxis([-1, 1]);
imagesc(poly,[0,1]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(poly_name_antag));
set(gca,'XTickLabel',poly_name_antag);
% % % xticklabel_rotate([], 45);
xlabel('Target Node(s) of Zanubrutinib');
set(gca,'YTick',1:length(poly_label));
set(gca,'YTickLabel',poly_label,'fontsize',10);
ylabel('Biochemical Inputs');
title(strcat(drugsToSimulate.Drug{i},' Polypharmacology'));
xtickangle(90)
hcb=colorbar;
set(get(hcb,'label'),'string','Hypertrophy Inhibition');