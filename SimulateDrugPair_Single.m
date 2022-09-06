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
% % % [status, result] = generateDrugODEs(model);

%% Create DrugsToSimulate.csv
% % % [status, result] = generateDrugsToSimulate(model);

%% Inputs for simulations
opts = detectImportOptions('DrugsToSimulate.csv');
opts = setvartype(opts,{'AgonistTargetIndex','AntagonistTargetIndex'},'char');
drugsToSimulate = readtable('DrugsToSimulate.csv',opts);

formattedReactions = load('formattedReactions.mat');
formattedReactions = formattedReactions.formattedReactions; % Extract from struct
% formattedReactions(1,:) = []; % remove empty (first) row

%% Set drug dose or doses (as a vector) & Sensitivity analysis


alteration_antag = [0.1] ;%%% [0,0.5,0,0.5,0,0.5]; % variable 0.5
inputs = ["AngII","ANPi","BNPi","CT1","EGF","ET1","FGF","IGF1","IL6","ISO","LIF","NE","NRG1","PE","Stretch","TGFB","TNFa"];
alteration_ag = alteration_antag;
[padding,list_numDrugs, list_drugAction, list_drugType, list_drugTarget, rowLabels_drugs] = createRowLabels(drugsToSimulate);
drug1 = 41; drug1=drug1-1;
drug2 = 27; drug2=drug2-1; %27 or 112

% Node parameters
[params,y0] = tempDrugODE_params;
[rpar,tau,ymax,speciesNames]=params{:};
w = rpar(1,:);
%New input value
for i = 1:length(inputs)
    w(i) = 0.03; 
end

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
for i = [drug1,drug2] % Iterate through drugs
    for j = [drug1,drug2]
        
        doseNew1 = dose; doseNew2 = dose;
        drugBindingNew1 = drugBinding; drugBindingNew2 = drugBinding;
        drugAgonismNew1 = drugAgonism; drugAgonismNew2 = drugAgonism;
        
        % drug1 is an agonist
        if strcmp(drugsToSimulate.IsAgonist{i}, 'Yes') == 1
            if isempty(find(drugsToSimulate.AgonistTarget{i} == ';', 1)) % Drug has one agonist target
                locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AgonistTarget{i});
                targetnode = drugsToSimulate.AgonistTarget{i};
                if strcmp('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, agonist
                    drugBindingNew1(locationOfReactions) = 1;
                    drugAgonismNew1(locationOfReactions) = 1;
                    doseNew1(locationOfReactions) = -1*alteration_antag(1);
                else % Non-Competitive, agonist
                    drugBindingNew1(locationOfReactions) = -1;
                    drugAgonismNew1(locationOfReactions) = 1;
                    doseNew1(locationOfReactions) = alteration_antag(1);
                end
            else % Drug has multiple targets
                geneIDsOfTargets = strsplit(drugsToSimulate.AgonistTarget{i}, ';');
                targetnode = geneIDsOfTargets{1};
                for m = 1:length(geneIDsOfTargets)
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{m});
                    if strcmp('Competitive', drugsToSimulate.DrugAction{i})
                        drugBindingNew1(locationOfReactions) = 1;
                        drugAgonismNew1(locationOfReactions) = 1;
                        doseNew1(locationOfReactions) = -1*alteration_antag(1);
                    else 
                        drugBindingNew1(locationOfReactions) = -1;
                        drugAgonismNew1(locationOfReactions) = 1;
                        doseNew1(locationOfReactions) = alteration_antag(1);
                    end
                end
            end
        end

        % drug1 is an antagonist
        if strcmp(drugsToSimulate.IsAntagonist{i}, 'Yes') == 1
            if isempty(find(drugsToSimulate.AntagonistTarget{i} == ';', 1)) % Drug has one antagonist target
                locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AntagonistTarget{i});
                targetnode = drugsToSimulate.AntagonistTarget{i};
                if strcmp('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, antagonist
                    drugBindingNew1(locationOfReactions) = 1;
                    drugAgonismNew1(locationOfReactions) = -1;
                    doseNew1(locationOfReactions) = alteration_antag(1);
                else
                    drugBindingNew1(locationOfReactions) = -1;
                    drugAgonismNew1(locationOfReactions) = -1;
                    doseNew1(locationOfReactions) = alteration_antag(1);
                end
            else
                geneIDsOfTargets = strsplit(drugsToSimulate.AntagonistTarget{i}, ';');
                targetnode = geneIDsOfTargets{2};
                for p = 1:length(geneIDsOfTargets)
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{p});
                    if strcmp('Competitive', drugsToSimulate.DrugAction{i})
                        drugBindingNew1(locationOfReactions) = 1;
                        drugAgonismNew1(locationOfReactions) = -1;
                        doseNew1(locationOfReactions) = alteration_antag(1);
                    else 
                        drugBindingNew1(locationOfReactions) = -1;
                        drugAgonismNew1(locationOfReactions) = -1;
                        doseNew1(locationOfReactions) = alteration_antag(1);
                    end
                end
            end
        end

        % drug2 is an agonist
        if strcmp(drugsToSimulate.IsAgonist{j}, 'Yes') == 1
            if isempty(find(drugsToSimulate.AgonistTarget{j} == ';', 1)) % Drug has one agonist target
                locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AgonistTarget{j});
                targetnode = drugsToSimulate.AgonistTarget{j};
                if strcmp('Competitive', drugsToSimulate.DrugAction{j}) % Competitive, agonist
                    drugBindingNew2(locationOfReactions) = 1;
                    drugAgonismNew2(locationOfReactions) = 1;
                    doseNew2(locationOfReactions) = -1*alteration_antag(1);
                else % Non-Competitive, agonist
                    drugBindingNew2(locationOfReactions) = -1;
                    drugAgonismNew2(locationOfReactions) = 1;
                    doseNew2(locationOfReactions) = alteration_antag(1);
                end
            else % Drug has multiple targets
                geneIDsOfTargets = strsplit(drugsToSimulate.AgonistTarget{j}, ';');
                targetnode = geneIDsOfTargets{1};
                for m = 1:length(geneIDsOfTargets)
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{m});
                    if strcmp('Competitive', drugsToSimulate.DrugAction{j})
                        drugBindingNew2(locationOfReactions) = 1;
                        drugAgonismNew2(locationOfReactions) = 1;
                        doseNew2(locationOfReactions) = -1*alteration_antag(1);
                    else 
                        drugBindingNew2(locationOfReactions) = -1;
                        drugAgonismNew2(locationOfReactions) = 1;
                        doseNew2(locationOfReactions) = alteration_antag(1);
                    end
                end
            end
        end

        % drug2 is an antagonist
        if strcmp(drugsToSimulate.IsAntagonist{j}, 'Yes') == 1
            if isempty(find(drugsToSimulate.AntagonistTarget{j} == ';', 1)) % Drug has one antagonist target
                locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AntagonistTarget{j});
                targetnode = drugsToSimulate.AntagonistTarget{j};
                if strcmp('Competitive', drugsToSimulate.DrugAction{j}) % Competitive, antagonist
                    drugBindingNew2(locationOfReactions) = 1;
                    drugAgonismNew2(locationOfReactions) = -1;
                    doseNew2(locationOfReactions) = alteration_antag(1);
                else
                    drugBindingNew2(locationOfReactions) = -1;
                    drugAgonismNew2(locationOfReactions) = -1;
                    doseNew2(locationOfReactions) = alteration_antag(1);
                end
            else
                geneIDsOfTargets = strsplit(drugsToSimulate.AntagonistTarget{j}, ';');
                targetnode = geneIDsOfTargets{2};
                for p = 1:length(geneIDsOfTargets)
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{p});
                    if strcmp('Competitive', drugsToSimulate.DrugAction{j})
                        drugBindingNew2(locationOfReactions) = 1;
                        drugAgonismNew2(locationOfReactions) = -1;
                        doseNew2(locationOfReactions) = alteration_antag(1);
                    else 
                        drugBindingNew2(locationOfReactions) = -1;
                        drugAgonismNew2(locationOfReactions) = -1;
                        doseNew2(locationOfReactions) = alteration_antag(1);
                    end
                end
            end
        end
        if i==j
            doseNew = doseNew1;
        else
            doseNew = doseNew1 + doseNew2;
        end
        drugBindingNew = drugBindingNew1 + drugBindingNew2;
        drugAgonismNew = drugAgonismNew1 + drugAgonismNew2;
        for ii = 1:length(doseNew)
            if drugBindingNew(ii) < -1
                drugBindingNew(ii) = -1;
            elseif drugBindingNew(ii) > 1
                drugBindingNew(ii) = 1;
            end
            if drugAgonismNew(ii) < -1
                drugAgonismNew(ii) = -1;
            elseif drugAgonismNew(ii) > 1
                drugAgonismNew(ii) = 1;
            end
        end
        
        rparNew = [w;n;EC50;doseNew;drugBindingNew;drugAgonismNew];
        paramsNew = {rparNew,tau,ymax,speciesNames};
        tspan = [0 50]; options = []; 
        [t2,y2] = ode15s(@tempDrugODE,tspan,y0,options,paramsNew); 
               
        ySimEnd = y2(end,:)';    
        drugpair(i,j) = ySimEnd(19);
    end   
end
initial = 0.159274221840270;
nodrug = y0(19);
outdrug1 = drugpair(drug1,drug1);
outdrug2 = drugpair(drug2,drug2);
outdrugpair = drugpair(drug1,drug2);

%% Figures

x=categorical({'Control', 'No Drug', drugsToSimulate.Drug{drug1},drugsToSimulate.Drug{drug2},strcat(drugsToSimulate.Drug{drug1},' + ',drugsToSimulate.Drug{drug2})});
x= reordercats(x,{'Control', 'No Drug', drugsToSimulate.Drug{drug1},drugsToSimulate.Drug{drug2},strcat(drugsToSimulate.Drug{drug1},' + ',drugsToSimulate.Drug{drug2})});
y=[initial,nodrug,outdrug1,outdrug2,outdrugpair];
figure
b=bar(x,y);
b.FaceColor = '#3182bd';
ylim([0,1])
ylabel('Cell Area')
xtickangle(0)

