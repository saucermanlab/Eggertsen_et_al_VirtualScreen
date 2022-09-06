%% Run Validation Script + Excel Table
close all; clear; clc;

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

%% Generate the ODE and parameter files from NETFLUX
[status, result] = generateDrugODEs(model);

%% Validate the drug model
threshold = [0.001:0.0005:0.05]; %Can be performed for range of thresholds
for jj = 1:length(threshold)
    opts = detectImportOptions('DrugsToSimulate.csv');
    opts = setvartype(opts,{'AgonistTargetIndex','AntagonistTargetIndex'},'char');
    drugsToSimulate = readtable('DrugsToSimulate.csv',opts);

    formattedReactions = load('formattedReactions.mat');
    formattedReactions = formattedReactions.formattedReactions; % Extract from struct

    validationfname = 'HypertrophyPharm_Validation.xlsx'; % Experimental sheet to validate against

    [~, txt, raw] = xlsread(validationfname);
    validID = txt(2:end,1);
    dataName = [];

    % read the validation sheet
    [~, txt, raw] = xlsread(validationfname);
    % remove rows without data
    noData = cellfun(@(x)isequal(x,'No Data'), txt(1:end, 7));
    txt(noData, :) = []; 
    noData = cellfun(@isempty, txt(1:end, 7));
    txt(noData, :) = [];
    assignin('base', 'txt', txt);
    input1 = txt(2:end, 2);
    input2 = txt(2:end, 3);
    inputCode = txt(2:end, 4);
    measurement = txt(2:end,7); 
    outputSpec = txt(2:end, 5);
    control = txt(2:end, 6);
    validationIDs = txt(2:end, 1);
    validationTags = txt(2:end, 8);
    celltype = txt(2:end, 12);

    % set validation threshold change
    thresh = threshold(jj); % threshold, Ryall et al., 2012 set to 0.001 for sensitivity analysis
    inc = {'Increase'};
    dec = {'Decrease'};
    noc = {'No Change'};

    % determination of number of predictions matching references
    numMatching = 0; % number of predictions consistent with the qualitative literature species behavior

    % loop over all validation simulations read in from the excel sheet
    for i = 1:length(inputCode)
        disp(['Simulation # ', num2str(i), ' of ',num2str(length(inputCode))]) % write the simulation number to the command line to track loop progress

        % evaluate validation conditions from excel sheet
        eval(inputCode{i});
        alteration_ag = 0.8; alteration_antag = 0.8;

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

        for ii = drugofinterest-1 %input drug of interest
            for j = 1:length(alteration_antag) 
                [t2,y2] = drugtreatment(dose,drugBinding,drugAgonism,drugsToSimulate,formattedReactions,alteration_antag,w,n,EC50,tau,ymax,speciesNames,y0,ii,j);
                ySimEnd = y2(end,:)';
            end
        end

        % subtract final species activation from initial to determine the effect
        % the validation simulation has had on the species' activation with
        % respect to the defined threshold, then write out the qualitative
        % change of the species' activation state to the excel file
        yStartL{i} = yEnd;
        yEndL{i} = ySimEnd;

        % find indices of output species
        outputSpeciesIndex = zeros(1, length(measurement));
        for k = 1:length(outputSpec)
            [~,outputSpeciesIndex(k)] = ismember(outputSpec{k},speciesNames);
        end

        activityChange = real(yEndL{i}(outputSpeciesIndex(i)))-real(yStartL{i}(outputSpeciesIndex(i)));   

        if activityChange > thresh % increase
            prediction{i} = 'Increase';
            predChange{i} = num2str(activityChange);
            if isequal(inc,measurement(i))
                numMatching = numMatching + 1;
                match(i) = 1; %if the simulation matches the experimental validation put a 1 in the vector
            else
                match(i) = 0; %if the simulation does not match put a 0 in the matrix
            end
        elseif activityChange < -thresh % decrease
            prediction{i} = 'Decrease';
            predChange{i} = num2str(activityChange);
            if isequal(dec,measurement(i))
                numMatching = numMatching + 1;
                match(i) = 1;
            else
                match(i) = 0;
            end
        else % no change
            prediction{i} = 'No Change';
            predChange{i} = num2str(activityChange);
            if isequal(noc,measurement(i))
                numMatching = numMatching + 1;
                match(i) = 1;
            else
                match(i) = 0;
            end
        end
    end

    for j = 1:length(match)
        if match(j) == 1
            match2{j} = 'yes';
            match3(j) = 1;
        else
            match2{j} = 'no';
            match3(j) = 0;
        end
    end
    match = match2;

    BMatch = match3';
    ind = find(strcmp('In-Int', validationTags));
    byClass.inint = sum(BMatch(ind))/length(BMatch(ind))*100;

    ind = find(strcmp('comb-out', validationTags));
    byClass.comb = sum(BMatch(ind))/length(BMatch(ind))*100;

    ind = find(strcmp('in-out', validationTags));
    byClass.inout = sum(BMatch(ind))/length(BMatch(ind))*100;

    ind = find(strcmp('KO Int', validationTags));
    byClass.KO = sum(BMatch(ind))/length(BMatch(ind))*100;

    % output the KO'd rxn indices and the percent agreement
    resultChart = {validationIDs, input1, input2, outputSpec, control, measurement, prediction', predChange', match', validationTags, celltype}; %create a cell array showing input, output, and whether they matched validation
    header = {'ID', 'input' , 'input2',  'output', 'control', 'measurement', 'prediction', 'predicted change', 'match', 'tag', 'celltype'};
    resultChart = horzcat(resultChart{:});
    resultChart = vertcat(header, resultChart);
    delete(dataName);
    %xlswrite(dataName, resultChart);
    %resultChart;
    %csvwrite(dataName, resultChart);

    disp(['wrote ', which('Validation Results.xlsx')]);
    percentMatch = numMatching/length(measurement)*100
    z(jj) = percentMatch;
    match = []; match2 = []; match3 = [];
end
figure;
scatter(threshold,z);
ylabel('Validation Percent')
xlabel('Activity Threshold')

% % % ExpQual = [resultChart(2:end,6)];
% % % ModelQual = [resultChart(2:end,7)];
% % % ValQual = [resultChart(2:end,9)];
% % % CellType = [resultChart(2:end,11)];
% % % 
% % % YesInd = strmatch('yes',ValQual);
% % % 
% % % ExpQualInc = strmatch('Increase',ExpQual);
% % % ExpQualDec = strmatch('Decrease',ExpQual);
% % % 
% % % MQualInc = strmatch('Increase',ModelQual);
% % % MQualDec = strmatch('Decrease',ModelQual);
% % % 
% % % ExpVal = zeros(1, numel(ExpQual));
% % % ExpVal(ExpQualInc) = 1;
% % % ExpVal(ExpQualDec) = -1; 
% % % 
% % % MVal = zeros(1, numel(ModelQual));
% % % MVal(MQualInc) = 1;
% % % MVal(MQualDec) = -1; 
% % % 
% % % % Write validation results to the Validation spreadsheet and open
% % % A = [MVal' ExpVal'];
% % % xlswrite(validationfname,ModelQual,1,'F2')
% % % winopen(validationfname)

