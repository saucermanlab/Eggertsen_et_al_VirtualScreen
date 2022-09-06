% Written by Taylor Eggertsen, modified from Anirudha Chandrabhatla
% Last Updated 7/12/2022
% Version 4.0
%
% Always run the code simulateDrugs.m before running this
% code to ensure you are running the right network, input, and drug
% combination

clear all; close all; clc
%% Inputs for simulations

opts = detectImportOptions('DrugsToSimulate.csv');
opts = setvartype(opts,{'AgonistTargetIndex','AntagonistTargetIndex'},'char');
drugsToSimulate = readtable('DrugsToSimulate.csv',opts);

formattedReactions = load('formattedReactions.mat');
formattedReactions = formattedReactions.formattedReactions; % Extract from struct

% Set drug dose 
alteration_antag = [0.8]; 
alteration_ag = alteration_antag;
inputs = ["AngII","ANPi","BNPi","CT1","EGF","ET1","FGF","IGF1","IL6","ISO","LIF","NE","NRG1","PE","Stretch","TGF\beta","TNFa"];
drugs = ["Midostaurin"];%,"Acetylsalicylic acid","Adapalene","Alpelisib","Aprindine","Baricitinib","Bosutinib","Celecoxib","Cobimetinib","Dabrafenib","Ellagic acid","Midostaurin","Purvalanol","Regorafenib","Resveratrol","Zanubrutinib"];

for kk = 1%:length(drugs)
    inp = 16 ; % input for model
    [~,drugindex] = ismember(drugs(kk),drugsToSimulate.Drug);
    drugofinterest = drugindex+1; 
    phen = 19; % position of phenotype of interest (e.g. 19 = Cell Area)
    KD_thresh = 0.03; % threshold determines stringency of knockdown screen

    % Parameters and initial conditions
    [params,y0] = tempDrugODE_params;
    tspan = [0 50]; options = [];

    % Pull parameters out in order to alter them in the loop without negating the previous alteration
    [rpar,tau,ymax,speciesNames]=params{:}; 
    w = rpar(1,:);
    n = rpar(2,:);
    EC50 = rpar(3,:);
    dose = rpar(4,:);
    drugBinding = rpar(5,:); 
    drugAgonism = rpar(6,:); 
    % Input Stimulus
    w(inp) = 0.1;
    
    % Repack parameters
    rpar = [w;n;EC50;dose;drugBinding;drugAgonism];
    inputNodeW = num2cell(1:length(speciesNames)); % Nodes to test drug against
    params = {rpar,tau,ymax,speciesNames};
    [t,y] = ode15s(@tempDrugODE,tspan,y0,options,params);
    yEnd = y(end,:)';

    % Reset initial y values
    y0 = real(yEnd); 
 
    knockoutData_drug = table('Size', [(length(speciesNames)+1) 2], 'VariableTypes', {'string', 'double'});
    knockoutData_control = table('Size', [(length(speciesNames)+1) 2], 'VariableTypes', {'string', 'double'});

    %% Drug knockdowns
    % Iterate through network nodes and simulate the knockdown of each species after addition of each drug into the network.
    for i = drugofinterest-1 % Input the row number (in 'drugsToSimulate.mat') of the drug you want to analyze.

        disp(drugsToSimulate{i, 1}); % Prints drug name to the command window 
        % simulate knockdowns
        for r = 0:length(speciesNames)
            disp(num2str(r));
            ymax_Knockout = ymax;
            if r > 0
                ymax_Knockout(r) = 0.01;
            end
            for j = 1:length(alteration_ag)
                [t2,y2] = drugtreatment(dose,drugBinding,drugAgonism,drugsToSimulate,formattedReactions,alteration_antag,w,n,EC50,tau,ymax_Knockout,speciesNames,y0,i,j);
                yEnd1 = y2(end,:);

            end 
            % Calculate final cell area
            if r == 0
                T_Fin = y2(end,:)'; T0 = y2(1,:)';
                T_Diff = T_Fin-T0;
                ySim = real(y2(:,phen))';
                [c1] = ySim;
                c1End_control = c1(end);
                c1End_knockout = c1(end);
                label = 'Control';
            else
                ySim = real(y2(:,phen))';
                [c1] = ySim;
                c1End_knockout = c1(end);
                label = speciesNames{r};
            end
            
            % Calculates the change in cell area from the knockout (in the
            % presence of stimulus and drug) and control (also in presence
            % of stimulus and drug)
            dataVector = {label, c1End_knockout-c1End_control};
            knockoutData_drug(r+1, :) = dataVector;
        end

        X = categorical(knockoutData_drug{:,1}); 
        X = reordercats(X,string(X));

        figure
        bar(X,knockoutData_drug{:,2})
        % % % % % set(gca, 'xtick', 1:(length(speciesNames)+1));
        set(gca, 'xticklabel', knockoutData_drug{:,1},'FontSize', 14);
        xtickangle(90)
        ylabel(strcat('Change in ',{' '},char(speciesNames(phen)),' ((Stimulus+Drug+KD) - (Stimulus+Drug))'))
        title(strcat(drugsToSimulate{i, 1},'; ',char(inputs(inp))))
        xlabel('Knockdowns')
    end

    %% Control knockouts

    % Set drug dose or doses (as a vector)
    dose_antag = [0.8]; 
    dose_ag = dose_antag;

    % Parameters and initial conditions
    [params,y0] = tempDrugODE_params;
    tspan = [0 50]; options = [];
    % Pull parameters out in order to alter them in the loop without negating the previous alteration
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
    params = {rpar,tau,ymax,speciesNames};
    [t,y] = ode15s(@tempDrugODE,tspan,y0,options,params);
    yEnd = y(end,:)';

    % Reset initial y values
    y0 = real(yEnd); 
    
    % Perturbed simulations
    sens = zeros(height(drugsToSimulate),length(inputNodeW));

    drugTargetResponse = zeros(1, height(drugsToSimulate)); 
    drugTargetResponseID = zeros(1, height(drugsToSimulate));

    % simulate knockdowns
    for r = 0:length(speciesNames)
        disp(num2str(r))
        ymax_Knockout = ymax;
        if r > 0
            ymax_Knockout(r) = 0.01;
        end

        params = {rpar,tau,ymax_Knockout,speciesNames};
        tspan = [0 50]; 
        options = []; 
        [t2,y2] = ode15s(@tempDrugODE,tspan,y0,options,params); 
        yEnd1 = y2(end,:);

        % Calculate final cell areas 
        if r == 0
            ySim = real(y2(:,phen))';
            [c2] = ySim;
            c2End_control = c2(end);
            c2End_knockout = c2(end);
            label = 'Control';
        else
            ySim = real(y2(:,phen))';
            [c2] = ySim;
            c2End_knockout = c2(end);
            label = speciesNames{r};
        end

        % Calculates the change in cell area from the knockout (in the
        % presence of stimulus alone) and control (also in presence
        % of stimulus alone)
        dataVector_control = {label, c2End_knockout-c2End_control};
        knockoutData_control(r+1, :) = dataVector_control;
    end

    figure
    bar(X,knockoutData_control{:,2})
    % % % % % set(gca, 'xtick', 1:(length(speciesNames)+1));
    set(gca, 'xticklabel', knockoutData_drug{:,1},'FontSize', 14);
    xtickangle(90)
    ylabel(strcat('Change in ',{' '},char(speciesNames(phen)),' ((Stimulus+KD) - Stimulus)'))
    title(strcat('No Drug; ',char(inputs(inp))))
    xlabel('Knockdowns')

    % Calculate the 'difference of the difference' - This finds how the
    % changes in cell area in each analysis differ from each other,
    % resulting in a measurement of how each node contributes to the
    % mechanism of the drug
    knockoutData = knockoutData_drug{:,2}-knockoutData_control{:,2};
   
    %Filter by nodes that were affected by drug action
    for ii=1:length(T_Diff)
        if abs(T_Diff(ii)) < 0.01
            knockoutData(ii+1) = 0;
        end
    end
    
    %Filter by effect size (identifies the most important nodes)
    thresh = 0.1*max(abs(knockoutData));
    for ii = 1:length(knockoutData)
        if abs(knockoutData(ii)) < thresh
            knockoutData(ii) = 0;
        end
    end

    figure
    bar(X,knockoutData)
    % % % % % set(gca, 'xtick', 1:(length(speciesNames)+1));
    set(gca, 'xticklabel', knockoutData_drug{:,1},'FontSize', 14);
    xtickangle(90)
    ylabel(strcat('Difference of KD Effect on',{' '},char(speciesNames(phen))))
    title(strcat(drugsToSimulate{i, 1},' - No Drug; ',char(inputs(inp))))
    xlabel('Knockdowns')
    
    knockout(:,kk) = knockoutData/max(abs(knockoutData)); %check this
end

%write data to text file
T = table(X,knockoutData,'VariableNames',{'Species_name','KDd'}); 
T(1,:) = []; filename = strcat(char(drugsToSimulate{i, 1}),'_',char(inputs(inp)),'_KDdata.txt');
writetable(T, filename)
%write more data to text file
U = table(X,knockoutData_control{:,2},'VariableNames',{'Species_name','KD_cellarea'});
U(1,:) = []; filename = strcat(char(drugsToSimulate{i, 1}),'_',char(inputs(inp)),'_KD_cellarea');
writetable(U,filename)

%% Select only certain data points based off of combined data structure
knockoutData_controlslim = [0]; knockoutData_drugslim = [0]; X_slim = ['Control'];
for ii=1:length(knockoutData)
    if abs(knockoutData(ii))>KD_thresh
        knockoutData_controlslim = [knockoutData_controlslim;knockoutData_control{ii,2}];
        knockoutData_drugslim = [knockoutData_drugslim;knockoutData_drug{ii,2}];
        X_slim = [X_slim;knockoutData_drug{ii,1}];
    end
end

%plot drug and control graphs for Cell Area
figure
subplot(2,1,2)
bar(knockoutData_drugslim+c1End_control)
set(gca, 'xtick', 1:length(X_slim));
set(gca, 'xticklabel', X_slim,'FontSize', 10);
xtickangle(90)
ylabel(char(speciesNames(phen)))
xlabel('Drug + Knockdowns')
ylim([0 1]);
title(drugsToSimulate{i, 1})

%%%figure
subplot(2,1,1)
bar(knockoutData_controlslim+c2End_control)
set(gca, 'xtick', 1:length(X_slim));
set(gca, 'xticklabel', X_slim,'FontSize', 10);
xtickangle(90)
ylabel(char(speciesNames(phen)))
xlabel('Knockdowns')
ylim([0 1]);
title('No Drug')

%eliminate empty rows
specNames = [' ', speciesNames];
for col = 1:size(knockout,2)
    for ro = 1:size(knockout,1)
        if abs(knockout(ro,col)) < 10^-3
            knockout(ro,col) = 0;
        end
    end
end
clean = knockout;
knockout(~any(clean,2),:)=[];
specNames(~any(clean,2))=[];

figure;
set(gca, 'Visible', 'on');
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
% % %maxVal = 2;
% % caxis([-1, 1]);
imagesc(knockout,[-1,1]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:kk);
set(gca,'XTickLabel',drugs);
% % % xticklabel_rotate([], 45);
xlabel('Drugs');
set(gca,'YTick',1:length(specNames));
set(gca,'YTickLabel',specNames,'fontsize',10);
ylabel('Species');
title(strcat('Mechanistic Screen: ',inputs(inp)));
xtickangle(90)
colorbar;


