% Written by Taylor Eggertsen, modified from Anirudha Chandrabhatla
% Last Updated 7/12/2022
% Version 6.0
%
% Performs a virtual drug screen on DrugBank list -
% outputs a heatmap of drug predictions in the given network model

%% If any inputs have changed, run the following section: 
close all; clear all; clc;

% Input the network mode you are testing
model='Hypertrophy_Model.xlsx';

warning off;
formattedReactions = table;
% Species/Reaction information from network, 'species'/'reactions' tab 
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

%% Generate the drug modified  ODE and parameter files from Netflux
[status, result] = generateDrugODEs(model);

%% Create DrugsToSimulate.csv which contains drugs and targets from DrugBank
% % % [status, result] = generateDrugsToSimulate(model); % can be commented out once drug file is obtained

%% Inputs for simulations - uses the drug file created above
opts = detectImportOptions('DrugsToSimulate.csv');
opts = setvartype(opts,{'AgonistTargetIndex','AntagonistTargetIndex'},'char');
drugsToSimulate = readtable('DrugsToSimulate_52.csv',opts);

formattedReactions = load('formattedReactions.mat');
formattedReactions = formattedReactions.formattedReactions; % Extract from struct
% % % formattedReactions(1,:) = []; % remove empty (first) row
% input 1 for drug predictions; input 0 to directy see target activity
sim_output = 0; 

%% Run the Virtual Screen

% Set drug dose or doses (as a vector) & specify model inputs
if sim_output == 1
    alteration_antag = [0.8] ;
    inputs = ["Control","AngII","ANPi","BNPi","CT1","EGF","ET1","FGF","IGF1","IL6","ISO","LIF","NE","NRG1","PE","Stretch","TGFB","TNFa"];
elseif sim_output == 0
    alteration_antag = [0:0.05:1];
    inputs = ["NE"];
end
alteration_ag = alteration_antag;
% % % [padding,list_numDrugs, list_drugAction, list_drugType, list_drugTarget, rowLabels_drugs] = createRowLabels(drugsToSimulate);

for inp = 1:length(inputs)
    disp(['Input # ', num2str(inp), ' of ',num2str(length(inputs))]) 
    % Unpack node parameters
    [params,y0] = tempDrugODE_params;
    [rpar,tau,ymax,speciesNames]=params{:};
    w = rpar(1,:);
    n = rpar(2,:);
    EC50 = rpar(3,:);
    dose = rpar(4,:);
    drugBinding = rpar(5,:); 
    drugAgonism = rpar(6,:); 
    %New input value
    if inp>1
        w(inp-1) = 0.1; 
    elseif sim_output == 0
        w(12) = 0.1;
    end
    
    % Repack node parameters
    rpar = [w;n;EC50;dose;drugBinding;drugAgonism];
    inputNodeW = num2cell(1:length(speciesNames)); % Nodes to test drug against
    params=[rpar,tau,ymax,speciesNames];
    % Steady-state Control Simulation
    tspan = [0 50]; options = [];
    [t,y] = ode15s(@tempDrugODE,tspan,y0,options,params);
    yEnd = y(end,:)';

    % Reset initial y values
    y0 = real(yEnd); 
    initArea(inp)=y0(19);

    sens = zeros(height(drugsToSimulate),length(inputNodeW));
    for i = 1:height(drugsToSimulate) % Iterate through drugs
        for j = 1:length(alteration_antag) 
            [t2,y2,targetnode] = drugtreatment(dose,drugBinding,drugAgonism,drugsToSimulate,formattedReactions,alteration_antag,w,n,EC50,tau,ymax,speciesNames,y0,i,j);
            ySimEnd = y2(end,:)';

            for k = 1:length(inputNodeW)
                sens(i,k) = real(ySimEnd(k)) - real(yEnd(k)); %expressing sensitivity as fold change
            end

            % Sensitivity Analysis
            drugs = string(drugsToSimulate{:,1});
            agonisttargets = string(drugsToSimulate{:,4});
            antagonisttargets = string(drugsToSimulate{:,8});
            specs = speciesNames;
            drugResponse = real(sens);

% % %             rowLabels=rowLabels_drugs; 
            alltargets = agonisttargets;

% % %             rowLabels_clustered = [];
% % %             for m = 1:length(rowLabels)
% % %                 rowLabels_clustered = vertcat(rowLabels_clustered, rowLabels(m));
% % %             end


            %whittling drugResponse to outputs only (hypertrophy model outputs)
            if sim_output == 1
                outspecs = ["CellArea","aMHC","bMHC","ANP","BNP","sACT","SERCA"];   
            elseif sim_output == 0
                outspecs = speciesNames;
            end
            drugResponse2 = zeros(size(drugResponse,1),length(outspecs)); 
            specs2 = cell(1,length(outspecs));
            for jj = 1:length(outspecs)
                outcol = find(strcmp(specs,outspecs(jj)));
                drugResponse2(:,jj) = drugResponse(:,outcol);
                specs2(jj) = specs(outcol);
            end

            % Dose Response
            if sim_output == 1
                dpos = 1; 
                doseresponse(i,j) = drugResponse2(i,dpos);
            elseif sim_output == 0
                if class(targetnode)=='cell'
                    for kk = 1:length(targetnode)
                        dpos(kk) = find(strcmp(targetnode(kk),outspecs));
                    end
                    if sum(drugResponse2(i,dpos)) < 0
                        doseresponse(i,j) = min(drugResponse2(i,dpos));
                    else
                        doseresponse(i,j) = max(drugResponse2(i,dpos));
                    end
                elseif class(targetnode)=='char'
                    dpos = find(strcmp(targetnode,outspecs));
                    doseresponse(i,j) = drugResponse2(i,dpos);
                end
            end
        end   
    end

    % Group drug responses by mechanism of action
    for rep = 1:length(alltargets)
        if strcmp(alltargets(rep),"")==1
            alltargets(rep)=antagonisttargets(rep);
        end
    end
    aglist=[""];antaglist=[""];count=1;
    for ii = 1:length(alltargets)
        if agonisttargets(ii)==""
            if ismember(alltargets(ii),antaglist)==0
                antaglist=[antaglist,alltargets(ii)];
                DG(count) = drugs(ii);
                DR(count,:) = doseresponse(ii,:);
                Dg2(count,:) = drugResponse2(ii,:);
% % %                 RL(count) = rowLabels(ii);
                count = count+1;
            end
        elseif antagonisttargets(ii)==""
            if ismember(alltargets(ii),aglist)==0
                aglist=[aglist,alltargets(ii)];
                DG(count) = drugs(ii);
                DR(count,:) = doseresponse(ii,:);
                Dg2(count,:) = drugResponse2(ii,:);
% % %                 RL(count) = rowLabels(ii);
                count = count+1;
            end
        end
    end
    doseresponse = DR;
    drugs = DG;
    drugResponse2 = Dg2;
% % %     rowLabels = RL';
    
    %drugs by input
    if sim_output == 1
        drugbyinput(:,inp) = doseresponse;
    elseif sim_output == 0
        drugbyinput(:,inp) = doseresponse(:,end);
    end
end

%% Figures

figure
x = categorical(inputs);
x = reordercats(x,inputs);
bar(x,initArea);
xlabel("Biochemical inputs");
ylabel("Predicted Cell Area (a.u.)");

%clean up drugbyinput 
if sim_output == 1
    clean = drugbyinput;
elseif sim_output == 0
    clean = doseresponse;
end
for col = 1:size(clean,2)
    for ro = 1:size(clean,1)
        if abs(clean(ro,col)) < 10^-3
            clean(ro,col) = 0;
        end
    end
end
drugbyinput(~any(clean,2),:)=[];
doseresponse(~any(clean,2),:)=[];
drugResponse2(~any(clean,2),:)=[];
% % % rowLabels(~any(clean,2),:)=[];
drugs(~any(clean,2))=[]; 

figure;
set(gca, 'Visible', 'on');
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
% % %maxVal = 2;
% % caxis([-1, 1]);
if sim_output == 1
    imagesc(drugbyinput,[-1,1]);
    set(gca,'XAxisLocation','bottom');
    set(gca,'XTick',1:inp);
    set(gca,'XTickLabel',inputs,'fontsize',10);
    xlabel('Biochemical Inputs','fontsize',20);
elseif sim_output == 0
    imagesc(doseresponse,[-1,1]);
    set(gca,'XAxisLocation','bottom');
    set(gca,'XTick',1:length(alteration_antag));
    set(gca,'XTickLabel',alteration_antag,'fontsize',10);
    xlabel('Drug Dose','fontsize',20);
end
% % % xticklabel_rotate([], 45);
set(gca,'YTick',1:length(drugs));
set(gca,'YTickLabel',drugs,'fontsize',10);
ylabel('Representative Drugs','fontsize',16);
% % % title([strcat('Change in ',{' '},specs2{dpos})]);
xtickangle(90)
hcb=colorbar;
% % % % % % set(get(hcb,'label'),'string',strcat('Change in ',{' '},specs2{dpos}));
if sim_output == 1
    set(get(hcb,'label'),'string','Change in Cell Area');
elseif sim_output == 0
    set(get(hcb,'label'),'string','Change in Target Node Activity');
end
set(get(hcb,'label'),'fontsize',16); set(get(hcb,'label'),'rotation',90);
% % % colorbar;

%Cluster the above heatmap
if sim_output == 1
    cgo = clustergram(drugbyinput,'Colormap',myrgbcmap,'RowLabels',drugs,'ColumnLabels',inputs);
end

figure;
set(gca, 'Visible', 'on');
cmaprange = [0.5:0.005:1];
blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(myrgbcmap);
% % %maxVal = 2;
% % caxis([-1, 1]);
imagesc(drugResponse2,[-1,1]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(outspecs));
set(gca,'XTickLabel',outspecs,'fontsize',10);
xlabel('Phenotypic Outputs','fontsize',20);
% % % xticklabel_rotate([], 45);
set(gca,'YTick',1:length(drugs));
set(gca,'YTickLabel',drugs,'fontsize',10);
ylabel('Representative Drugs','fontsize',16);
% % % title([strcat('Change in ',{' '},specs2{dpos})]);
xtickangle(90)
hcb=colorbar;
% % % % % % set(get(hcb,'label'),'string',strcat('Change in ',{' '},specs2{dpos}));
set(get(hcb,'label'),'string','Change in Activity');
set(get(hcb,'label'),'fontsize',16); set(get(hcb,'label'),'rotation',90);
% % % colorbar;

% % % figure;
% % % set(gca, 'Visible', 'on');
% % % cmaprange = [0.5:0.005:1];
% % % blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
% % % myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
% % % colormap(myrgbcmap);
% % % % % %maxVal = 2;
% % % % % caxis([-1, 1]);
% % % imagesc(drugbyinput,[-1,1]);
% % % set(gca,'XAxisLocation','bottom');
% % % set(gca,'XTick',1:inp);
% % % set(gca,'XTickLabel',inputs,'fontsize',10);
% % % % % % xticklabel_rotate([], 45);
% % % xlabel('Biochemical Inputs');
% % % set(gca,'YTick',1:length(rowLabels));
% % % set(gca,'YTickLabel',rowLabels,'fontsize',10);
% % % set(gca,'fontname','FixedWidth');
% % % ylabel('Drug Targets');
% % % title([strcat('Change in Node Activity: ',specs2{dpos})]);
% % % xtickangle(90)
% % % colorbar;
% % % ylabel(colorbar,{'Change in Activity', '(drugged - non drugged)'}, 'Rotation', 270, 'FontWeight', 'bold', 'FontSize', 10);

%% Additional stats

for ii = 1:length(drugs)
    frac(ii) = 0;
    for jj = 1:length(inputs)
        if drugbyinput(ii,jj)<-0.05
            frac(ii) = frac(ii) + 1;
        end
    end
end
Drug = drugs'; Inhibited = frac';
efficacy = table(Drug,Inhibited); %Shows how many hypertrophic conditions are inhibited by drug
xxx = categorical(efficacy{:,1});
xxx = reordercats(xxx,efficacy{:,1});
figure
bar(xxx,efficacy{:,2})

robust = efficacy{:,2}>12;
robustDrugs = efficacy{robust,1};
noeffect = efficacy{:,2}==0;
noeffectDrugs = efficacy{noeffect,1};
context = logical(ones(length(efficacy{:,2}),1) - noeffect - robust);
contextDrugs = efficacy{context,1};