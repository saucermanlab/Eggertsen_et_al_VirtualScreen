function [status, result] = generateDrugsToSimulate(model)
% Uses a webscraping python code to extract drug data for list of DrugBank
% drugs

% Match full target name and drug ID for compatability with drug action webscraper
allPharm='"allPharm.csv"';
model = ['"' model '"'];
systemCommand = ['py drugID_targetName_match.py ',allPharm,' ',model];
[status, result] = system(systemCommand);

% Get Drug Actions, DrugActionsTest_Output.csv 
drugData='"drugID_names_matched.csv"'; %need to change this so it's just the drug ID's
systemCommand = ['py webScrapeDrugAction.py ',drugData];
[status, result] = system(systemCommand);

% Create DrugsToSimulate.csv
actionCSV='"DrugActionsTest_Output.csv"';
targetCSV='"allPharm.csv"';
drugNames='"drugIDs.csv"';
%call function, drugsToSimulate should then appear in working directory
systemCommand = ['py DrugID_Pull_ANTE2.py ', targetCSV,' ', actionCSV,' ', model,' ', drugNames];
[status, result] = system(systemCommand);
end

