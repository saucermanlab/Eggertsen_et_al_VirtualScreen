function [status, result] = generateDrugODEs(model)
% Generates the ODE and parameter files from NETFLUX for the given model,
% and modifies them to allow drug simulation

if exist([pwd '\ODEfun_fromNetflux.m'],'file') == 2
    delete('ODEfun_fromNetflux.m');
end
if exist([pwd '\ODEfun_loadParams_fromNetflux.m'],'file') == 2
    delete('ODEfun_loadParams_fromNetflux.m');
end
namepos = findstr('.xls', model); namestr = cellstr(model(1:namepos-1));

[specID,reactionIDs,~,paramList,ODElist,~, error] = util.xls2Netflux(namestr,model);
commandODE = util.exportODE2(specID,paramList,ODElist);
[a,commandPARAM,b] = util.exportODE(specID,paramList,ODElist,'ODEfun');
util.textwrite('ODEfun_fromNetflux.m',commandODE);
util.textwrite('ODEfun_loadParams_fromNetflux.m',commandPARAM);

% Generate the drug-modified ODE and parameter files
ODEfile='ODEfun_fromNetflux.m';
ParamFile='ODEfun_loadParams_fromNetflux.m';
[status, result] = modifyMatlabDrugs(ODEfile,ParamFile);
end

