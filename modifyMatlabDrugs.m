function [status, result] = modifyMatlabDrugs(ODEfile,Paramfile)
% Modifies the existing ODE files (generated from NETFLUX) to allow drug
% simulations

% New lines to be inserted into the ODEfile
ins1 = ["%%%%% start new code from python";
    "w = rpar(1,:);";
    "n = rpar(2,:);";
    "EC50 = rpar(3,:);";
    "dose = rpar(4,:);";
    "drugBinding = rpar(5,:);";
    "drugAgonism = rpar(6,:);";
    "rpar = [w;n;EC50;dose;drugBinding;drugAgonism];";
    "%%%%% end new code from python"];
ins2 = ["%%%%% start new code from python";
    "function calcFact = calculateFact(EC50, n, w, x)";
    "beta = (EC50.^n - 1)./(2*EC50.^n - 1); ";
    "K = (beta - 1).^(1./n); ";
    "calcFact = w.*(beta.*x.^n)./(K.^n + x.^n);";
    "";
    "if calcFact > w ";
    "    calcFact = w; % cap fact(x)<= 1 ";
    "elseif calcFact > 1";
    "    calcFact = 1;";
    "elseif calcFact < 0";
    "    calcFact = 0;";
    "end  ";
    "";
    "function fact = act(x,rpar) ;"
    "% hill activation function with parameters w (weight), n (Hill coeff), EC50 ";
    "    w = rpar(1); ";
    "    n = rpar(2); ";
    "    EC50 = rpar(3); ";
    "    dose = rpar(4);";
    "    drugBinding = rpar(5);";
    "    drugAgonism = rpar(6);";
    "";
    "    if drugBinding == 1";
    "        x = x - dose;";
    "    end";
    "    if x < 0";
    "        x = 0;";
    "    end";
    "    if x > 1";
    "        x = 1;";
    "    end";
    "";
    "    fact = calculateFact(EC50, n, w, x);";
    "";
    "    %if strcmpi(drugType, \'Non-Competitive Agonist\') ";
    "    if drugBinding == -1 && drugAgonism == 1";
    "        fact = fact * (1-dose) + dose;";
    "    %elseif strcmpi(drugType, \'Non-Competitive Antagonist\')";
    "    elseif drugBinding == -1 && drugAgonism == -1";
    "        fact = fact * (1-dose);";
    "    end";
    "";
    "function finhib = inhib(x,rpar) ";
    "% inverse hill function with parameters w (weight), n (Hill coeff), EC50 ";
    "    w = rpar(1); ";
    "    n = rpar(2); ";
    "    EC50 = rpar(3); ";
    "    dose = rpar(4);";
    "    drugBinding = rpar(5);";
    "    drugAgonism = rpar(6);";
    "";
    "%     if drugType == -1";
    "%         w = alteration .* w;";
    "%     elseif drugType == 1";
    "%         x = x + alteration; % THIS IS THE OPPOSITE (\'+\', \'-\') FOR INHIBITION RXNS";
    "%     end";
    "    %if strcmpi(drugType, \'Competitive Agonist\') || strcmpi(drugType, \'Competitive Antagonist\') ";
    "    if drugBinding == 1";
    "        x = x + dose;";
    "    end";
    "    if x < 0";
    "        x = 0;";
    "    end";
    "    if x > 1";
    "        x = 1;";
    "    end";
    "";
    "    fact = calculateFact(EC50, n, w, x);";
    "    finhib = w - fact;";
    "";
    "    %if strcmpi(drugType, \'Non-Competitive Agonist\') ";
    "    if drugBinding == -1 && drugAgonism == 1";
    "        finhib = finhib * (1-dose) + dose;";
    "    %elseif strcmpi(drugType, \'Non-Competitive Antagonist\')";
    "    elseif drugBinding == -1 && drugAgonism == -1";
    "        finhib = finhib * (1-dose);";
    "    end";
    "";
    "function z = OR(x,y)  ";
    "% OR logic gate ";
    "    z = x + y - x*y;";
    "";
    "function z = AND(rpar,varargin) ";
    "% AND logic gate, multiplying all of the reactants together ";
    "    w = rpar(1); ";
    "    if w == 0 ";
    "        z = 0; ";
    "    else ";
    "        v = cell2mat(varargin); ";
    "        z = prod(v)/w^(nargin-2);  ";
    "    end ";
    "%%%%% end new code from python"];

%read in ODEfile
ofile = readlines(ODEfile);
p1="[rpar,tau,ymax,speciesNames]=params{:}; ";
res1 = find(strcmp(p1,ofile),1);
p2="% utility functions ";
res2 = find(strcmp(p2,ofile),1);
%insert new lines
newofile = vertcat(ofile(1:res1),ins1,ofile(res1+1:res2),ins2);
newofile(1)="function dydt=tempDrugODE(t,y,params)";

% New lines to be inserted into ParamFile
ins3 = ["%%%%% start new code from python";
    "dose = ones(1, length(w));";
    "drugBinding = zeros(1, length(w));";
    "drugAgonism = zeros(1, length(w));";
    "rpar = [w;n;EC50;dose;drugBinding;drugAgonism];";
    "%%%%% end new code from python"];

%read in ParamFile
pfile = readlines(Paramfile);
p3="rpar = [w;n;EC50];";
res3 = find(strcmp(p3,pfile),1);

%insert new lines
newpfile = vertcat(pfile(1:res3),ins3,pfile(res3+1:end));
newpfile(1) = "function [params,y0] = tempDrugODE_params()";

%write to new functions
writelines(newofile,"tempDrugODE.m")
writelines(newpfile,"tempDrugODE_params.m")

status = "done"; result = [size(newofile);size(newpfile)];
end