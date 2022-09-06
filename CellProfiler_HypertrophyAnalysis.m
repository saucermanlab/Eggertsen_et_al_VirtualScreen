%Cell Profiler - Hypertrophy Analysis
% Taylor Eggertsen
%Initial version - 9/30/2021

close all; clear;clc;
%import area data
myocytes = readtable('hypertrophyCellProfilerDatamyocytes_ident.csv');
image = myocytes.ImageNumber;
area = myocytes.AreaShape_Area;

%Filter out the 'cells' that have area suspiciously equivalent to nuclei
nuclei = readtable('hypertrophyCellProfilerDatanucMyo_ident.csv');
nucarea = nuclei.AreaShape_Area;
for i = 1:height(myocytes)
    if area(i)>nucarea(i)
        myocytes.AreaShape_Area(i) = area(i);
    else
        myocytes.AreaShape_Area(i) = 0;
    end
end

%Define conditions
control = [1,10,51,60];
controlname = ['Control'];
c1 = [2,3,4,5]; c2 = [12,13,14,15]; c3 = [22,23,24,25]; 
c4 = [6,7,8,9]; c5 = [16,17,18,19]; c6 = [26,27,28,29]; 
c7 = [32,33,34,35]; c8 = [42,43,44,45]; c9 = [52,53,54,55]; 
c10 = [36,37,38,39]; c11 = [46,47,48,49]; c12 = [56,57,58,59];
c13 = [11,21,31,41]; c14 = [20,30,40,50];
c1name = ['320 nM in AngII']; c2name = ['80 nM in AngII']; c3name = ['20 nM in AngII']; 
c4name = ['320 nM in NE']; c5name = ['80 nM in NE']; c6name = ['20 nM in NE'];
c7name = ['320 nM in TGFb']; c8name = ['80 nM in TGFb']; c9name = ['20 nM in TGFb'];
c10name = ['AngII']; c11name = ['NE']; c12name = ['TGFb'];
c13name = ['Serum']; c14name = ['320 nM Midostaurin'];

% % % % % %% Find Mean Area & SEM
% % % % % 
% % % % % %control
% % % % % control_cells=[];
% % % % % for i=control
% % % % %     for j=1:length(image)
% % % % %         if image(j)==i
% % % % %             control_cells = [control_cells,area(j)];
% % % % %         end
% % % % %     end
% % % % % end
% % % % % controlmean = mean(control_cells);
% % % % % controlsem = std(control_cells)/sqrt(length(control_cells));
% % % % % %conditions
% % % % % c1_cells = []; c2_cells = []; c3_cells = []; 
% % % % % c4_cells = []; c5_cells = []; c6_cells = [];
% % % % % c7_cells = []; c8_cells = []; c9_cells = [];
% % % % % c10_cells = []; c11_cells = []; c12_cells = [];
% % % % % 

%% Alternate find Mean Area & SEM

%Break table into 60 subtables
G = findgroups(myocytes{:, 1});
T_split = splitapply( @(varargin) varargin, myocytes , G);

%Average each well's cell area
for i = 1:60
    marea(i) = mean(nonzeros(T_split{i,3}));
    count(i) = length(T_split{i,3});
% % %     figure
% % %     histogram(T_split{i,3},20)
end
control_cells=[]; control_count=[];
c1_cells = []; c2_cells = []; c3_cells = []; 
c4_cells = []; c5_cells = []; c6_cells = [];
c7_cells = []; c8_cells = []; c9_cells = [];
c10_cells = []; c11_cells = []; c12_cells = [];
c13_cells = []; c14_cells = [];
c1_count = []; c2_count = []; c3_count = []; 
c4_count = []; c5_count = []; c6_count = [];
c7_count = []; c8_count = []; c9_count = [];
c10_count = []; c11_count = []; c12_count = [];
c13_count = []; c14_count = [];
for i = 1:60
    if ismember(i,control)==1
        control_cells = [control_cells,marea(i)];
        control_count = [control_count,count(i)];
    elseif ismember(i,c1)==1
        c1_cells = [c1_cells,marea(i)];
        c1_count = [c1_count,count(i)];
    elseif ismember(i,c2)==1
        c2_cells = [c2_cells,marea(i)];
        c2_count = [c2_count,count(i)];
    elseif ismember(i,c3)==1
        c3_cells = [c3_cells,marea(i)];
        c3_count = [c3_count,count(i)];
    elseif ismember(i,c4)==1
        c4_cells = [c4_cells,marea(i)];
        c4_count = [c4_count,count(i)];
    elseif ismember(i,c5)==1
        c5_cells = [c5_cells,marea(i)];
        c5_count = [c5_count,count(i)];
    elseif ismember(i,c6)==1
        c6_cells = [c6_cells,marea(i)];
        c6_count = [c6_count,count(i)];
    elseif ismember(i,c7)==1
        c7_cells = [c7_cells,marea(i)];
        c7_count = [c7_count,count(i)];
    elseif ismember(i,c8)==1
        c8_cells = [c8_cells,marea(i)];
        c8_count = [c8_count,count(i)];
    elseif ismember(i,c9)==1
        c9_cells = [c9_cells,marea(i)];
        c9_count = [c9_count,count(i)];
    elseif ismember(i,c10)==1
        c10_cells = [c10_cells,marea(i)];
        c10_count = [c10_count,count(i)];
    elseif ismember(i,c11)==1
        c11_cells = [c11_cells,marea(i)];
        c11_count = [c11_count,count(i)];
    elseif ismember(i,c12)==1
        c12_cells = [c12_cells,marea(i)];
        c12_count = [c12_count,count(i)];
    elseif ismember(i,c13)==1
        c13_cells = [c13_cells,marea(i)];
        c13_count = [c13_count,count(i)];
    elseif ismember(i,c14)==1
        c14_cells = [c14_cells,marea(i)];
        c14_count = [c14_count,count(i)];
    end
end
controlmean = mean(control_cells);
controlsem = std(control_cells)/sqrt(length(control_cells));

%%
c1mean = mean(c1_cells);
c1sem = std(c1_cells)/sqrt(length(c1_cells));
c2mean = mean(c2_cells);
c2sem = std(c2_cells)/sqrt(length(c2_cells));
c3mean = mean(c3_cells);
c3sem = std(c3_cells)/sqrt(length(c3_cells));
c4mean = mean(c4_cells);
c4sem = std(c4_cells)/sqrt(length(c4_cells));
c5mean = mean(c5_cells);
c5sem = std(c5_cells)/sqrt(length(c5_cells));
c6mean = mean(c6_cells);
c6sem = std(c6_cells)/sqrt(length(c6_cells));
c7mean = mean(c7_cells);
c7sem = std(c7_cells)/sqrt(length(c7_cells));
c8mean = mean(c8_cells);
c8sem = std(c8_cells)/sqrt(length(c8_cells));
c9mean = mean(c9_cells);
c9sem = std(c9_cells)/sqrt(length(c9_cells));
c10mean = mean(c10_cells);
c10sem = std(c10_cells)/sqrt(length(c10_cells));
c11mean = mean(c11_cells);
c11sem = std(c11_cells)/sqrt(length(c11_cells));
c12mean = mean(c12_cells);
c12sem = std(c12_cells)/sqrt(length(c12_cells));
c13mean = mean(c13_cells);
c13sem = std(c13_cells)/sqrt(length(c13_cells));
c14mean = mean(c14_cells);
c14sem = std(c14_cells)/sqrt(length(c14_cells));

x = categorical({controlname,c1name,c2name,c3name,c4name,c5name,c6name,c7name,c8name,c9name,c10name,c11name,c12name,c13name,c14name});
x = reordercats(x,{controlname,c1name,c2name,c3name,c4name,c5name,c6name,c7name,c8name,c9name,c10name,c11name,c12name,c13name,c14name});
cellarea = [controlmean,c1mean,c2mean,c3mean,c4mean,c5mean,c6mean,c7mean,c8mean,c9mean,c10mean,c11mean,c12mean,c13mean,c14mean];
cellsem = [controlsem,c1sem,c2sem,c3sem,c4sem,c5sem,c6sem,c7sem,c8sem,c9sem,c10sem,c11sem,c12sem,c13sem,c14sem];
%p = zeros(1,6); err = zeros(1,6); coloc = zeros(6,10); cms = zeros(6,10);


figure
bar(x,cellarea)
hold on
er = errorbar(x,cellarea,cellsem,cellsem);
er.LineStyle = 'none';
ylabel('Cardiomyocyte Cell Area')

outputtoR = [controlmean,c10mean,c11mean,c12mean,c14mean,c1mean,c4mean,c7mean];
outputtoR_sem = [controlsem,c10sem,c11sem,c12sem,c14sem,c1sem,c4sem,c7sem];

%% Toxicity

control_count = mean(count(control)); controlstd = std(count(control));
c1_count = mean(count(c1)); c1std = std(count(c1));
c2_count = mean(count(c2)); c2std = std(count(c2));
c3_count = mean(count(c3)); c3std = std(count(c3));
c4_count = mean(count(c4)); c4std = std(count(c4));
c5_count = mean(count(c5)); c5std = std(count(c5));
c6_count = mean(count(c6)); c6std = std(count(c6));
c7_count = mean(count(c7)); c7std = std(count(c7));
c8_count = mean(count(c8)); c8std = std(count(c8));
c9_count = mean(count(c9)); c9std = std(count(c9));
c10_count = mean(count(c10)); c10std = std(count(c10));
c11_count = mean(count(c11)); c11std = std(count(c11));
c12_count = mean(count(c12)); c12std = std(count(c12));
c13_count = mean(count(c13)); c13std = std(count(c13));
c14_count = mean(count(c14)); c14std = std(count(c14));

xx = categorical({c12name,c9name,c8name,c7name,c14name});
xx = reordercats(xx,{c12name,c9name,c8name,c7name,c14name});
yy = [c12mean,c9mean,c8mean,c7mean,c14mean];
ss = [c12std,c9std,c8std,c7std,c14std];

%"TGFb","TGFb + 20nM Midostaurin","TGFb + 80nM Midostaurin","TGFb + 320nM Midostaurin","320nM Midostaurin"
cellcount = [c12_count,c9_count,c8_count,c7_count,c14_count];
cellstd = [c12std,c9std,c8std,c7std,c14std];
%p = zeros(1,6); err = zeros(1,6); coloc = zeros(6,10); cms = zeros(6,10);

figure
bar(xx,cellcount)
hold on
er = errorbar(xx,cellcount,cellstd,cellstd);
er.LineStyle = 'none';
ylabel('Cardiomyocyte Count')

xxx = categorical({c11name,c6name,c5name,c4name,c14name});
xxx = reordercats(xxx,{c11name,c6name,c5name,c4name,c14name});
yyy = [c11mean,c6mean,c5mean,c4mean,c14mean];
sss = [c11std,c6std,c5std,c4std,c14std];

%"TGFb","TGFb + 20nM Midostaurin","TGFb + 80nM Midostaurin","TGFb + 320nM Midostaurin","320nM Midostaurin"
cellcount = [c11_count,c6_count,c5_count,c4_count,c14_count];
cellstd = [c11std,c6std,c5std,c4std,c14std];
%p = zeros(1,6); err = zeros(1,6); coloc = zeros(6,10); cms = zeros(6,10);

figure
bar(xxx,cellcount)
hold on
er = errorbar(xxx,cellcount,cellstd,cellstd);
er.LineStyle = 'none';
ylabel('Cardiomyocyte Count')


%% statistics

area_stats = vertcat(c12_cells',c9_cells',c8_cells',c7_cells',c14_cells');
area_stats_std = [c12std,c9std,c8std,c7std,c14std];

%% plate controls summary
p_control = control_cells;
p_Ang = c10_cells;
p_NE = c11_cells;
p_TGFb = c12_cells;
plate = [p_control; p_Ang; p_NE; p_TGFb];

tot_control = [2920.69545454545,2583.87297297297,4231.48305084746,2722.98784194529,3551.73333333333,2340.73271889401,2937.40909090909,263.250000000000,2470.66666666667,826.700000000000,4772.71428571429,1056.28571428571];
tot_Ang = [1874.65761689291,2454.15146299484,2645.98790322581,2245.13084112150,559.166666666667,3460.88888888889,2430.29090909091,2355.73897058824,2380.23834196891,2563.09480122324,2444.59183673469,2898.26123595506];
tot_NE = [3004.41935483871,2199.72285251216,3453.30935251799,2422.92448979592,4017.60000000000,3664.55555555556,3956.63829787234,2393.15284974093,2348.74579831933,2889.48750000000,2022.04411764706,2599.06951871658];
tot_TGFb = [2657.19175627240,3366.69921875000,2132.78521617852,3332.97484276730,1682.98888888889,3014.16239316239,2466.49414062500,3829.64417177914,2684.54046242775,2794.65963060686,2839.08609271523,407.810810810811];

tot_control_mean = mean(tot_control); 
tot_Ang_mean = mean(tot_Ang);
tot_NE_mean = mean(tot_NE);
tot_TGFb_mean = mean(tot_TGFb) ;

tot_control_sem = std(tot_control)/sqrt(length(tot_control));
tot_Ang_sem = std(tot_Ang)/sqrt(length(tot_Ang));
tot_NE_sem = std(tot_NE)/sqrt(length(tot_NE));
tot_TGFb_sem = std(tot_TGFb)/sqrt(length(tot_TGFb));