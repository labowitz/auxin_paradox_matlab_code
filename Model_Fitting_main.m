%% 
% Loading data, normalization and converting to growth rate

clear all; close all;
addpath('./MatlabFuncs/')

xdata = [0, .2, .5, 1, 2, 5, 10, 20, 50, 100];

% dataTable = readtable('Data\MeanSTD_Fitting.xlsx');

cellLine = 1; 
%Cellline 1 => sc13; 2 => sc35

dataFile = sprintf('Data/MeanSTD_Fitting%d.xlsx', cellLine);
dataTable = readtable(dataFile); %%mac syntax
dataTable = table2array(dataTable);
y_data = dataTable(1:4, 2:end);
y_error = dataTable(5:8, 2:end);

% Converting y to growth rate
Y0 = 1/12; K = 1; t=96; nu = 2.69;
save('growthrate_fixedparam','nu','K','t','Y0','-mat');

r_data = revLogest(y_data, Y0, K, t, nu);
% disp(r_data);
r_diff = r_data(1, :) - mean(r_data(1, :));
r_norm = r_data - r_diff;
% disp(r_norm);
y_norm = endLogest(r_norm, Y0, 1, t, nu);
% disp(y_norm);
%% 
% Fitting the Blast and Casp9 arm together.

run('./MatlabFuncs/FuncGenerator.m');

options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',100000, 'MaxIterations',40000); 
x0_both = [0.03, 0.1, 0.5, 0.1, 1, 2, 2, 1];

lowerLim = repelem(0, length(x0_both));
upperLim = repelem(inf, length(x0_both));

if cellLine == 1
    phiG = 0.0395;
    lowerLim(1) = phiG;
    upperLim(1) = phiG;
end


[x_both, renorm_both, residual_both] = lsqcurvefit(rate_func2,x0_both,xdata,[r_norm(3, :); r_norm(2, :)],lowerLim,upperLim,options);
MSE_both = mean((residual_both).^2,2); % mean percentage error
disp(x_both);
%%
rpredict = [rate_func2(x_both, xdata); rate_comb(x_both, xdata)];
cpredict = endLogest(rpredict, Y0, 1, t, nu);

color1 = '#2CA02C';
color2 = '#D62728';
color3 = '#0D79f2';


figure(1); %clf(fig);

hold on;

e1 = errorbar([1:length(xdata)],y_norm(3, :),y_error(3, :),'o','Color',color1,'MarkerSize',7,'MarkerFaceColor',color1);
e2 = errorbar([1:length(xdata)],y_norm(2, :),y_error(2, :),'o','Color',color2,'MarkerSize',7,'MarkerFaceColor',color2);
e3 = errorbar([1:length(xdata)],y_norm(4, :),y_error(4, :),'o','Color',color3,'MarkerSize',7,'MarkerFaceColor',color3);

p1 = plot([1:length(xdata)],cpredict(1,:),'.-','Color',color1,'LineWidth',2);
p2 = plot([1:length(xdata)],cpredict(2,:),'.-','Color',color2,'LineWidth',2);
p3 = plot([1:length(xdata)],cpredict(3,:),'--','Color',color3,'LineWidth',2);

ylim([0, 1.1])
xticklabels({'0','.2','.5','1','2','5','10', '20', '50', '100'});
xlabel('NAA (uM)');
ylabel({'Confluence'});
x1 = gca;
x1.FontSize = 24;

%%
rate_syn_fit = @(phi_s, xdata) (rate_syn(x_both, phi_s, xdata));

options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',100000, 'MaxIterations',40000); 

[phi_s, renorm_all, residual_all] = lsqcurvefit(rate_syn_fit, [0.01], xdata, r_norm(4, :), 0, inf, options);
disp(phi_s);
save('./Data/Param_Fit2', "x_both", 'phi_s');
%%
rpredict_all = [rate_func2(x_both, xdata); rate_syn(x_both, phi_s, xdata)];
cpredict_all = endLogest(rpredict_all, Y0, 1, t, nu);

p3 = plot([1:length(xdata)],cpredict_all(3,:),':','Color',color3,'LineWidth',2);
% 
% ylim([0, 1.1])
% xticklabels({'0','.2','.5','1','2','5','10', '20', '50', '100'});
% xlabel('NAA (uM)');
% ylabel({'Confluence'});
% x2 = gca;
% x2.FontSize = 24;
if cellLine == 1
    title('Cell line #13');
else
    title('Cell line #35');
end
box off;