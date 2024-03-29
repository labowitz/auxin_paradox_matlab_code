clear all; close all; clf;
addpath('./MatlabFuncs/');

load('Data/Param_Fit2.mat', 'x_both', 'phi_s'); %%mac syntax

run('./MatlabFuncs/FuncGenerator.m');
% Set known values:

param = x_both;

% AP1903 and blasticidin concentration
B_ext0 = 50; %ug
AP0 = 20; %uM

% BLASTR and iCasp protein expression
P_R0 = param(6);
eta_I =.01 ;% assumed
P_C0 =  param(3)./((eta_I.* AP0).^(1./param(8)));

% Logistic function parameters
N_cap = 1;
Nu1 = 4.18; %Ncap cells secrete 100uM auxin a day 
Delta1 = log(2) / 24; %24 hour half life
A_cap = Nu1.*N_cap./Delta1;
nu = 2.69; %Correction for generalized logistic dynmaics

% Enumerate varied Parameters
% APs = logspace(0, 3, 201);
% Bs = logspace(0, log10(200), 201);

% Enumerate varied Parameters (Sub-sampling)
APs = logspace(0, 3, 501);
Bs = logspace(0, log10(200), 501);

% Linear space
% APs = linspace(0, 1e3, 31);
% Bs = linspace(0, 200, 31);

[AP,B_ext] = ndgrid(APs, Bs);
% 
% Parameter Sweep Set up: 

% screen of the AP x Blast map
grCodes = [-2, -1, 0, 1, 2, 3];
sortParams = {[], [], [], [], [], []};

tic
parfor i = 1:numel(AP)
    % Factor in the AP1903 and Blasticidin concentration, and protein expressions into the parameters
    theta = param;
    theta(3) = param(3) .* (AP(i) ./ AP0).^(1/param(8));
    theta(5) = param(5) .* (B_ext(i) ./ B_ext0);
        
    rateParamedFunc = @(A) rate_syn(theta, phi_s, A);
    grTypes(i) = GrowthCurveType(rateParamedFunc, A_cap);
end  

for i = 1:numel(AP)
    % Factor in the AP1903 and Blasticidin concentration, and protein expressions into the parameters
    theta = param;
    theta(3) = param(3) .* (AP(i) ./ AP0).^(1/param(8));
    theta(5) = param(5) .* (B_ext(i) ./ B_ext0);
    
    % Store different parameters by their type
    sortParams{grTypes(i) + 3} = [sortParams{grTypes(i) + 3}; theta];
end  

% Measure and display the time used
elapsed1 = toc;
disp(elapsed1./60.*minutes) 

grTypes = reshape(grTypes, [length(APs), length(Bs)]);

close all;

MarkerSize = 5;
figure(1);  

colors = containers.Map([-2, -1, 0, 1, 2, 3], ... 
                        {[192, 192, 192] / 255, ...
                         [183, 238, 217] / 255, ...
                         [254, 226, 187] / 255, ... 
                         [210, 195, 224] / 255, ... 
                         [245, 203, 213] / 255, ...
                         [187, 240, 254] / 255});

colorMap = [[192, 192, 192] / 255; ...
            [183, 238, 217] / 255; ...
            [254, 226, 187] / 255; ... 
            [210, 195, 224] / 255; ... 
            [245, 203, 213] / 255; ...
            [187, 240, 254] / 255];


[AP_sub_mesh, B_sub_mesh] = ndgrid(APs, Bs);

    
hold on;

% curTypes = grTypes(:, :, a);
% for grCode = grCodes
%     thisCode = find(grTypes == grCode);
%     drawY = AP_sub_mesh(thisCode);
%     drawX = B_sub_mesh(thisCode);
%     plot(drawX,drawY,'s','MarkerSize',MarkerSize,'MarkerFaceColor',colors(grCode),'MarkerEdgeColor',colors(grCode))
% end
contourf( B_sub_mesh, AP_sub_mesh, grTypes, grCodes)
colormap(colorMap)

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';

hold off;

%% Plot the "unclassified" group
auxins = logspace(-2, log10(150), 1000);
figure(2);
clf;
hold on;
for i = 1 : size(sortParams{1}, 1)
    theta = sortParams{1}(i, :); %Won't work in the linear case, set to 1 instead
    rates = rate_syn(theta, phi_s, auxins);
    plot(auxins, rates, 'color', [0.5, 0.5, 0.5], 'LineWidth', 3);
end


plot(auxins, rate_syn(sortParams{6}(2, :), phi_s, auxins), 'color', colors(3), 'LineWidth', 3);
plot(auxins, rate_syn(sortParams{5}(75, :), phi_s, auxins), 'color', colors(2), 'LineWidth', 3);


plot(auxins, auxins * 0, '--k', 'LineWidth', 3);
ax = gca;
ax.XScale = 'log';
hold off;

%% Plot the examples
auxins = logspace(-3, log10(150), 1000);
figure(2);
clf;
hold on;
for i = [2, 4, 5, 3, 6]
    theta = sortParams{i}(37, :); %Won't work in the linear case, set to 1 instead
    rates = rate_syn(theta, phi_s, auxins);
    plot(auxins, rates, 'color', colors(i - 3), 'LineWidth', 3);
end
plot(auxins, auxins * 0, '--k', 'LineWidth', 3);
ax = gca;
ax.XScale = 'log';
hold off;

%% Screen expression profiles for paradoxical percentage

RC_total = logspace(log10(P_R0 + P_C0)-1, log10(P_R0 + P_C0)+1, 15);
RC_ratio = logspace(-1, 1, 15);

% replace AP and B with a downsampled ones
APs = logspace(0, 3, 21);
Bs = logspace(0, log10(200), 21);

[AP,B_ext,tots,rats] = ndgrid(APs, Bs, RC_total, RC_ratio);

tic
parfor i = 1:numel(AP)
    % Factor in the AP1903 and Blasticidin concentration, and protein expressions into the parameters
    P_R_i = tots(i) * rats(i) / (rats(i) + 1);
    P_C_i = tots(i) / (rats(i) + 1);
    theta = param;
    theta(3) = param(3).* (P_C_i ./ P_C0) .* (AP(i) ./ AP0).^(1/param(8));
    theta(5) = param(5) .* (B_ext(i) ./ B_ext0);
    theta(6) = param(6) .* (P_R_i ./ P_R0); 
    
    rateParamedFunc = @(A) rate_syn(theta, phi_s, A);
    grTypes2(i) = GrowthCurveType(rateParamedFunc, A_cap);
end  

% Measure and display the time used
elapsed1 = toc;
disp(elapsed1./60.*minutes) 

grTypes2 = reshape(grTypes2, [length(APs), length(Bs), length(RC_total) * length(RC_ratio)]);

for i = 1 : length(RC_total) * length(RC_ratio)
    paraDens(i) = sum((grTypes2(:,:,i) == 3), 'all') / (length(APs) * length(Bs));
end
paraDens = reshape(paraDens, [length(RC_total), length(RC_ratio)]);

figure(3);
heatmap(paraDens, 'GridVisible','off');