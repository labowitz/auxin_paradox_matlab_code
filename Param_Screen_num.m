clear all; close all;
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
Nu1 = 100 ./ (N_cap * 24); %Ncap cells secrete 100uM auxin a day 
Delta1 = log(2) / 24; %24 hour half life
A_cap = Nu1.*N_cap./Delta1;
nu = (3.15 + 2.69) / 2; %Correction for generalized logistic dynmaics

%enumerate varied Parameters
APs = logspace(0, 4, 51);
Bs = logspace(0, 3, 51);
Rs = logspace(log10(P_R0)-1, log10(P_R0)+1, 3);
Cs = logspace(log10(P_C0)-1, log10(P_C0)+1, 3);
[AP,B_ext,P_R,P_C] = ndgrid(APs, Bs, Rs, Cs);
% 
% Parameter Sweep Set up: 


wiredParam = [];
tic
parfor i = 1:numel(AP)
    % Factor in the AP1903 and Blasticidin concentration, and protein expressions into the parameters
    theta = param;
    theta(3) = param(3).* (P_C(i) ./ P_C0) .* (AP(i) ./ AP0).^(1/param(8));
    theta(5) = param(5) .* (B_ext(i) ./ B_ext0);
    theta(6) = param(6) .* (P_R(i) ./ P_R0);   
        
    rateParamedFunc = @(A) rate_syn(theta, phi_s, A);
    grTypes(i) = GrowthCurveType(rateParamedFunc, A_cap);
    
    if grTypes(i) == -2
        wiredParam = [wiredParam; theta];
    end
end  

elapsed1 = toc;
disp(elapsed1./60.*minutes) 
% Solve for equilbrium points of rate function :


grTypes = reshape(grTypes, [length(APs), length(Bs), length(Rs), length(Cs)]);

close all;

MarkerSize = 8.5;
figure(1);

colors = containers.Map([-2, -1, 0, 1, 2, 3], ... 
                        {[.1, .1, .1], [.2, 0, 1], ...
                         [.2, 1, .2], [.9, 0, 0], ...
                         [.9, .9, .1], [.1, .9, 1]});


[AP_sub_mesh, B_sub_mesh] = ndgrid(APs, Bs);

for a = 1:(length(Rs)*length(Cs))
    subplot(length(Rs), length(Cs), a, 'replace');
    hold on;
    
    curTypes = grTypes(:, :, a);
    for grCode = [-2, -1, 0, 1, 2, 3]
        thisCode = find(curTypes == grCode);
        drawY = AP_sub_mesh(thisCode);
        drawX = B_sub_mesh(thisCode);
        plot(drawX,drawY,'o','MarkerSize',MarkerSize,'MarkerFaceColor',colors(grCode),'MarkerEdgeColor',colors(grCode))
    end

    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';
    
    hold off;
end

auxins = logspace(-3, log10(150), 1000);
figure(2);
hold on;
for i =1:size(wiredParam, 1)
    theta = wiredParam(i, :);
    rates = rate_syn(theta, phi_s, auxins);
    plot(auxins, rates);
end
ax = gca;
ax.XScale = 'log';
hold off;
    
    