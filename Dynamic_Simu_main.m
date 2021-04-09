%% Load data

clear all; close all;
addpath('./MatlabFuncs/')

load('Data/Param_Fit2.mat', 'x_both', 'phi_s'); %%mac syntax

xdata = [0, .2, .5, 1, 2, 5, 10, 20, 50, 100];

run('./MatlabFuncs/FuncGenerator.m');

auxins = logspace(-2,2,100);
param = x_both; 

%% Solve for threshold
%Set known values
B_ext = 50; %ug
AP = 20; %uM
Nu1 = 4.18; %Fitted auxin secretion rate

A01 = fzero(@(A) rate_syn(param, phi_s, A), 0);
A02 = fzero(@(A) rate_syn(param, phi_s, A), 100);
N01 = A01 * log(2) / 24 / Nu1;
N02 = A02 * log(2) / 24 / Nu1;

%% Time Trajectories 
% _Varied initial Conditions_

tspan = linspace(0, 80 * 24, 1000);

c0 = logspace(-4, -0.25, 10);
a0 = A01;
[A1, N1] = meshgrid(a0,c0);

figure(1);

cir_on = [];
cir_off = [];

for i = 1 : numel(A1)
    y0 = [A1(i), N1(i)];
    [t,f] = ode45(@(t,y) DynSys_syn(t,y,param,phi_s,B_ext,AP),tspan,y0');
    
    [tc,fc] = ode45(@(t,y) DynSys_syn(t,y,param,phi_s,0,0),tspan,y0');
    cir_on = [cir_on, f(:,2)];
    cir_off = [cir_off, fc(:,2)];
    
    hold on;
    plot(tspan./24,f(:,2),'k','Linewidth',5,'Color',[100,255,201]./255)
    plot(tspan./24,fc(:,2),'k','Linewidth',5,'Color',[1,0,0])
end

plot(tspan/24, repelem(N01, length(tspan)), 'k--');

hold off;

save("Data\Dynamcs.mat", "cir_on", "cir_off", "tspan", 'N01');

xlabel('Time (days)')
ylabel('Cell Population')

ax1 = gca;
ax1.FontSize = 25;
%% Delay Differential Equation
  
figure(4)
tau = [.25,8,24,48,72];
clist = colormap(cool(length(tau)));
for i = 1 : numel(tau)
    sol = dde23(@(t,y,Z) DynSys_syn_delay2(t,y,Z,param,phi_s,B_ext,AP), tau(i), @history, tspan);
    t_delay = sol.x';
    f_delay = sol.y';
    hold on
    plot(t_delay./24,f_delay(:,2),'k','Linewidth',3,'Color',clist(i,:))                         

end
c = colorbar('eastoutside','Ticks',[.1,.3,.5,.7,.9],...
         'TickLabels',{'.25','8 ','24','48','72'},'FontSize',18);
c.Label.String = '\tau (hr)';

y00 = [0.3300 .1];
[t0,f0] = ode45(@(t,y) DynSys_syn(t,y,param,phi_s,B_ext,AP),tspan,y00');
plot(tspan./24,f0(:,2),'k--','Linewidth',2,'Color',[50,50,201]./255)


ylim([0,1]);
xlim([0,t_delay(end)./24])
xlabel('Time (days)')
ylabel('Cell Population')

ax4 = gca;
ax4.FontSize = 25;
% ax4.YScale = 'log';




%% _Paradoxical Feedback Response to Mutation_

mutSpike = 0.01;

T = 25;
y0 = [10, 0.1, 0];

t_all1 = [];
f_all1 = [];
y0_t = y0;
for tp = [0,T,T*2,T*3]
    tspan1 = linspace(tp*24,(tp+T)*24,100);
    [t1,f1] = ode45(@(t,y) DynSys_syn_mut(t,y,param,phi_s,B_ext,AP),tspan1,y0_t');
    t_all1 = [t_all1; t1];
    f_all1 = [f_all1; f1];
    y0_t = f1(end,:) + [0, 0, mutSpike];
end


figure(2)
clf; hold on;

plot(t_all1./24, f_all1(:,2), 'LineWidth', 8, 'Color', [100,255,201]./255);
plot(t_all1./24, f_all1(:,3), 'LineWidth', 8, 'Color', [1,0,0]);

xlabel('Time (days)')
ylabel('Cell Population')
% xlim([0, 60])
ax2 = gca;
ax2.FontSize = 25;
box off;
%% _Simple Feedback Response to Mutation_

param_sf = param;
param_sf(2) = 0;


t_all2 = [];
f_all2 = [];
y0_t = y0;
for tp = [0,T,T*2,T*3]
    tspan1 = linspace(tp*24,(tp+T)*24,100);
    [t1,f1] = ode45(@(t,y) DynSys_syn_mut(t,y,param_sf,0,B_ext,AP),tspan1,y0_t');
    t_all2 = [t_all2; t1];
    f_all2 = [f_all2; f1];
    y0_t = f1(end,:) + [0, 0, mutSpike];
end


save("Data\Dynamcs_mut1.mat", 't_all1', 'f_all1', "t_all2", 'f_all2')

figure(3)
clf; hold on;

plot(t_all2./24, f_all2(:,2), 'LineWidth', 8, 'Color', [100,255,201]./255);
plot(t_all2./24, f_all2(:,3), 'LineWidth', 8, 'Color', [1,0,0]);

xlabel('Time (days)')
ylabel('Cell Population')
xlim([0, 60])
ax3 = gca;
ax3.FontSize = 25;
box off;
%% History Function for delay differential equation

function s = history(t) % history function for t <= 0
  s = ones(2,1);
  s(1) = 0.3300;
  s(2) = .1;
end