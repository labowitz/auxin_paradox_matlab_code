%% Fit for Line S1-13
% 
% Full System*:
% $$log(\hbox{ydata}) - log\left(\frac{N_0}{N_{confluent}}\right) =  (\psi_g 
% + \psi_d)\cdot\frac{(\eta_AA + 1)^4 + \left(\frac{P_R}{\delta}\right)^4}{(\eta_B 
% B_{ext} + 1)(\eta_AA+ 1)^4+ \left( \frac{P_R}{\delta} \right)^4 } - \psi_d-\psi_C 
% \cdot \frac{\eta_{I} I \cdot  \left(\frac{P_C}{\delta}\right)^2}{\eta_{I}I \cdot 
% \left(\frac{P_C}{\delta}\right)^2+(\eta_AA +1)^2} $$
% Blast =0:
% $$log(\hbox{ydata}) - log\left(\frac{N_0}{N_{confluent}}\right) = \psi_g-\psi_C 
% \frac{I^*\left(\frac{P_C}{\delta}\right)^2}{I^*\left(\frac{P_C}{\delta}\right)^2+(\eta_AA 
% +1)^2}$$
% AP1903 =0: 
% $$log(\hbox{ydata}) - log\left(\frac{N_0}{N_{confluent}}\right) =  (\psi_g 
% + \psi_d)\cdot\frac{(\eta_AA + 1)^4 + \left(\frac{P_R}{\delta}\right)^4}{(\eta_B 
% B_{ext} + 1)(\eta_AA+ 1)^4+ \left( \frac{P_R}{\delta} \right)^4 } - \psi_d$$
% 
% 
% 
% *Note that ydata is the ratio of cell number to confluent cell number 
% 
% Input the observation axuin levels:

clear all; close all

xdata = ...
 [ 0 .2 .5 1 2 5 10 20 50 100];
%% 
% Input the  responses Line S1-3 for just BLAST:

%%NoAP1903
ydata1 = ...
    [0.994351247	0.9916548	0.989389418	0.996206878	0.994587755	0.983118519	0.893697128	0.400989342	0.006447241	0.00148715];

% Fit the data to the model. 
% No AP1903 model:
% 
% AP1903 =0: 
% $$log(\hbox{ydata}) - log\left(\frac{N_0}{N_{confluent}}\right) =  (\psi_g 
% + \psi_d)\cdot\frac{(\eta_AA + 1)^4 + \left(\frac{P_R}{\delta}\right)^4}{(\eta_B 
% B_{ext} + 1)(\eta_AA+ 1)^4+ \left( \frac{P_R}{\delta} \right)^4 } - \psi_d$$
% 
% x(1) = $\psi_g =\delta$ (fixed at 0.0286 hr^-1 -->  24 hrs doubling time ) 
% , x(2) = $\eta_A$, x(3) = $\psi_C$, x(4) = $I^* P_C^2$, x(5) = $\psi_d$,  x(6) 
% = $B^* =\eta_B B_{\textrm{ext}}$, x(7) = $P_R$ 

rate1 = @(x,xdata)((((x(1) + x(5)).*((x(2).*xdata +1).^4 + (x(7)./x(1)).^4)./((x(6) + 1).*(x(2).*xdata +1).^4 + (x(7)./x(1)).^4)) - x(5)).*50);

fun1 = @(x,xdata)(exp(rate1(x,xdata) - log(4)));


 options = optimoptions('lsqcurvefit','MaxFunctionEvaluations',25000, 'MaxIterations',4000);
x01 = [0.0286 2 .01 .1 1 1 2];
% 
 [x1,renorm1,residual1] = lsqcurvefit(fun1,x01,xdata,ydata1,[0.0286 1.5 0 0 0 0 0],[0.0286  1.5 inf inf inf inf inf],options);
MSE1 = mean((residual1).^2); % mean percentage error 
MPE1 = mean((residual1)./ydata1).*100; % mean percentage error 
%% 
% 
% Plot the data and the fitted curve.

times =xdata; 

err1 = [0.002809344	0.005162316	0.004957646	0.001467934	0.002520465	0.019980023	0.033964324	0.054158992	0.002579971	0.000759584];
figure(1)
errorbar([1:length(xdata)],ydata1,err1,'o','Color',[0,.7,0],'MarkerSize',7,...
    'MarkerEdgeColor',[0,.7,0],'MarkerFaceColor',[0,.7,0])
hold on 
plot([1:length(xdata)],fun1(x1,times),'.-','Color',[0,.7,0],'MarkerFaceColor',[0,.7,0],'LineWidth',2)
xticklabels({'0','.2','.5','1','2','5','10', '20', '50', '100'})
xlabel('NAA (\mu M)')
ylabel({'Confluence'})
ax1 = gca;
ax1.FontSize = 24;
ylim([0,1.1])
% ax.FontWeight = 'bold';
box off
legend('Experimental Data','Model Fit','Location','SouthWest')

 text(.5,.4,['Mean Squared Error $ \approx$ ',num2str(round(MSE1,3))],'Interpreter','Latex','FontSize',16)

 text(.5,.55,'$(\psi_g + \psi_d)\cdot\frac{(\eta_AA + 1)^4 + (\frac{P_R}{\delta})^4}{(\eta_BB_{ext} +1)(\eta_AA + 1)^4 + (\frac{P_R}{\delta})^4 }-\psi_d$','Interpreter','Latex','FontSize',18)

% Input the  responses Line S1-3 for just BLAST:

%%NoBlast
 ydata2 = ...
    [0.020884656	0.038708919	0.053608692	0.0795322	0.106651625	0.426742328	0.801811791	0.927911111	0.978334618	0.916049206];

% Fit the data to the model. 
% Blast =0:
% $$log(\hbox{ydata}) - log\left(\frac{N_0}{N_{confluent}}\right) = \psi_g-\psi_C 
% \frac{I^*\left(\frac{P_C}{\delta}\right)^2}{I^*\left(\frac{P_C}{\delta}\right)^2+(\eta_AA 
% +1)^2}$$
% 
% x(1) = $\psi_g =\delta$, x(2) = $\eta_A$, x(3) = $\psi_C$, x(4) = $I^* P_C^2$, 
% x(5) = $\psi_d$,  x(6) = $B^*$, x(7) = $P_R$ 

rate2 = @(x,xdata)((x(1)- x(3).*((x(4)./x(1)).^2./((x(4)./x(1)).^2 + (x(2).*xdata + 1).^2))).*50);
 
fun2 = @(x,xdata)(exp(rate2(x,xdata) - log(4)));

x02 = [x1(1) x1(2) .01 .1 x1(5) x1(6) x1(7)];

[x2,renorm2,residual2] = lsqcurvefit(fun2,x02,xdata,ydata2,[x1(1) x1(2) 0 0 x1(5) x1(6) x1(7)],[x1(1) x1(2)  inf inf x1(5) x1(6) x1(7)],options);
 MSE2 = mean((residual2).^2);  % mean squared error 
%  MPE2 = mean((residual2)./ydata2).*100; % mean percentage error 
% Plot the data and the fitted curve.


err2 = [0.010813322	0.018252794	0.026223576	0.012246634	0.031479622	0.16338727	0.171613373	0.059727441	0.01959887	0.03884396];
figure(2)
errorbar([1:length(xdata)],ydata2,err2,'o','Color',[1,.7,0],'MarkerSize',7,...
    'MarkerEdgeColor',[1,.7,0],'MarkerFaceColor',[1,.7,0])
hold on 
plot([1:length(xdata)],fun2(x2,times),'.-','Color',[1,.7,0],'MarkerFaceColor',[1,.7,0],'LineWidth',2)
xticklabels({'0','.2','.5','1','2','5','10', '20', '50', '100'})
xlabel('NAA (\mu M)')
ylabel({'Confluence'})
ax1 = gca;
ax1.FontSize = 24;
ylim([0,1.1])
% ax.FontWeight = 'bold';
box off
legend('Experimental Data','Model Fit','Location','NorthWest')
 text(6.3,.2,['Mean Squared Error $ \approx$ ',num2str(round(MSE2,3))],'Interpreter','Latex','FontSize',16)
 
 text(6.3,.4,'$\psi_g-\delta_N \frac{\eta_II\left(\frac{P_C}{\delta}\right)^2}{\eta_II\left(\frac{P_C}{\delta}\right)^2+(\eta_{A}A +1)^2}$','Interpreter','Latex','FontSize',18)
% Input the  responses Line S1-3 for just BLAST:

%Full
ydata3 = ...
[0.028152532	0.050824187	0.055106803	0.078311489	0.12464966	0.216782237	0.052270824	0.020610733	0.004709373	0.000569237];

% Fit the data to the model. 
% Full System:
% $$log(\hbox{ydata}) - log\left(\frac{N_0}{N_{confluent}}\right) =  (\psi_g 
% + \psi_d)\cdot\frac{(\eta_AA + 1)^4 + \left(\frac{P_R}{\delta}\right)^4}{(\eta_B 
% B_{ext} + 1)(\eta_AA+ 1)^4+ \left( \frac{P_R}{\delta} \right)^4 } - \psi_d-\psi_C 
% \cdot \frac{\eta_{I} I \cdot  \left(\frac{P_C}{\delta}\right)^2}{\eta_{I}I \cdot 
% \left(\frac{P_C}{\delta}\right)^2+(\eta_AA +1)^2} $$
% 
% x(1) = $\psi_g =\delta$, x(2) = $\eta_A$, x(3) = $\psi_C$, x(4) = $I^* P_C^2$, 
% x(5) = $\psi_d$,  x(6) = $B^*$, x(7) = $P_R$ 


rate3 = @(x,xdata)((((x(1) + x(5)).*((x(2).*xdata +1).^4 + (x(7)./x(1)).^4)./((x(6) + 1).*(x(2).*xdata +1).^4 + (x(7)./x(1)).^4)) - x(5) - x(3).*((x(4)./x(1)).^2./((x(4)./x(1)).^2 + (x(2).*xdata + 1).^2))  ).*50);
 
fun3 = @(x,xdata)(exp(rate3(x,xdata) - log(4)));

x03 = x2;

[x3,renorm3,residual3] = lsqcurvefit(fun3,x03,xdata,ydata3,[x2(1) x2(2) x2(3) x2(4) x2(5) x2(6) x2(7)],[x2(1) x2(2) x2(3) x2(4) x2(5) x2(6) x2(7)],options);
 MSE3 = mean((residual3).^2);  % mean squared error 
%  MPE3 = mean((residual3)./ydata3).*100; % mean percentage error 
% Plot the data and the fitted curve.


err3 = [0.014081056	0.012390267	0.020057457	0.012582435	0.020161014	0.025438707	0.021837773	0.007677594	0.004816866	0.00058251];
figure(3)
errorbar([1:length(xdata)],ydata3,err3,'o','Color',[1,0,0],'MarkerSize',7,...
    'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])
hold on 
plot([1:length(xdata)],fun3(x3,times),'.-','Color',[1,0,0],'MarkerFaceColor',[1,0,0],'LineWidth',2)
xticklabels({'0','.2','.5','1','2','5','10', '20', '50', '100'})
xlabel('NAA (\mu M)')
ylabel({'Confluence'})
ax1 = gca;
ax1.FontSize = 24;
ylim([0,1.1])
% ax.FontWeight = 'bold';
box off
legend('Experimental Data','Model Prediction','Location','Best')
text(1,.75,['Mean Squared Error $ \approx$ ',num2str(round(MSE3,3))],'Interpreter','Latex','FontSize',14)
 
 text(.4,.65,'$(\psi_g + \psi_d)\cdot\frac{(\eta_AA + 1)^4 + (\frac{P_R}{\delta})^4}{(\eta_BB_{ext} +1)(\eta_AA + 1)^4 + (\frac{P_R}{\delta})^4 }-\psi_d $','Interpreter','Latex','FontSize',15)
 text(2,.5,'$-\delta_N \frac{\eta_II\left(\frac{P_C}{\delta}\right)^2}{\eta_II\left(\frac{P_C}{\delta}\right)^2+(\eta_{A}A +1)^2}$','Interpreter','Latex','FontSize',15)
% title('AP1903')
%%
figure(4)
errorbar([1:length(xdata)],ydata3,err3,'o','Color',[1,0,0],'MarkerSize',7,...
    'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])
hold on 
plot([1:length(xdata)],fun3(x2,times),'.-','Color',[1,0,0],'MarkerFaceColor',[1,0,0],'LineWidth',2)

errorbar([1:length(xdata)],ydata2,err2,'o','Color',[1,.7,0],'MarkerSize',7,...
    'MarkerEdgeColor',[1,.7,0],'MarkerFaceColor',[1,.7,0])
hold on 
plot([1:length(xdata)],fun2(x2,times),'.-','Color',[1,.7,0],'MarkerFaceColor',[1,.7,0],'LineWidth',2)


errorbar([1:length(xdata)],ydata1,err1,'o','Color',[0,.7,0],'MarkerSize',7,...
    'MarkerEdgeColor',[0,.7,0],'MarkerFaceColor',[0,.7,0])
hold on 
plot([1:length(xdata)],fun1(x1,times),'.-','Color',[0,.7,0],'MarkerFaceColor',[0,.7,0],'LineWidth',2)


xticklabels({'0','.2','.5','1','2','5','10', '20', '50', '100'})
xlabel('NAA (\mu M)')
ylabel({'Confluence'})
ax1 = gca;
ax1.FontSize = 24;
ylim([0,1.1])
% ax.FontWeight = 'bold';
box off
% Change the Sensitivity: 
% 
% Fit the data to the model. 
% Full System:
% $$log(\hbox{ydata}) - log\left(\frac{N_0}{N_{confluent}}\right) =  (\psi_g 
% + \psi_d)\cdot\frac{(\eta_AA + 1)^4 + \left(\frac{P_R}{\delta}\right)^4}{(\eta_B 
% B_{ext} + 1)(\eta_AA+ 1)^4+ \left( \frac{P_R}{\delta} \right)^4 } - \psi_d-\psi_C 
% \cdot \frac{\eta_{I} I \cdot  \left(\frac{P_C}{\delta}\right)^2}{\eta_{I}I \cdot 
% \left(\frac{P_C}{\delta}\right)^2+(\eta_AA +1)^2} $$
% 
% x(1) = $\psi_g =\delta$, x(2) = $\eta_A$, x(3) = $\psi_C$, x(4) = $\eta_I 
% \cdot I{\cdot P}_C^2$, x(5) = $\psi_d$,  x(6) = $\eta_B B_{\textrm{ext}}$, x(7) 
% = $P_R$ 
% 
% Rewrite x(4) and x(6) to expose death rate sensitivities $\kappa$and $\kappa_C$.**
% 
% x(4) = $\eta_I \cdot I{\cdot P}_C^2 =\frac{I\cdot \lambda_C^2 }{\kappa_C }$, 
% x(6) = $\eta_B B_{\textrm{ext}} =\frac{B_{\textrm{ext}} }{\kappa }$
% 
% reassign x*(4) = r and x*(6) = \frac{1}{ \kappa} and constrain $\kappa =r\cdot 
% \kappa_C$ and 
% 
% x(6) = $\eta_B B_{\textrm{ext}} =\frac{B_{\textrm{ext}} }{\kappa }=x^* \left(6\right)\cdot 
% B_{\textrm{ext}\;}$, x(4) = $r\frac{I\cdot \lambda_C^2 }{\kappa }$=$x^* \left(4\right){\cdot 
% x}^* \left(6\right)$$\cdot I\cdot \lambda_C^2$
% 
% where $x^* \left(4\right)=r$ and  $x^* \left(6\right)=\frac{1}{\kappa }$
% 
% **These are from origonal equations: 
% 
% $$(\psi_g + \psi_d)  \frac{\kappa}{B + \kappa } - \psi_d - \psi_C \frac{C^2I}{C^2I 
% + \kappa_C }$$
% 
% which lead to the above equations in the Full System after appropriate substitutions. 
% 
% 

%Set known values
B_ext = .0473; %uM
I = 5; %uM
x4star(6)= x2(6)./B_ext;

% Assume \kappa_C = .01
x4star(4) = (1./.01)./x4star(6);
lambdasq = x2(4)./(x4star(4).*x4star(6).*I);


x4(6) = B_ext.*x4star(6);

x4(4) = lambdasq.*I.*x4star(4).*x4star(6);

rate3p = @(x,xdata)((((x(1) + x(5)).*((x(2).*xdata +1).^4 + (x(7)./x(1)).^4)./(((B_ext.*x(6)) + 1).*(x(2).*xdata +1).^4 + (x(7)./x(1)).^4)) - x(5) - x(3).*((lambdasq.*I.*x(4).*x(6)./x(1)).^2./((lambdasq.*I.*x(4).*x(6)./x(1)).^2 + (x(2).*xdata + 1).^2))  ).*50);
 
fun3p = @(x,xdata)(exp(rate3p(x,xdata) - log(4)));

% Least squares fit:


x04 = x3;
[x4n,renorm4,residual4] = lsqcurvefit(fun3p,x04,xdata,ydata3,[x2(1) x2(2) x2(3) x4star(4) x2(5)  x4star(6) x2(7)],[x2(1) x2(2) x2(3) x4star(4) x2(5)  x4star(6) x2(7)],options);
 MSE4 = mean((residual4).^2);  % mean squared error 
 MPE4 = mean((residual4)./ydata3).*100; % mean percentage error 
 
P_Sensitive1 = (abs(x4star(6) - x4n(6))./x4star(6)).*100; % mean percentage error 

close all
figure(5)  
errorbar([1:length(xdata)],ydata3,err3,'o','Color',[1,0,0],'MarkerSize',7,...
    'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])
hold on 
plot([1:length(xdata)],fun3p(x4n,times),'.-','Color',[1,0,0],'MarkerFaceColor',[1,0,0],'LineWidth',2)

xticklabels({'0','.2','.5','1','2','5','10', '20', '50', '100'})
xlabel('NAA (\mu M)')
ylabel({'Confluence'})
ax1 = gca;
ax1.FontSize = 24;
ylim([0,1.1])
% ax.FontWeight = 'bold';
box off
legend('Experimental Data','Model Prediction','Location','Best')
text(.8,.76,[' Death "sensitivities" $\kappa_C$ and $\kappa$  increase by $\approx$ ',num2str(round(P_Sensitive1,0)),'\%'],'Interpreter','Latex','FontSize',14)
text(.8,.68,[' Mean Squared Error $\approx$ ',num2str(round(MSE4,3))],'Interpreter','Latex','FontSize',14)
 
% 
%  Minimize wieghted sum of squares (alternate option)


 
 err3 = [0.014081056	0.012390267	0.020057457	0.012582435	0.020161014	0.025438707	0.021837773	0.007677594	0.004816866	0.00058251];
 err3mod = [0.014081056	0.012390267	0.020057457	0.012582435	0.015	0.003	0.008	0.007677594	0.004816866	0.00058251];

 W = 1./(err3mod.^2);
 wobj = @(x,xdata,y,w) sum(w .* ((fun3p(x,xdata) - y) .^ 2)./std(ydata3).^2);
 
%  x_p = fminsearch(@(x) wobj(x,xdata,Fy3,W),[x0noB,1,1,5,4])
 x_p = fmincon(@(x) wobj(x,xdata,ydata3,W),x04,[],[],[],[],[x2(1) x2(2) x2(3) x4star(4) x2(5)  0 x2(7)],[x2(1) x2(2) x2(3) x4star(4) x2(5)  inf x2(7)]);
%% 
% 



MSE5 = mean((fun3p(x_p,times) - ydata3).^2);
P_Sensitive2 = (abs(x4star(6) - x_p(6))./x4star(6)).*100; 

close all
figure(6)  
errorbar([1:length(xdata)],ydata3,err3,'o','Color',[1,0,0],'MarkerSize',7,...
    'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])
hold on 
plot([1:length(xdata)],fun3p(x_p,times),'.-','Color',[1,0,0],'MarkerFaceColor',[1,0,0],'LineWidth',2)
xticklabels({'0','.2','.5','1','2','5','10', '20', '50', '100'})
xlabel('NAA (\mu M)')
ylabel({'Confluence'})
ax1 = gca;
ax1.FontSize = 24;
ylim([0,1.1])
% ax.FontWeight = 'bold';
box off
legend('Experimental Data','Model Prediction','Location','Best')
text(.8,.76,[' Death "sensitivities" $\kappa_C$ and $\kappa$  increase by $\approx$ ',num2str(round(P_Sensitive2,0)),'$\%$'],'Interpreter','Latex','FontSize',14)
text(.8,.68,[' Mean Squared Error $\approx$ ',num2str(round(MSE5,3))],'Interpreter','Latex','FontSize',14)
save('paramnew','x_p','x3','x4n','lambdasq','I','B_ext')
%% _Discussion Points_
%% 
% * Is the wieghted fit with modified death "senstivities" valid? 
% * Are there other mechanisms that could explain the data? 
%% 
% 
%%  Growth Rate Curves
% 
% Net Growth Unmodified parameters:


rate_nofix = @(x,xdata)(((x(1) + x(5)).*((x(2).*xdata +1).^4 + (x(7)./x(1)).^4)./((x(6) + 1).*(x(2).*xdata +1).^4 + (x(7)./x(1)).^4)) - x(5) - x(3).*((x(4)./x(1)).^2./((x(4)./x(1)).^2 + (x(2).*xdata + 1).^2))  );
times =linspace(xdata(1),200,1000);

figure(7)
dy1 = rate_nofix(x3,times);

semilogx(times,dy1,'r-','LineWidth',3) 


 
xlh = xlabel('NAA (\muM)');
xlh.Position(2) = xlh.Position(2) - 0.001;
xlh.Position(1) = xlh.Position(1) + 5;

ylabel({'Net Growth Rate (hr^{-1})'})
ax1 = gca;
ax1.FontSize = 24;
ax1.XAxisLocation = 'origin';
ax1.YAxisLocation = 'origin';
% ax.FontWeight = 'bold';
box off

 xlim([.2,100])
 ylim([-.04,.03])
xticks([.2 .5 1 2 5 10 20 50 100])
%  Net Growth Modified parameters. Death "sensitivites":

rate_sensew = @(x,xdata)(((x(1) + x(5)).*((x(2).*xdata +1).^4 + (x(7)./x(1)).^4)./(((B_ext.*x(6)) + 1).*(x(2).*xdata +1).^4 + (x(7)./x(1)).^4)) - x(5) - x(3).*((lambdasq.*I.*x(4).*x(6)./x(1)).^2./((lambdasq.*I.*x(4).*x(6)./x(1)).^2 + (x(2).*xdata + 1).^2))  );

figure(8)
dy2 = rate_sensew(x_p,times);

semilogx(times,dy2,'r-','LineWidth',3) 


 
xlh = xlabel('NAA (\muM)');
xlh.Position(2) = xlh.Position(2) - 0.001;
xlh.Position(1) = xlh.Position(1) + 5;

ylabel({'Net Growth Rate (hr^{-1})'})
ax2 = gca;
ax2.FontSize = 24;
ax2.XAxisLocation = 'origin';
ax2.YAxisLocation = 'origin';
% ax.FontWeight = 'bold';
box off

 xlim([1,100])
 ylim([-.04,.01])
xticks([.2 .5 1 2 5 10 20 50 100])
% Death Rate: 

figure(9)
 
 ratedeath = @(x,xdata)((x(1)- x(3).*((x(4)./x(1)).^2./((x(4)./x(1)).^2 + (x(2).*xdata + 1).^2))));
 dy3 = ratedeath(x2,times); 

 semilogx(times,dy3,'-','LineWidth',3,'Color',[1,.7,0]) 
 xlim([.2,100])
 

 
xlh2 = xlabel('NAA (\muM)');
xlh2.Position(2) = xlh2.Position(2) - 0.001;
xlh2.Position(1) = xlh2.Position(1) + 5;

ylabel({'Death Rate (hr^{-1})'})
ax = gca;
ax.FontSize = 24;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% ax.FontWeight = 'bold';
box off

text(3,-.01,'$\psi_g-\delta_N \frac{I^*\left(\frac{P_C}{\delta}\right)^2}{I^*\left(\frac{P_C}{\delta}\right)^2+(\eta_AA +1)^2}$','Interpreter','Latex','FontSize',20)
 xticks([.2 .5 1 2 5 10 20 50 100])

 
 ylim([-.04,.03])
%%

 
figure(10)
 ratedeath = @(x,xdata)(( x(3).*((x(4)./x(1)).^2./((x(4)./x(1)).^2 + (x(2).*xdata + 1).^2))));

 dy4 = ratedeath(x2,times); 

 semilogx(times,dy4,'-','LineWidth',3,'Color',[1,.7,0]) 
 xlim([.2,100])
 
 
xlh3 = xlabel('NAA (\muM)');
% xlh3.Position(2) = xlh3.Position(2) - 0.02;
% xlh3.Position(1) = xlh3.Position(1) + 5;

ylabel({'Death Rate (hr^{-1})'})
ax = gca;
ax.FontSize = 24;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% ax.FontWeight = 'bold';
box off

text(.5,.01,'$\delta_N \frac{I^*\left(\frac{P_C}{\delta}\right)^2}{I^*\left(\frac{P_C}{\delta}\right)^2+(\eta_AA +1)^2}$','Interpreter','Latex','FontSize',20)
% title('AP1903')

 xticks([.2 .5 1 2 5 10 20 50 100])

 ylim([0,.07])
 
 
 
%%
 
 figure(11)
 rateprolif = @(x,xdata)(((x(1) + x(5)).*((x(2).*xdata +1).^4 + (x(7)./x(1)).^4)./((x(6) + 1).*(x(2).*xdata +1).^4 + (x(7)./x(1)).^4)) - x(5));
 dy5 = rateprolif(x1,times); 
 semilogx(times,dy5,'-','LineWidth',3,'Color',[0,.7,0]) 
 
% xlh4 = xlabel('NAA (\muM)');
%   xlh4.Position(2) = xlh4.Position(2) ;
%   xlh4.Position(1) = xlh4.Position(1)  ;

ylabel({'Proliferation Rate (hr^{-1})'})
ax = gca;
ax.FontSize = 24;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% ax.FontWeight = 'bold';
box off

xticks([.2 .5 1 2 5 10 20 50 100])
 xticklabels({'.2' ,'.5' ,'1 ','2' ,'5' ,'10' ,'20' ,'50','100'})


 ylim([-.02,.03])
 
  xlim([.2,100])
%% Time Trajectories 
%% _Varied initial Conditions_

tspan = linspace(0,60*24,10000);
  
c1= logspace(3,7,15);
a1 = 10;%linspace(5,20,5);
[A1,N1] = meshgrid(a1,c1);
 figure(12)
for i = 1:5:numel(A1)
y0 = [A1(i),N1(i)];
[t,f] = ode45(@(t,y,x) parodox_Elowitz1_fitnewnew(t,y,x4n,B_ext,lambdasq,I),tspan,y0');
[tc,fc] = ode45(@(t,y,x) parodox_Elowitz1_fitnewnew(t,y,x4n,0,lambdasq,0 ),tspan,y0');
 hold on
 plot(tspan./24,f(:,2),'k','Linewidth',5,'Color',[100,255,201]./255)
  plot(tspan./24,fc(:,2),'k','Linewidth',5,'Color',[1,0,0])
end

xlabel('Time (days)')
ylabel('Cell Population')
xlim([0,30])
ax = gca;
ax.FontSize = 25;

%% 
%% _Paradoxical Feedback Response to Mutation_

T1 = 25;
tspan1 = linspace(0,T1*24,10000);
y01 = [A1(numel(A1)),N1(numel(A1)),0];
[t1,f1] = ode45(@(t,y) parodox_Elowitz2_fitnewnew(t,y,x4n,B_ext,lambdasq,I),tspan1,y01');

figure(13)
plot(t1./24,(f1(:,2)),'LineWidth',8,'Color',[100,255,201]./255)
hold on
plot(t1./24,(f1(:,3)),'LineWidth',8,'Color',[1,0,0])

%1st mutation
T2 = 50;
tspan2 = linspace(T1*24,T2*24,10000);
auxin_e2 = f1(end,1);
N_e2 = f1(end,2);
N_e1 = 4e6;

y02 = [auxin_e2,N_e2,N_e1];
[t2,f2] = ode45(@(t,y) parodox_Elowitz2_fitnewnew(t,y,x4n,B_ext,lambdasq,I),tspan2,y02');
plot(t2./24,(f2(:,2)),'LineWidth',8,'Color',[100,255,201]./255)
hold on
plot([t1(end);t2]./24,([f1(end,3);f2(:,3)]),'LineWidth',8,'Color',[1,0,0])

%2nd mutation
T3 = 75;
tspan3 = linspace(T2*24,T3*24,10000);

auxin_e3 = f2(end,1);
N_e3 = f2(end,2);

y03 = [auxin_e3,N_e3,N_e1];
[t3,f3] = ode45(@(t,y) parodox_Elowitz2_fitnewnew(t,y,x4n,B_ext,lambdasq,I),tspan3,y03');
h1 = plot(t3./24,(f3(:,2)),'LineWidth',8,'Color',[100,255,201]./255);
hold on
h2 = plot([t2(end);t3]./24,([f2(end,3);f3(:,3)]),'LineWidth',8,'Color',[1,0,0]);


xlabel('Time (days)')
ylabel('Cell Population')
xlim([0,75])
ax = gca;
ax.FontSize = 25;
box off
%% _Simple Feedback Response to Mutation_

[t1,f1] = ode45(@(t,y) parodox_Elowitz3_fitnewnew(t,y,x4n,B_ext,lambdasq,I),tspan1,y01');
figure(14)
plot(t1./24,(f1(:,2)),'LineWidth',8,'Color',[100,255,201]./255)
hold on
plot(t1./24,(f1(:,3)),'LineWidth',8,'Color',[1,0,0])

%1st mutation
T2 = 50;
tspan2 = linspace(T1*24,T3*24,10000);
auxin_e2 = f1(end,1);
N_e2 = f1(end,2);


y02 = [auxin_e2,N_e2,N_e1];
[t2,f2] = ode45(@(t,y) parodox_Elowitz3_fitnewnew(t,y,x4n,B_ext,lambdasq,I),tspan2,y02');
plot(t2./24,(f2(:,2)),'LineWidth',8,'Color',[100,255,201]./255)
hold on
plot([t1(end);t2]./24,([f1(end,3);f2(:,3)]),'LineWidth',8,'Color',[1,0,0])



xlabel('Time (days)')
ylabel('Cell Population')
xlim([0,75])
ax = gca;
ax.FontSize = 25;
box off


%% 
% 
%% _Discussion Points_
%% 
% * Do the time scales seem correct?  
%% 
%