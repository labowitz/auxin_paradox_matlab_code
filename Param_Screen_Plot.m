load('Data\ParamMatrix.mat')

allExp = [ratePure_Death; ratePure_Growth; ratemono_stable; ratemono_unstable; ratebi_stable];

close all;
Size = 8.5;

%%
figure(1);

colorB = [.2, 0, 1];
colorG = [.2, 1, .2];
colorR = [.9, 0, 0];
colorY = [.9, .9, .1];
colorP = [.1, .9, 1];

for a = 1:9
    subplot(3, 3, a, 'replace');
    hold on;

    plot(BLAST_PD{a},AP1903_PD{a},'o','MarkerSize',Size,'MarkerFaceColor',colorR,'MarkerEdgeColor',colorR)
    plot(BLAST_PG{a},AP1903_PG{a},'o','MarkerSize',Size,'MarkerFaceColor',colorG,'MarkerEdgeColor',colorG)   
    
    plot(BLAST_MU{a},AP1903_MU{a},'o','MarkerSize',Size,'MarkerFaceColor',colorB,'MarkerEdgeColor',colorB)    
    plot(BLAST_MS{a},AP1903_MS{a},'o','MarkerSize',Size,'MarkerFaceColor',colorY,'MarkerEdgeColor',colorY)
    
    plot(BLAST_BIS{a},AP1903_BIS{a},'o','MarkerSize',Size,'MarkerFaceColor',colorP,'MarkerEdgeColor',colorP)
 
    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';
    
    hold off;
end


%%
auxins = logspace(-2, 2, 100);

rVa_G = double(subs(ratePure_Growth(10),'A',auxins));

figure(2)
plot(auxins, rVa_G, 'r-', 'LineWidth', 3, 'Color', colorG);

xlabel('NAA (uM)');
ylabel({'Net Growth Rate (hr^{-1})'})

box off;

ylim([0 ,.05])
xticks([.01, .1, 1, 10, 100])
xtickformat('%.1e')

ax = gca;
ax.XScale = 'log';
ax.FontSize = 18;
ax.XAxisLocation = 'origin';


%%
rVa_PD = double(subs(ratePure_Death(10),'A',auxins));

figure(3)
plot(auxins, rVa_PD, '-', 'LineWidth', 3, 'Color', colorR);

xlabel('NAA (uM)');
ylabel({'Net Growth Rate (hr^{-1})'})

box off;

ylim([-0.06, .0])
xticks([.01, .1, 1, 10, 100])
xtickformat('%.1e')

ax = gca;
ax.XScale = 'log';
ax.FontSize = 18;
ax.XAxisLocation = 'origin';

%%
rVa_MU = double(subs(ratemono_unstable(80),'A',auxins));

figure(4)
plot(auxins, rVa_MU, '-', 'LineWidth', 3, 'Color', colorB);

xlabel('NAA (uM)');
ylabel({'Net Growth Rate (hr^{-1})'})

box off;

ylim([-0.04, .04])
xticks([.01, .1, 1, 10, 100])
xtickformat('%.1e')

ax = gca;
ax.XScale = 'log';
ax.FontSize = 18;
ax.XAxisLocation = 'origin';

%%
rVa_MS = double(subs(ratemono_stable(80),'A',auxins));

figure(5)
plot(auxins, rVa_MS, '-', 'LineWidth', 3, 'Color', colorY);

xlabel('NAA (uM)');
ylabel({'Net Growth Rate (hr^{-1})'})

box off;

ylim([-0.1, .04])
xticks([.01, .1, 1, 10, 100])
xtickformat('%.1e')

ax = gca;
ax.XScale = 'log';
ax.FontSize = 18;
ax.XAxisLocation = 'origin';
rVa_BS = double(subs(ratebi_stable(30),'A',auxins));

%%
figure(6)
plot(auxins, rVa_BS, '-', 'LineWidth', 3, 'Color', colorP);

xlabel('NAA (uM)');
ylabel({'Net Growth Rate (hr^{-1})'})

box off;

ylim([-0.06, .04])
xticks([.01, .1, 1, 10, 100])
xtickformat('%.1e')

ax = gca;
ax.XScale = 'log';
ax.FontSize = 18;
ax.XAxisLocation = 'origin';