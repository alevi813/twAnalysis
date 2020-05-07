function [] = twagSummaryPlot(pp, cp, maxPTAm)


figure
subplot(2,3,1);
errorbar(1:7, pp.kernel.b, pp.kernel.s.se, ...
    'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
title('Psychophysical kernel')
xlabel('Pulse')
ylabel('Weight')
set(gca, 'Xtick', 1:7)
xlim([.5 7.5]);

subplot(2,3,2);
errorbar(1:7, cp.session.allRevco.pulse.m, cp.session.allRevco.pulse.se, ...
    'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
title('Choice probability by Pulse')
h = refline(0,0.5);
set(h, 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Pulse')
ylabel('CP')
set(gca, 'Xtick', 1:7)
xlim([.5 7.5]);

subplot(2,3,3);
plot(1:7, maxPTAm,...
    'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
title('Max PTA value by pulse')
xlabel('Pulse')
ylabel('Pta max value')
set(gca, 'Xtick', 1:7)
xlim([.5 7.5]);

% plot them against each other
subplot(2,3,4);
scatter(pp.kernel.b, cp.session.allRevco.pulse.m);
xlabel('Behavioral weight');
ylabel('Pulse CP');
title('PPK v. CP')

% plot them against each other
subplot(2,3,5);
scatter(pp.kernel.b, maxPTAm);
xlabel('Behavioral weight');
ylabel('PTA max value');
title('PPK v. PTA')

% plot them against each other
subplot(2,3,6);
scatter(cp.session.allRevco.pulse.m, maxPTAm);
xlabel('Pulse CP');
ylabel('PTA max value');
title('CP v. PTA')

%% summary plot using frozen only
%plot kernel and pulse cp
figure
subplot(2,3,1);
errorbar(1:7, pp.kernel.b, pp.kernel.s.se, ...
    'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
title('Psychophysical kernel')
xlabel('Pulse')
ylabel('Weight')
set(gca, 'Xtick', 1:7)
xlim([.5 7.5]);

subplot(2,3,2);
errorbar(1:7, cp.session.frozen.pulse.m, cp.session.frozen.pulse.se, ...
    'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
title('Choice probability by Pulse')
h = refline(0,0.5);
set(h, 'Color', [0 0 0], 'LineStyle', '--');
xlabel('Pulse')
ylabel('CP')
set(gca, 'Xtick', 1:7)
xlim([.5 7.5]);

subplot(2,3,3);
plot(1:7, maxPTAm,...
    'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
title('Max PTA value by pulse')
xlabel('Pulse')
ylabel('Pta max value')
set(gca, 'Xtick', 1:7)
xlim([.5 7.5]);

% plot them against each other
subplot(2,3,4);
scatter(pp.kernel.b, cp.session.frozen.pulse.m);
xlabel('Behavioral weight');
ylabel('Pulse CP');
title('PPK v. CP')

% plot them against each other
subplot(2,3,5);
scatter(pp.kernel.b, maxPTAm);
xlabel('Behavioral weight');
ylabel('PTA max value');
title('PPK v. PTA')

% plot them against each other
subplot(2,3,6);
scatter(cp.session.frozen.pulse.m, maxPTAm);
xlabel('Pulse CP');
ylabel('PTA max value');
title('CP v. PTA')