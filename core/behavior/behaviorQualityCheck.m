eList   = [58 60 61 63 64 65 66 67 68 69 70 71 72 74 75 76 77 78 79 80 81 82 83 85 86 87 88 89];


for ii = 1:length(eList)
    temp = ppLoader('nancy', eList(ii));
    
    slope(ii)      = temp{1,1}.pp.pmfUnfolded.theta(2);
    inflection(ii) = temp{1,1}.pp.pmfUnfolded.theta(1);
    
end

figure
subplot(1,2,1); hold on
plot(1:length(eList), slope,'-o')
xlabel('Session')
ylabel('slope of PMF')
line([10.5 10.5], [0.04 0.13])
line([20.5 20.5], [0.04 0.13])

subplot(1,2,2)
plot(1:length(eList), inflection, '-o')
xlabel('Session')
ylabel('Inflection pt of PMF')

line([10.5 10.5], [-25 10])
line([20.5 20.5], [-25 10])