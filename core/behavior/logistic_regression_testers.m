
%%

for iSubject = 1:length(bSubject)
    
    num_trials(iSubject) = length(bSubject(iSubject).stimulus.pulses_zScore);
    num_pulses = 7;
    revco_trials = bSubject(iSubject).reverseCorrelation;
    
    pulses = bSubject(iSubject).stimulus.pulses_zScore;
    
    % figure(1); clf
    % histogram(pulses(revco_trials,:)); hold on
    % histogram(pulses(~revco_trials,:));
    
    kernel_true = bSubject(iSubject).ppkRidge.w_norm;
    
    sigmoid = @(x) 1./(1 + exp(-x));
    dv = pulses*kernel_true;
    pRight = sigmoid(dv);
    choices = rand(num_trials(iSubject),1) < pRight;
    
    [kernel_revco, ~, s_revco] = glmfit(pulses(revco_trials,:), choices(revco_trials), 'binomial');
    [kernel_all, ~, s_all] = glmfit(pulses, choices, 'binomial');
    
    kernel_rVal(iSubject) = corrcoef(kernel_all(2:end), kernel_revco(2:end));
    
    num_revco_trials(iSubject) = length(pulses(revco_trials,:));
    
    figure(1);
    if iSubject <= 4
        subplot(3, 5, iSubject)
    else
        subplot(3, 5, iSubject+1)
    end
    plot(kernel_true, 'k'); hold on
    errorbar(1:num_pulses, kernel_revco(2:end), s_revco.se(2:end), 'b')
    errorbar(1:num_pulses, kernel_all(2:end), s_all.se(2:end), 'r')
    xlim([.5 7.5])
    set(gca, 'XTick', [ 1 2 3 4 5 6 7 ] );
    xlabel('Pulse #')
    ylabel('Weight')
end

legend({'True', 'Zero-mean', 'All'})

set(gcf, 'color', 'w');






