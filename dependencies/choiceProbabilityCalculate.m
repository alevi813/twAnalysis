function [m,s] = choiceProbabilityCalculate(X, Y)
% calculate choice probability over range of values
% [m,s] = calculateChoiceProbability(X, Y)
% Inputs:
% 		X [m x n] m trials by n bins spikes
% 		Y [m x 1] m trials by 1 (binary) true classes
% Outputs:
% 		m [n x 1] mean area under the ROC curve for X and Y
% 		s [n x 1] standard error of m

% jly  20140929 wrote it

dims = size(X);
m = zeros(dims(2),1);
s = zeros(dims(2),1);

for ii = 1:dims(2)
	tmp = roc(X(:,ii), Y);
	m(ii) = tmp.AUC;
	s(ii) = tmp.serror;
end
