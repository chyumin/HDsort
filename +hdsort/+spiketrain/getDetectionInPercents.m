function [TPR, PPV, ERR, FCR, out] = getDetectionInPercents(R)
% Input: R object that is created in hdsort.spiketrain.alignGDFs()

x = R.table;
idx = x(:, 1) == -1;
x(idx,:) = [];

out.N_true = x(:, 3);
out.N_true_overlaps = x(:, 4);
out.N_true_non_overlaps = x(:, 5);
out.N_detected = x(:, 6);
out.FA = x(:, 7);
out.FP = x(:, 8);
out.FN = x(:, 9);
out.TP = x(:, 12);
out.FC = x(:, 15);
out.E_total = x(:, 21);

TPR = 100*out.TP ./ out.N_true;
TPR(isnan(TPR))=0;

PPV = 100*out.TP ./ out.N_detected;
PPV(isnan(PPV))=0;

ERR = 100*out.E_total ./ out.N_true;
ERR(isnan(ERR)) = 0;

FCR = 100*out.FC ./ out.N_true;
FCR(isnan(FCR))=0;

%  1 'True Unit'
%  2 'assigned Unit'
%  3 'True Spikes'
%  4 'Non-Overlaps'
%  5 'Overlaps'
%  6 'Detected'
%  7 'Falsely Assigned from other Neuron'
%  8 'False Positive'
%  9  'False Negative'
% 10 'False Negative Non-Overlaps'
% 11 'False Negative Overlap'
% 12 'True Positive'
% 13 'True Positive Non-Overlaps'
% 14 'True Positive Overlaps'
% 15 'Classification Errors'
% 16 'Classification Errors Non-Overlaps'
% 17 'Classification Errors Overlaps'
% 18 'Detection Errors Overlaps'
% 19 'total Detection Errors'
% 20 'total Errors Overlaps'
% 21 'total Errors'

end
