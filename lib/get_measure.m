% -------------------------------------------------------------------------
% GET_MEASURE calculates the number of true positives, false positives,
% true negatives, and false negatives between two matricies excluding
% self edges.
% -------------------------------------------------------------------------

function M = get_measure(A_true, A_pred)

% Identify cases.
tp = (A_pred == 1) & (A_true == 1);
fp = (A_pred == 1) & (A_true == 0);
tn = (A_pred == 0) & (A_true == 0);
fn = (A_pred == 0) & (A_true == 1);

A = nan(size(A_true));

% Replace with case.
A(tp) = 1;
A(fp) = 2;
A(tn) = 3;
A(fn) = 4;

% Remove self edges.
B = triu(A, 1) + tril(A, -1);
B = reshape(B, 1, []);

% Count number of each case.
TP = sum(B == 1);
FP = sum(B == 2);
TN = sum(B == 3);
FN = sum(B == 4);

N = size(A_true,1);
n = N*(N -1);
M = [TP/n FP/n TN/n FN/n];

end
