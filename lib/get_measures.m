% -------------------------------------------------------------------------
% GET_MEASURES calculates the measures using six different thresholding
% methods:
%
%   1   connectivitiy threshold on IW
%   2   elbow rule on IW
%   3   ES > 0.5
%   4   elbow rule on IW and ERS > 0.5
%   5   ES > 0.5 and ERS > 0.5
%   6   IW with same number of edges as (3)
%   7   ERS > 0.5
% -------------------------------------------------------------------------

function M = get_measures(A, edges, subset)

A_true = A(subset,subset);

% Rank order IW and get elbow.
[~, inds] = sort(edges(:,3), 'descend');
thresh = get_elbow(edges(inds,3)');

% Use IW as filter to get connectivity of 0.05
N = length(subset);
n = round(0.05*N*N);
inds_N = inds(1:n);
B_N = make_matrix(A, edges, inds_N);
A_pred_N = B_N(subset,subset);
M.N = get_measure(A_true, A_pred_N);

% Use IW as filter.
inds_IW = inds(1:thresh);
B_IW = make_matrix(A, edges, inds_IW);
A_pred_IW = B_IW(subset,subset);
M.IW = get_measure(A_true, A_pred_IW);

% Use ES as filter.
inds_ES = find(edges(:,4) > 0.5);
B_ES = make_matrix(A, edges, inds_ES);
A_pred_ES = B_ES(subset,subset);
M.ES = get_measure(A_true, A_pred_ES);

% Use IW as filter with ERS.
inds_IW_ERS = inds_IW(edges(inds_IW, 5) > 0.5);
B_IW_ERS = make_matrix(A, edges, inds_IW_ERS);
A_pred_IW_ERS = B_IW_ERS(subset,subset);
M.IW_ERS = get_measure(A_true, A_pred_IW_ERS);

% Use ES as filter with ERS.
inds_ES_ERS = find(edges(:,4) > 0.5 & edges(:,5) > 0.5);
B_ES_ERS = make_matrix(A, edges, inds_ES_ERS);
A_pred_ES_ERS = B_ES_ERS(subset,subset);
M.ES_ERS = get_measure(A_true, A_pred_ES_ERS);

% Use ERS alone.
inds_ERS = find(edges(:,5) > 0.5);
B_ERS = make_matrix(A, edges, inds_ERS);
A_pred_ERS = B_ERS(subset,subset);
M.ERS = get_measure(A_true, A_pred_ERS);

% Use IW with same number as ES.
n_IWes = length(inds_ES);
inds_IWes = inds(1:n_IWes);
B_IWes = make_matrix(A, edges, inds_IWes);
A_pred_IWes = B_IWes(subset,subset);
M.IWes = get_measure(A_true, A_pred_IWes);

% Use IW with same number as ERS.
n_IWers = length(inds_ERS);
inds_IWers = inds(1:n_IWers);
B_IWers = make_matrix(A, edges, inds_IWers);
A_pred_IWers = B_IWers(subset,subset);
M.IWers = get_measure(A_true, A_pred_IWers);

% Use IW with same number as IW with ERS.
n_IWiwers = length(inds_IW_ERS);
inds_IWiwers = inds(1:n_IWiwers);
B_IWiwers = make_matrix(A, edges, inds_IWiwers);
A_pred_IWiwers = B_IWiwers(subset,subset);
M.IWiwers = get_measure(A_true, A_pred_IWiwers);

% Use IW with same number as ES with ERS
n_IWesers = length(inds_ES_ERS);
inds_IWesers = inds(1:n_IWesers);
B_IWesers = make_matrix(A, edges, inds_IWesers);
A_pred_IWesers = B_IWesers(subset,subset);
M.IWesers = get_measure(A_true, A_pred_IWesers);

end

function B = make_matrix(A, edges, inds)
    B = zeros(size(A));
    for i = inds'
        B(edges(i,1), edges(i,2)) = 1;
    end
end