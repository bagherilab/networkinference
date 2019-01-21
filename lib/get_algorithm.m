% -------------------------------------------------------------------------
% GET_ALGORITHM returns a function handle to the selected algorithm
% parameterized by the given settings.
% -------------------------------------------------------------------------

function F = get_algorithm(ALGORITHM, SETTINGS)

switch ALGORITHM
    case 'CORR'
        F = @(d, id) corrcoef(d');
    case 'GENIE3'
        F = @(d, id) GENIE3(d', SETTINGS.regulatorIDs, ...
            SETTINGS.treeMethod, SETTINGS.K);
    case 'MIDER'
        F = @(d, id) mider(d', SETTINGS.miderOptions);
    case 'TIGRESS'
        F = @(d, id) score_edges(tigress(struct('expdata', d'), ...
            'R', SETTINGS.resampling, 'verbose', false, ...
            'alpha', SETTINGS.randomization, 'L', SETTINGS.lars));
    case 'BANJO'
        F = @(d, id) wrap_banjo_jar(d', id, ...
            SETTINGS.banjoTemplate, SETTINGS.filePath);
end

end