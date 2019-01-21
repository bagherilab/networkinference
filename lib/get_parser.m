% -------------------------------------------------------------------------
% GET_PARSER returns a function handle to the selected algorithm-specific
% parser given a nxn matrix of inferred weights for n nodes.
% -------------------------------------------------------------------------

function F = get_parser(ALGORITHM)

switch ALGORITHM
    case 'CORR'
        F = @(mat) abs(mat);
    case 'GENIE3'
        F = @(mat) threshold(mat);
    case 'MIDER'
        F = @(mat) threshold(mat);
    case 'TIGRESS'
        F = @(mat) mat;
    case 'BANJO'
        F = @(mat) mat;
end

end

function B = threshold(A)
    B = A; B(A < 0) = 0;
end