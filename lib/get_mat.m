% -------------------------------------------------------------------------
% GET_MAT returns an nP x nP consolidated matrix of values for nP
% parameters pulled from the given data object.
% -------------------------------------------------------------------------

function A = get_mat(D, from, to, nP, iSlice)

A = zeros(nP);

for iParamA = 1:nP
    for iParamB = 1:nP
        A(iParamA, iParamB) = D{iParamA, iParamB, iSlice}(from, to);
    end
end

end