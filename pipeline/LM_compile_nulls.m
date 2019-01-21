% -------------------------------------------------------------------------
% LM_COMPILE_NULLS compiles results for a given algorithm across the nP
% parameter A subruns into a single cell matrix for a given noise type
% and stimulus type.
% -------------------------------------------------------------------------

function [output, code] = LM_compile_nulls(INDEX, ALGORITHM, SUFFIX)

SETTINGS = LM_SETTINGS();
nP = SETTINGS.nParams; % number of parameters
[code, iNoise, ~, ~] = get_code(INDEX, 2); % convert array index
path = ['Results_' ALGORITHM '_'];
noise = SETTINGS.noiseNames{iNoise};
results = cell(nP);

for iParamA = 1:nP
    d = load([path code(1:4) 'A' num2str(iParamA) SUFFIX '.mat']);
    
    for iParamB = 1:nP
        results{iParamA, iParamB} = d.(noise)(iParamB, :);
    end
end

output.(noise) = results;

end