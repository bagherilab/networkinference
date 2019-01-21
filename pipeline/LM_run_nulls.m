% -------------------------------------------------------------------------
% LM_RUN_NULLS runs selected algorithms for the null models and
% saves the results to a structure of cell arrays.
%
% This function is designed specifically for use with a job array for
% external parallelization for the (nN noise types) x (nS stimulus types)
% x (nP parameter A values) as well as internal Matlab built-in parfor
% parallelization for the (nP parameter B values) across (N nulls) x
% (nT timeslices) x (nR replicates).
% -------------------------------------------------------------------------

function [output, code] = LM_run_nulls(INDEX, ALGORITHM, HIDDEN, INPUT)

SETTINGS = LM_SETTINGS();
nP = SETTINGS.nParams; % number of parameters
nR = SETTINGS.nReplicates; % number of replicates
nT = SETTINGS.nTimeslices; % number of slices
N = SETTINGS.nNulls; % number of null models

% Get suffix from selected INPUT.
suffix = SETTINGS.inputNames{INPUT};

% If including HIDDEN, only the full timeslice is used.
if HIDDEN
    nT = 1;
end

[code, iNoise, iStim, iParamA] = get_code(INDEX, 2); % convert array index
jarpath = SETTINGS.banjoJar; % path to BANJO jar file

% Setup parallelization.
results = cell(nP, N);
warning('off');

% Load data.
D = load(['Simulations_NULL_S' num2str(iStim) suffix '.mat']);
noise = SETTINGS.noiseNames{iNoise};
D = D.(noise);

% Run algorithm.
parfor iParamB = 1:nP
    LOCAL = LM_SETTINGS();
    algorithm = get_algorithm(ALGORITHM, LOCAL);
    d = D{iParamA, iParamB};

    for i = 1:N
        fprintf('Running [%s] on NOISE[%d] STIM[%d] PA[%d] PB[%d] for [%d]\n', ...
            ALGORITHM, iNoise, iStim, iParamA, iParamB, i);
        weights = cell(nT, nR);
        
        for iTime = 1:nT
            slice = LOCAL.slices{iTime};
            dT = d{i}(:, slice);
            
            if HIDDEN
                stimFun = get_stim_function(INPUT, LOCAL.tStim, LOCAL.stimAmount);
                dT = [dT; stimFun(LOCAL.timePoints)];
            end
            
            for iRep = 1:nR
                if strcmp(ALGORITHM, 'BANJO')
                    javaaddpath(jarpath);
                end
                subcode = [code '_' sprintf('B%dT%dR%d%s-%d', ...
                    iParamB, iTime, iRep, suffix, i)];
                weights{iTime, iRep} = algorithm(dT, subcode);
            end
        end
        
        % Store raw weights.
        results{iParamB, i} = weights;
    end
end

output.(noise) = results;

end