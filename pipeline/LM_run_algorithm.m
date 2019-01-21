% -------------------------------------------------------------------------
% LM_RUN_ALGORITHM runs selected algorithms for the true model and
% saves results to a structure of cell arrays.
%
% This function is designed specifically for use with a job array for
% external parallelization for the (nM motifs) x (nL logic gates) x
% (nS stimulus types) as well as internal Matlab built-in parfor
% parallelization for the (nP parameter A values) across (nP parameter
% B values) x (nN noise types) x (nT timeslices) x (nR replicates).
%
% The HIDDEN flag is used to indicate if data for a hidden "stimulus" node
% should be included in the inference.
% -------------------------------------------------------------------------

function [output, code] = LM_run_algorithm(INDEX, ALGORITHM, HIDDEN, INPUT)

SETTINGS = LM_SETTINGS();
nP = SETTINGS.nParams; % number of parameters
nR = SETTINGS.nReplicates; % number of replicates
nT = SETTINGS.nTimeslices; % number of slices
nN = SETTINGS.nNoises; % number of noise variations

% Get suffix from selected INPUT.
suffix = SETTINGS.inputNames{INPUT};

% If including HIDDEN, only the full timeslice is used.
if HIDDEN
    nT = 1; nodes = 6;
else
    nodes = 5;
end

[code, iMotif, iLogic, iStim] = get_code(INDEX, 1); % convert array index
jarpath = SETTINGS.banjoJar; % path to BANJO jar file

% Check for special cases.
nanCodeCheck = sum(~cellfun('isempty', strfind(SETTINGS.nanCodes, code)));
subsetCodeCheck = sum(~cellfun('isempty', strfind(SETTINGS.subsetCodes, code)));

if nanCodeCheck == 1
    fprintf('*** NAN CASE FOR [%s] on MOTIF[%d] GATE[%d] S[%d] ***\n', ...
        ALGORITHM, iMotif, iLogic, iStim);
        if HIDDEN
            output = SETTINGS.empty6;
        else
            output = SETTINGS.empty5;
        end
    return;
end

% Setup parallelization.
results = cell(nP);
warning('off');

% Adjust data for subset cases.
if subsetCodeCheck == 1
    if iStim == 2
        subset = [1 3 4 5];
    elseif iStim == 3
        subset = [2 3 4 5];
    end
else
    subset = 1:5;
end

% Run algorithm.
parfor iParamA = 1:nP
    LOCAL = LM_SETTINGS();
    
    if subsetCodeCheck == 1
        LOCAL.regulatorIDs = 1:(nodes - 1);
    end

    D = load(['Simulations_FULL_S' num2str(iStim) suffix '.mat']);
    algorithm = get_algorithm(ALGORITHM, LOCAL);

    for iParamB = 1:nP
        fprintf('Running [%s] on MOTIF[%d] GATE[%d] STIM[%d] PA[%d] PB[%d]\n', ...
            ALGORITHM, iMotif, iLogic, iStim, iParamA, iParamB);
        weights = cell(nT, nR, nN);
        
        for iNoise = 1:nN
            noise = LOCAL.noiseNames{iNoise};
            d = D.(noise){iMotif, iLogic}{iParamA, iParamB};
        
            for iTime = 1:nT
                slice = LOCAL.slices{iTime};
                dT = d(subset, slice);
                
                if HIDDEN
                    stimFun = get_stim_function(INPUT, LOCAL.tStim, LOCAL.stimAmount);
                    dT = [dT; stimFun(LOCAL.timePoints)];
                    inds = [subset 6];
                else
                    inds = subset;
                end
                
                for iRep = 1:nR
                    if strcmp(ALGORITHM, 'BANJO')
                        javaaddpath(jarpath);
                    end

                    subcode = [code '_' sprintf('A%dB%dT%dR%dN%d%s', ...
                        iParamA, iParamB, iTime, iRep, iNoise, suffix)];
                    w = -Inf*ones(nodes);
                    w(inds,inds) = algorithm(dT, subcode);
                    weights{iTime, iRep, iNoise} = w;
                end
            end
        end

        results{iParamA, iParamB} = weights;
    end
end

% Organize simulations by noise.
for iNoise = 1:nN
    noise = SETTINGS.noiseNames{iNoise};
    output.(noise) = cell(nP);
    
    for iParamA = 1:nP
        for iParamB = 1:nP
            output.(noise){iParamA, iParamB} = results{iParamA, iParamB}(:,:,iNoise);
        end
    end
end

end