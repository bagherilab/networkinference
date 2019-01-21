% -------------------------------------------------------------------------
% LM_CALC_METRICS calculates edge score (ES) and edge rank score (ERS)
% for parameter pairs of given motif/logic gate/stimulus.
%
% The function iterates through each pair of parameter values and
% calculates the ES and ERS for the true against all nulls models. The
% average and standard deviation across replicates are calculated as well
% as the mean and standard deviation of the raw weights.
% -------------------------------------------------------------------------

function [output, code] = LM_calc_metrics(INDEX, ALGORITHM, HIDDEN, SUFFIX)

SETTINGS = LM_SETTINGS();
nN = SETTINGS.nNoises; % number of noise variations
nP = SETTINGS.nParams; % number of parameters
nT = SETTINGS.nTimeslices; % number of slices
nR = SETTINGS.nReplicates; % number of replicates
N = SETTINGS.nNulls; % number of null models
[code, iMotif, iLogic, iStim] = get_code(INDEX, 1);
path = ['Results_' ALGORITHM];
parser = get_parser(ALGORITHM);

% If including HIDDEN, only the full timeslice is used.
if HIDDEN
    nT = 1;
end

% Check for special cases.
nanCodeCheck = sum(~cellfun('isempty', strfind(SETTINGS.nanCodes, code)));
subsetCodeCheck = sum(~cellfun('isempty', strfind(SETTINGS.subsetCodes, code)));

if nanCodeCheck == 1
    fprintf('*** NAN CASE FOR [%s] on MOTIF[%d] GATE[%d] S[%d] ***\n', ...
        ALGORITHM, iMotif, iLogic, iStim);
        output = SETTINGS.empty2;
    return;
end

% Setup parallelization.
analysis = cell(nN);

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

parfor iNoise = 1:nN
    LOCAL = LM_SETTINGS();
    noise = LOCAL.noiseNames{iNoise};
    
    % Load data.
    trueData = load([path '_' code SUFFIX '.mat']);
    nullData = load([path '_N' num2str(iNoise) code(5:6) SUFFIX '.mat']);
    dTrue = trueData.(noise);
    dNull = nullData.(noise);
    
    % Initialize output container.
    s = struct();

    for iParamA = 1:nP
        for iParamB = 1:nP
            dTrueAB = dTrue{iParamA, iParamB};
            dNullAB = dNull{iParamA, iParamB};

            for iTime = 1:nT
                fprintf('Analyzing [%s] results for PA[%d] PB[%d] T[%d]\n', ...
                    noise, iParamA, iParamB, iTime);
                ES = zeros(5,5,nR);
                ERS = zeros(5,5,nR);

                for iRep = 1:nR
                    % Containers for metrics.
                    es = -Inf*ones(5); es(subset, subset) = 0;
                    ers = -Inf*ones(5); ers(subset, subset) = 0;
                    
                    % Get true weights and ranking.
                    [trueWeights, trueRanks] = get_vects(dTrueAB{iTime, iRep}, subset, parser);
                    
                    for i = 1:N
                        [nullWeights, nullRanks] = get_vects(dNullAB{i}{iTime, iRep}, subset, parser);
                        es(subset, subset) = es(subset, subset) + ...
                            (trueWeights > nullWeights) + 0.5*(trueWeights == nullWeights);
                        ers(subset, subset) = ers(subset, subset) + ...
                            (trueRanks < nullRanks) + 0.5*(trueRanks == nullRanks);
                    end

                    ES(:,:,iRep) = es/N;
                    ERS(:,:,iRep) = ers/N;
                end

                % Calculate mean and standard deviation of true weights.
                trueMat = cat(3, dTrue{iParamA, iParamB}{iTime, :});
                trueMat(subset, subset, :) = parser(trueMat(subset, subset, :));
                s.IW.mean{iParamA, iParamB, iTime} = mean(trueMat, 3);
                s.IW.stdev{iParamA, iParamB, iTime} = std(trueMat, 1, 3);
                
                % Iterate through nulls.
                nullMeans = cell(N,1);
                nullStdevs = cell(N,1);
                for i = 1:N
                    nullMat = parser(cat(3, dNull{iParamA, iParamB}{i}{iTime, :}));
                    nullMat(subset, subset, :) = parser(nullMat(subset, subset, :));
                    nullMat(isinf(trueMat)) = -Inf;
                    
                    % Fix for NaN cases in correlation of all zeros.
                    if strcmp(ALGORITHM, 'CORR')
                        nullMat(isnan(nullMat)) = 0;
                    end
                    
                    nullMeans{i} = mean(nullMat, 3);
                    nullStdevs{i} = std(nullMat, 1, 3);
                end
                
                % Calculate mean and standard deviation of nulls across replicates.
                s.NW.mean{iParamA, iParamB, iTime} = mean(cat(3, nullMeans{:}), 3);
                s.NW.stdev{iParamA, iParamB, iTime} = std(cat(3, nullStdevs{:}), 1, 3);
        
                % Store values into structure.
                s.ES.mean{iParamA, iParamB, iTime} = mean(ES, 3);
                s.ES.stdev{iParamA, iParamB, iTime} = std(ES, 1, 3);
                s.ERS.mean{iParamA, iParamB, iTime} = mean(ERS, 3);
                s.ERS.stdev{iParamA, iParamB, iTime} = std(ERS, 1, 3);
            end
        end
    end
    
    analysis{iNoise} = s;
end

% Organize simulations by noise.
for iNoise = 1:nN
    noise = SETTINGS.noiseNames{iNoise};
    output.(noise) = analysis{iNoise};
end

end

function [weights, ranks] = get_vects(A, subset, parser)
    n = length(subset);
    weight_vect = reshape(parser(A(subset, subset)), [], 1);
    [~, inds] = sort(weight_vect, 'descend');
    rank_vect(inds) = 1:(n*n);
    ranks = reshape(rank_vect, n, []);
    weights = reshape(weight_vect, n, []);
end