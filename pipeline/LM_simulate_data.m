% -------------------------------------------------------------------------
% LM_SIMULATE_DATA simulates a five node ODE system for every combination
% of motif structure, logic gate, parameter A value, parameter B value,
% and stimulus type.
%
% This function is designed specifically for use with a job array for
% external parallelization for the (nM motifs) x (nL logic gates) as well 
% as internal Matlab built-in parfor parallelization for the 
% (nP parameter A values) x (nP parameter B values) x (nS stimulus types).
% -------------------------------------------------------------------------

function [output, code] = LM_simulate_data(INDEX, INPUT)

SETTINGS = LM_SETTINGS();
nN = SETTINGS.nNoises; % number of noise variations
nP = SETTINGS.nParams; % number of parameters
nS = SETTINGS.nStims; % number of stimulation types
[code, iMotif, iLogic, ~] = get_code(INDEX, 1); % convert array index

% Select motifs and gates.
motif = SETTINGS.motifMatrices{iMotif};
logic = SETTINGS.logicGates{iLogic};

% Setup parallelization.
simulations = cell(nN, nP, nP, nS);

parfor iParamA = 1:nP
    % Local copy of settings to reduce overhead.
    LOCAL = LM_SETTINGS(); 
    timeSpan = LOCAL.timeSpan;
    stim = LOCAL.stimVector;
    paramA = LOCAL.paramVals(iParamA); % get parameter A value
    container = cell(nN, nP, nS);
    
    for iParamB = 1:nP
        fprintf('Simulating MOTIF[%d] GATE[%d] PA[%d] PB[%d]\n', ...
            iMotif, iLogic, iParamA, iParamB);
        paramB = LOCAL.paramVals(iParamB); % get parameter B value

        % Get activation and stimulation functions.
        gateFun = get_gate_function(logic, paramA, paramB);
        hillFun = get_hill_function(LOCAL.k0);
        stimFun = get_stim_function(INPUT, LOCAL.tStim, LOCAL.stimAmount);

        % Simulate data to steady state to get initial conditions.
        odeFunSteady = @(t, y) five_node_network(t, y, motif, [0 0 0 0 0], ...
            LOCAL.decayConstants, gateFun, hillFun, stimFun);
        solSteady = ode45(odeFunSteady, timeSpan*10, [0 0 0 0 0], LOCAL.solverOptions);
        initial = solSteady.y(:,end);
        
        % Simulate data from steady state with different stimulation types.
        for iStim = 1:nS
            odeFun = @(t, y) five_node_network(t, y, motif, stim{iStim}, ...
                LOCAL.decayConstants, gateFun, hillFun, stimFun);
            sol = ode45(odeFun, timeSpan, initial, LOCAL.solverOptions);
            
            % Sample from data at selected time points.
            sample = deval(sol, LOCAL.timePoints);
            
            % Add noise.
            for iNoise = 1:nN
                rng((INDEX - 1)*nP*nP*nS + (iParamA - 1)*nP*nS + ...
                    (iParamB - 1)*nS + iStim);
                perc = LOCAL.noisePercents(iNoise);
                container{iNoise, iParamB, iStim} = add_noise(sample, perc);
            end
        end
    end
    
    simulations(:, iParamA, :, :) = container;
end

% Organize simulations by noise.
for iNoise = 1:nN
    noise = SETTINGS.noiseNames{iNoise};
    output.(noise) = cell(nP,nP,nS);
    output.(noise)(:,:,:) = simulations(iNoise,:,:,:);
end

end