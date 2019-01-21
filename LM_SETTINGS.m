function s = LM_SETTINGS()

%% OUTPUT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.path = '/home/jessica_yu/PRJ_NetworkInference/';
s.codePath = [s.path 'Code'];
s.filePath = [s.path 'results/'];
s.banjoTemplate = [s.path '_banjo/banjo_settings_template.txt'];
s.banjoJar = [s.path '_banjo/banjo.jar'];

%% INFERENCE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
s.nEdges = 2;
s.edgeInds = {[1, 3], [2, 3]};
s.edgeNames = {'AC', 'BC'};
s.nValues = 4;
s.values = {'IW', 'NW', 'ES', 'ERS'};
 
%% MODEL SELECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.nParams = 17;
s.paramVals = logspace(-2, 2, s.nParams);
s.nLogics = 6;
s.logicGates = {'AND', 'OR', 'SUM', 'NAND', 'NOR', 'SUB'};
s.nMotifs = 6;
s.motifs = {'FI', 'FF', 'UFB', 'SFB', 'DFB', 'FBFF'};
s.motifMatrices = {
	[0 0 1 0 0; 0 0 1 0 0; 0 0 0 1 1; 0 0 0 0 0; 0 0 0 0 0] % two node fan-in
	[0 1 1 0 0; 0 0 1 0 0; 0 0 0 1 1; 0 0 0 0 0; 0 0 0 0 0] % feed forward
    [0 1 1 0 0; 1 0 1 0 0; 0 0 0 1 1; 0 0 0 0 0; 0 0 0 0 0] % upsteam feedback
	[0 0 1 0 0; 0 0 1 0 0; 0 1 0 1 1; 0 0 0 0 0; 0 0 0 0 0] % single feedback
	[0 0 1 0 0; 0 0 1 0 0; 1 1 0 1 1; 0 0 0 0 0; 0 0 0 0 0] % double feedback
	[0 1 1 0 0; 0 0 1 0 0; 1 0 0 1 1; 0 0 0 0 0; 0 0 0 0 0] % feedback w. feed forward
};

%% ODE SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.solverOptions = odeset('RelTol', 1e-6, 'NonNegative', 1:5);
t0 = 0; % experiment starting time
tF = 10; % experiment ending time
interval = 0.5; % sampling interval
s.tStim = tF/2 + 0.0001; % stimulation end time
s.timeSpan = [t0 tF];
s.timePoints = t0:interval:tF;

%% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.nNoises = 5;
s.noiseNames = {'noise_000', 'noise_005', 'noise_010', 'noise_020', 'noise_050'};
s.noisePercents = [0, 0.05, 0.1, 0.2, 0.5];
s.k0 = 0.5; % unitless rate constant for non-fan-in edges
s.decayConstants = [0.5 0.5 0.5 0.5 0.5];
s.stimAmount = 1;
s.nStims = 3;
s.stimTypes = {'AB', 'A', 'B'};
s.stimVector = {
    [1 1 0 0 0]
    [1 0 0 0 0]
    [0 1 0 0 0]
};

%% STIMULUS INPUT TYPES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.nInputs = 7;
s.inputNames = {'', '_RAMP_UP', '_RAMP_DOWN', '_STEP_UP', '_STEP_DOWN', ...
    '_PULSE_TWO', '_PULSE_THREE' };

%% NETWORK INFERENCE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.algorithms = {'GENIE3', 'TIGRESS', 'BANJO', 'CORR', 'MIDER'};
s.nAlgorithms = 5;
s.nReplicates = 3;
s.nTimeslices = 3;
s.nNulls = 100;

s.slices = {
    find(s.timePoints >= 0)
    find(s.timePoints <= tF/2)
    find(s.timePoints >= tF/2)
};

s.nanCodes = {
    'M1L1S2', 'M4L1S2', 'M5L1S2', 'M2L1S3', ...
    'M1L1S3', 'M4L1S3', 'M5L1S3', 'M6L1S3' };
s.subsetCodes = {
    'M1L2S2', 'M1L3S2', 'M1L4S2', 'M1L5S2', 'M1L6S2', ...
    'M1L2S3', 'M1L3S3', 'M1L4S3', 'M1L5S3', 'M1L6S3', ...
    'M2L2S3', 'M2L3S3', 'M2L4S3', 'M2L5S3', 'M2L6S3', ...
    'M4L2S3', 'M4L3S3', 'M4L4S3', 'M4L5S3', 'M4L6S3' };

% Build matrix of nans for no network cases.
nans5 = {nan(5, 5)};
weights5 = cell(s.nTimeslices, s.nReplicates);
weights5(:) = nans5;
empty5 = cell(s.nParams);
empty5(:) = {weights5};

nans6 = {nan(6, 6)};
weights6 = cell(s.nTimeslices, s.nReplicates);
weights6(:) = nans6;
empty6 = cell(s.nParams);
empty6(:) = {weights6};

params = cell(s.nParams, s.nParams, s.nTimeslices);
params(:) = nans5;
stats.mean = params;
stats.stdev = params;

for i = 1:s.nNoises
    noise = s.noiseNames{i};
    s.empty5.(noise) = empty5;
    s.empty6.(noise) = empty6;
    
    s.empty2.(noise).IW = stats;
    s.empty2.(noise).NW = stats;
    s.empty2.(noise).ES = stats;
    s.empty2.(noise).ERS = stats;
end

%% GENIE3 PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.regulatorIDs = 1:5;
s.treeMethod = 'RF';
s.K = 'all';

%% TIGRESS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.resampling = 1000;
s.randomization = 0.3;
s.lars = 3;

%% MIDER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mider.useStatistics = 0;
mider.correctOutliers = 0;
mider.q = 1;
mider.MItype = 'MI';
mider.fraction = max(0.1*(log10(length(s.timePoints)) - 1), 0.01);
mider.taumax = 3;
mider.ert_crit = 2;
mider.threshold = 1;
mider.plotMI = 0;
s.miderOptions = mider;

end