function s = LM_GNW_SETTINGS()


%% OUTPUT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 50;
s.id = 'C';
s.name = ['GNW' s.id num2str(N)];
s.path = '/home/jessica_yu/PRJ_NetworkInference/';
s.filePath = [s.path 'Data/' s.name '/'];
s.banjoTemplate = [s.path 'Code/_banjo/banjo_settings_template.txt'];
s.banjoJar = [s.path 'Code/_banjo/banjo.jar'];

%% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.timePoints = 1:21;
s.nTimepoints = length(s.timePoints);

%% NETWORK INFERENCE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if s.id == 'C'
    s.replicates = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J' };
else
    s.replicates = { 'A', 'B', 'C' };
end

s.nReplicates = length(s.replicates);
s.nNodes = N;
s.nNulls = 100;

s.nIntervals = 3;
s.intervals = { 1:21, 1:11, 11:21 };

%% GENIE3 PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.regulatorIDs = 1:s.nNodes;
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