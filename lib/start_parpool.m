% -------------------------------------------------------------------------
% START_PARPOOL starts a parallel pool of the requested number of workers
% and creates a temporary directory for project files.
% -------------------------------------------------------------------------

function dir = start_parpool(SETTINGS, nCores)

% Create a temporary directory for project files.
dir = [SETTINGS.path tempname()];
mkdir(dir);

% Create cluster.
c = parallel.cluster.Local('JobStorageLocation',dir,'NumWorkers',nCores);

% Open parallel pool.
parpool(c, nCores);

end