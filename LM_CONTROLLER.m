% -------------------------------------------------------------------------
% LM_CONTROLLER acts as the entry point from which all other code can be
% run based on TASK ids.
% -------------------------------------------------------------------------

function LM_CONTROLLER(TASK, INDEX, ALGORITHM)

SETTINGS = LM_SETTINGS();
path = SETTINGS.filePath;

nL = SETTINGS.nLogics; % number of logics
nM = SETTINGS.nMotifs; % number of motifs
nP = SETTINGS.nParams; % number of parameter values
nS = SETTINGS.nStims;  % number of stimulus types
nN = SETTINGS.nNoises; % number of noise variations

% Add folders.
addpath(genpath(SETTINGS.codePath));
addpath(SETTINGS.filePath);

switch TASK
    case 0
        s = [ '\n' ...
            '  TASK [1]: Simulate data for given motif and logic\n' ...
            '     INDEX : [1-' num2str(nM*nL) ']\n' ...
            '     Simulations stored as Simulations_M#L#\n\n' ...
            '  TASK [2]: Consolidate simulations\n' ...
            '     Cell array stored as Simulations_FULL_S#\n\n' ...
            '  TASK [3]: Generate null models\n' ...
            '     Cell array stored as Simulations_NULL_S#\n\n' ...
            '  TASK [4]: Calculate true edge weights using selected ALGORITHM \n' ...
            '     INDEX : [1-' num2str(nM*nL*nS) ']\n' ...
            '     Edge weights stored as Results_[ALGORITHM]_M#L#S#\n\n' ...
            '  TASK [5]: Calculate null edge weights using selected ALGORITHM \n' ...
            '     INDEX : [1-' num2str(nN*nS*nP) ']\n' ...
            '     Edge weights stored as Results_[ALGORITHM]_N#S#A#\n\n' ...
            '  TASK [6]: Compile null results for selected ALGORITHM \n' ...
            '     INDEX : [1-' num2str(nN*nS) ']\n' ...
            '     Edge weights stored as Results_[ALGORITHM]_N#S#\n\n' ...
            '  TASK [7]: Analyzes edge weights for selected ALGORITHM \n' ...
            '     INDEX : [1-' num2str(nM*nL*nS) ']\n' ...
            '     Edge weights stored as Analysis_[ALGORITHM]_M#L#S#\n\n' ...
            '  TASK [8]: Creates CSVs for simulation data\n' ...
            '     CSVs are stored as M#L#S#N#.csv\n\n' ...
            '  TASK [9]: Creates CSVs for inference data\n' ...
            '     CSVs are stored as [ALGORITHM]_M#L#S#N#T#.csv\n\n' ...
            '  TASK [10]: Creates CSVs for summarized data\n' ...
            '     CSVs are stored as [ALGORITHM]_[EDGE].csv\n' ...
        '\n'];
        fprintf(s);
    case 1
        dir = start_parpool(SETTINGS, nP);
        index = (INDEX - 1)*nS + 1; % adjust index for stimulation types
        [output, code] = LM_simulate_data(index, 1);
        save([path 'Simulations_' code(1:4) '.mat'], '-struct', 'output');
        stop_parpool(dir);
    case 2
        for iStim = 1:nS
            for i = 1:(nM*nL)
                [code, iMotif, iLogic, ~] = get_code((i - 1)*nS + 1, 1);
                fprintf('Processing STIM[%d] MOTIF[%d] GATE[%d]\n', iStim, iMotif, iLogic);
                d = load([path 'Simulations_' code(1:4)]);
                
                for iNoise = 1:nN
                    noise = SETTINGS.noiseNames{iNoise};
                    output.(noise){iMotif, iLogic} = d.(noise)(:, :, iStim);
                end
            end
            
            save([path 'Simulations_FULL_S' num2str(iStim) '.mat'], '-struct', 'output');
        end
    case 3
        dir = start_parpool(SETTINGS, nP);
        for iStim = 1:nS
            output = LM_generate_nulls(iStim, '');
            save([path 'Simulations_NULL_S' num2str(iStim) '.mat'], '-struct', 'output');
        end
        stop_parpool(dir);
    case 4
        dir = start_parpool(SETTINGS, nP);
        [output, code] = LM_run_algorithm(INDEX, ALGORITHM, false, 1);
        save([path 'Results_' ALGORITHM '_' code '.mat'], '-struct', 'output');
        stop_parpool(dir);
    case 5
        dir = start_parpool(SETTINGS, nP);
        [output, code] = LM_run_nulls(INDEX, ALGORITHM, false, 1);
        save([path 'Results_' ALGORITHM '_' code '.mat'], '-struct', 'output');
        stop_parpool(dir);
    case 6
        index = (INDEX - 1)*nP + 1; % adjust index for parameter A values
        [output, code] = LM_compile_nulls(index, ALGORITHM, '');
        save([path 'Results_' ALGORITHM '_' code(1:4) '.mat'], '-struct', 'output');
    case 7
        dir = start_parpool(SETTINGS, nP);
        [output, code] = LM_calc_metrics(INDEX, ALGORITHM, false, '');
        save([path 'Analysis_' ALGORITHM '_' code '.mat'], '-struct', 'output');
        stop_parpool(dir);
    case 8
        LM_save_simulations('', '_I1');
    case 9
        LM_save_inference(ALGORITHM, '');
    case 10
        LM_save_summary(ALGORITHM);
end

end