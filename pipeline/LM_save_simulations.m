% -------------------------------------------------------------------------
% LM_SAVE_SIMULATIONS prints out csvs of the simulation results.
%
% Results include (nM motifs) x (nL logics) x (nS stimulus types) x (nN
% noise types) files with (nP parameter A values) x (nP parameter B
% values) x (N nodes) x (T timepoints) rows.
% -------------------------------------------------------------------------

function LM_save_simulations(SUFFIX, CODE)

SETTINGS = LM_SETTINGS();
nM = SETTINGS.nMotifs;
nL = SETTINGS.nLogics;
nS = SETTINGS.nStims;
nN = SETTINGS.nNoises;
nP = SETTINGS.nParams;
N = 5;
T = length(SETTINGS.timePoints);

[X, Y] = meshgrid(1:N, 1:T);
node = reshape(X, 1, []);
time = reshape(Y, 1, []);

for iStim = 1:nS
    D = load(['Simulations_FULL_S' num2str(iStim) SUFFIX]);

    for iMotif = 1:nM
        for iLogic = 1:nL
            for iNoise = 1:nN
                fprintf('Creating CSV for MOTIF[%d] LOGIC[%d] NOISE[%d] STIM[%d]\n', ...
                    iMotif, iLogic, iNoise, iStim);
                noise = SETTINGS.noiseNames{iNoise};
                code = sprintf('M%d_L%d_S%d_N%d', iMotif, iLogic, iStim, iNoise);

                filename = sprintf(['%s' CODE '.csv'], code);
                fid = fopen(filename, 'wt');
                fprintf(fid, 'paramA,paramB,node,time,value\n');
                fclose(fid);

                for iParamA = 1:nP
                    for iParamB = 1:nP
                        a = iParamA*ones(1, N*T);
                        b = iParamB*ones(1, N*T);
                        c = D.(noise){iMotif, iLogic}{iParamA, iParamB};
                        output = [a; b; node; time; reshape(c', 1, [])];
                        dlmwrite(filename, output', 'delimiter', ',', '-append');
                    end
                end
            end
        end
    end
end

end