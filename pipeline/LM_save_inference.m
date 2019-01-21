% -------------------------------------------------------------------------
% LM_SAVE_INFERENCE prints out csvs of the inference results.
%
% Results include (nM motifs) x (nL logics) x (nS stimulus types) x (nN
% noise types) x (nT time slices) files with (nP parameter A values) x
% (nP parameter B values) x (N nodes from) x (N nodes to) rows with
% inferred weight (IW), null weight (NW), ES, and ERS.
% -------------------------------------------------------------------------

function LM_save_inference(ALGORITHM, SUFFIX)

SETTINGS = LM_SETTINGS();
nM = SETTINGS.nMotifs;
nL = SETTINGS.nLogics;
nS = SETTINGS.nStims;
nN = SETTINGS.nNoises;
nP = SETTINGS.nParams;
nT = SETTINGS.nTimeslices;
nV = SETTINGS.nValues;
N = 5;

[X, Y] = meshgrid(1:N, 1:N);
from = reshape(X, 1, []);
to = reshape(Y, 1, []);

for iMotif = 1:nM
    LOCAL = LM_SETTINGS();

    for iLogic = 1:nL
        for iStim = 1:nS
            mls = ['M' num2str(iMotif) 'L' num2str(iLogic) 'S' num2str(iStim)];
            
            try
                D = load(['Analysis_' ALGORITHM '_' mls SUFFIX]);

                for iNoise = 1:nN
                    fprintf('Creating CSV for MOTIF[%d] LOGIC[%d] NOISE[%d] STIM[%d]\n', ...
                        iMotif, iLogic, iNoise, iStim);
                    noise = LOCAL.noiseNames{iNoise};
                    code = sprintf('M%d_L%d_S%d_N%d', iMotif, iLogic, iStim, iNoise);

                    for iSlice = 1:nT
                        filename = sprintf('%s_%s_T%d%s.csv', ALGORITHM, code, iSlice, SUFFIX);
                        fid = fopen(filename, 'wt');
                        fprintf(fid, 'paramA,paramB,nodeFrom,nodeTo,IW,NW,ES,ERS\n');
                        fclose(fid);

                        for iParamA = 1:nP
                            for iParamB = 1:nP
                                a = iParamA*ones(1, N*N);
                                b = iParamB*ones(1, N*N);
                                output = [a; b; from; to; zeros(4,N*N)];

                                for iValue = 1:nV
                                    value = LOCAL.values{iValue};
                                    c = D.(noise).(value).mean{iParamA, iParamB, iSlice};
                                    output(4 + iValue, :) = reshape(c', 1, []);
                                end

                                dlmwrite(filename, output', 'delimiter', ',', '-append');
                            end
                        end
                    end
                end
            catch
                fprintf('Unable to save %s\n', mls);
            end
        end
    end
end

end