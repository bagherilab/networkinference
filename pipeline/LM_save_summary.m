% -------------------------------------------------------------------------
% LM_SAVE_SUMMARY prints out csvs of the summarized inference results.
%
% Results include (algorithm) x (nE edges) files with (nM motifs) x (nL
% logics) x (nS stimulus types) x (nN noise types) x (nT time slices) rows
% with inferred weight (IW), null weight (NW), ES, and ERS mean, standard
% deviation, and chaos.
% -------------------------------------------------------------------------

function LM_save_summary(ALGORITHM)

SETTINGS = LM_SETTINGS();
nM = SETTINGS.nMotifs;
nL = SETTINGS.nLogics;
nS = SETTINGS.nStims;
nN = SETTINGS.nNoises;
nP = SETTINGS.nParams;
nT = SETTINGS.nTimeslices;
nV = SETTINGS.nValues;
nE = SETTINGS.nEdges;

for iEdge = 1:nE
    fromEdge = SETTINGS.edgeInds{iEdge}(1);
    toEdge = SETTINGS.edgeInds{iEdge}(2);

    filename = sprintf('%s_%s.csv', ALGORITHM, SETTINGS.edgeNames{iEdge});
    fid = fopen(filename, 'wt');
    fprintf(fid, ['motif,logic,stimulus,noise,slice,' ...
        'IW_mean,NW_mean,ES_mean,ERS_mean,' ...
        'IW_stdev,NW_stdev,ES_stdev,ERS_stdev,' ...
        'IW_chaos,NW_chaos,ES_chaos,ERS_chaos']);

    for iMotif = 1:nM
        for iLogic = 1:nL
            for iStim = 1:nS
                mls = ['M' num2str(iMotif) 'L' num2str(iLogic) 'S' num2str(iStim)];
                D = load(['Analysis_' ALGORITHM '_' mls]);

                for iNoise = 1:nN
                    fprintf('Creating CSV for MOTIF[%d] LOGIC[%d] NOISE[%d] STIM[%d]\n', ...
                        iMotif, iLogic, iNoise, iStim);
                    noise = SETTINGS.noiseNames{iNoise};

                    for iSlice = 1:nT
                        output = zeros(nV*3, 1);

                        for iValue = 1:nV
                            value = SETTINGS.values{iValue};
                            A = get_mat(D.(noise).(value).mean, ...
                                fromEdge, toEdge, nP, iSlice);
                            a = reshape(A, 1, []);
                            output(iValue) = mean(a);
                            output(iValue + 4) = std(a);
                            output(iValue + 8) = get_chaos(A);
                        end

                        fprintf(fid, '\n%d,%d,%d,%d,%d', iMotif, iLogic, iStim, iNoise, iSlice);
                        for i = 1:length(output)
                            fprintf(fid, ',%f', output(i));
                        end
                    end
                end
            end
        end
    end

    fclose(fid);
end
        
end