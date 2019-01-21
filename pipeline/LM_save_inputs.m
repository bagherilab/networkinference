% -------------------------------------------------------------------------
% LM_SAVE_INPUTS prints out csvs of the simulation results with different
% inputs.
%
% Results include (algorithm) x (nE edges) files with (nI inputs) x (nM
% motifs) x (nL logics) x (nS stimulus types) x (nN noise types) x (nT time
% slices) rows with inferred weight (IW), null weight (NW), ES, and ERS.
% -------------------------------------------------------------------------

function LM_save_inputs(ALGORITHM)

SETTINGS = LM_SETTINGS();
nM = SETTINGS.nMotifs;
nL = SETTINGS.nLogics;
nS = SETTINGS.nStims;
nN = SETTINGS.nNoises;
nE = SETTINGS.nEdges;
nV = SETTINGS.nValues;
nP = SETTINGS.nParams;
nT = SETTINGS.nTimeslices;
nI = SETTINGS.nInputs;

for iEdge = 1:nE
    fromEdge = SETTINGS.edgeInds{iEdge}(1);
    toEdge = SETTINGS.edgeInds{iEdge}(2);

    filename = sprintf('%s_%s_INPUTS.csv', ALGORITHM, SETTINGS.edgeNames{iEdge});
    fid = fopen(filename, 'wt');
    fprintf(fid, 'input,motif,logic,stimulus,noise,slice,IW,NW,ES,ERS');
    
    for iMotif = 1:nM
        for iLogic = 1:nL
            for iStim = 1:nS
                mls = ['M' num2str(iMotif) 'L' num2str(iLogic) 'S' num2str(iStim)];
 
                try
                    for iInput = [2:nI 1]
                        input = SETTINGS.inputNames{iInput};
                        D = load(['Analysis_' ALGORITHM '_' mls input]);

                        for iNoise = 1:nN
                            noise = SETTINGS.noiseNames{iNoise};

                            for iSlice = 1:nT
                                output = zeros(nV, 1);

                                for iValue = 1:nV
                                    value = SETTINGS.values{iValue};
                                    A = get_mat(D.(noise).(value).mean, ...
                                        fromEdge, toEdge, nP, iSlice);
                                    a = reshape(A, 1, []);
                                    output(iValue) = mean(a);
                                end

                                fprintf(fid, '\n%d,%d,%d,%d,%d,%d', ...
                                    iInput, iMotif, iLogic, iStim, iNoise, iSlice);
                                for i = 1:length(output)
                                    fprintf(fid, ',%f', output(i));
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

    fclose(fid);
end
        
end