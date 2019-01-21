% -------------------------------------------------------------------------
% LM_SAVE_BOTH prints out csvs of the simulation results with different
% inputs and with and without a hidden node.
%
% Results include (algorithm) x (nE edges) files with (nH hidden) x (nI
% inputs) x (nM motifs) x (nL logics) x (nS stimulus types) x (nN noise
% types) rows with inferred weight (IW), null weight (NW), ES, and ERS.
% -------------------------------------------------------------------------

function LM_save_both(ALGORITHM)

SETTINGS = LM_SETTINGS();
nM = SETTINGS.nMotifs;
nL = SETTINGS.nLogics;
nS = SETTINGS.nStims;
nN = SETTINGS.nNoises;
nE = SETTINGS.nEdges;
nV = SETTINGS.nValues;
nP = SETTINGS.nParams;
nI = SETTINGS.nInputs;
nH = 2;
hiddenNames = {'', '_HIDDEN'};

for iEdge = 1:nE
    fromEdge = SETTINGS.edgeInds{iEdge}(1);
    toEdge = SETTINGS.edgeInds{iEdge}(2);

    filename = sprintf('%s_%s_BOTH.csv', ALGORITHM, SETTINGS.edgeNames{iEdge});
    fid = fopen(filename, 'wt');
    fprintf(fid, 'input,hidden,motif,logic,stimulus,noise,IW,NW,ES,ERS');
    
    for iMotif = 1:nM
        for iLogic = 1:nL
            for iStim = 1:nS
                mls = ['M' num2str(iMotif) 'L' num2str(iLogic) 'S' num2str(iStim)];
                
                try
                    for iInput = [2:nI 1]
                        input = SETTINGS.inputNames{iInput};
                        
                        for iHidden = 1:nH 
                            hidden = hiddenNames{iHidden};
                            D = load(['Analysis_' ALGORITHM '_' mls input hidden]);

                            for iNoise = 1:nN
                                noise = SETTINGS.noiseNames{iNoise};
                                output = zeros(nV, 1);

                                for iValue = 1:nV
                                    value = SETTINGS.values{iValue};
                                    A = get_mat(D.(noise).(value).mean, ...
                                        fromEdge, toEdge, nP, 1);
                                    a = reshape(A, 1, []);
                                    output(iValue) = mean(a);
                                end

                                fprintf(fid, '\n%d,%d,%d,%d,%d,%d', iInput, (iHidden - 1), iMotif, iLogic, iStim, iNoise);
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