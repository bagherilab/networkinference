% -------------------------------------------------------------------------
% LM_SAVE_HIDDEN prints out csvs of the simulation results with and without
% a hidden node.
%
% Results include (algorithm) x (nE edges) files with (nH hidden) x (nM
% motifs) x (nL logics) x (nS stimulus types) x (nN noise types) rows with
% inferred weight (IW), null weight (NW), ES, and ERS.
% -------------------------------------------------------------------------

function LM_save_hidden(ALGORITHM)

SETTINGS = LM_SETTINGS();
nM = SETTINGS.nMotifs;
nL = SETTINGS.nLogics;
nS = SETTINGS.nStims;
nN = SETTINGS.nNoises;
nE = SETTINGS.nEdges;
nV = SETTINGS.nValues;
nP = SETTINGS.nParams;
nH = 2;
hiddenNames = {'', '_HIDDEN'};

for iEdge = 1:nE
    fromEdge = SETTINGS.edgeInds{iEdge}(1);
    toEdge = SETTINGS.edgeInds{iEdge}(2);

    filename = sprintf('%s_%s_HIDDEN.csv', ALGORITHM, SETTINGS.edgeNames{iEdge});
    fid = fopen(filename, 'wt');
    fprintf(fid, 'hidden,motif,logic,stimulus,noise,IW,NW,ES,ERS');
    
    for iMotif = 1:nM
        for iLogic = 1:nL
            for iStim = 1:nS
                mls = ['M' num2str(iMotif) 'L' num2str(iLogic) 'S' num2str(iStim)];
                
                for iHidden = 1:nH
                    hidden = hiddenNames{iHidden};
                    
                    try
                        D = load(['Analysis_' ALGORITHM '_' mls hidden]);
                        
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
    
                            fprintf(fid, '\n%d,%d,%d,%d,%d', (iHidden - 1), iMotif, iLogic, iStim, iNoise);
                            for i = 1:length(output)
                                fprintf(fid, ',%f', output(i));
                            end
                        end
                    catch
                        fprintf('Unable to save %s for hidden %s\n', mls, hidden);
                    end
                end
            end
        end
    end

    fclose(fid);
end

% Save results for hidden node edge.
filename = sprintf('%s_MERGED_HIDDEN.csv', ALGORITHM);
fid = fopen(filename, 'wt');
fprintf(fid, 'node,direction,motif,logic,stimulus,noise,pA,pB,IW');
parser = get_parser(ALGORITHM);

for iDir = [1 2]
    for iMotif = [2 6]
        for iLogic = [1 2]
            for iStim = [2]
                mls = ['M' num2str(iMotif) 'L' num2str(iLogic) 'S' num2str(iStim)];

                try
                    D = load(['Analysis_' ALGORITHM '_' mls '_HIDDEN']);

                    for iNoise = 1:nN
                        noise = SETTINGS.noiseNames{iNoise};

                        for iTo = 1:5
                            for iParamA = 1:nP
                                for iParamB = 1:nP
                                    nanCodeCheck = sum(~cellfun('isempty', strfind(SETTINGS.nanCodes, mls)));
                                    if nanCodeCheck == 1
                                        iw = NaN;
                                    else
                                        if iDir == 1
                                            iw = parser(D.(noise).IW.mean{iParamA,iParamB}(6,iTo));
                                        else
                                            iw = parser(D.(noise).IW.mean{iParamA,iParamB}(iTo,6));
                                        end
                                    end

                                    fprintf(fid, '\n%d,%d,%d,%d,%d,%d,%d,%d,%d', iTo, iDir, ...
                                        iMotif, iLogic, iStim, iNoise, iParamA, iParamB, iw);
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

fclose(fid);
        
end