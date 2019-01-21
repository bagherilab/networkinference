% -------------------------------------------------------------------------
% LM_SAVE_STATES prints out pairwise statistical comparisons between
% (1) slices or (2) stimulus types, with multiple hypothesis testing
% correction.
%
% Results include (algorithm) x (nE edges) files with (nM motifs) x (nL
% logics) x (nS stimulus types or nT time slices) x (nN noise types) rows
% with inferred weight (IW), null weight (NW), ES, and ERS p-value and
% q-value.
% -------------------------------------------------------------------------

function LM_save_stats(CASE, ALGORITHM)

SETTINGS = LM_SETTINGS();
nM = SETTINGS.nMotifs;
nL = SETTINGS.nLogics;
nN = SETTINGS.nNoises;
nP = SETTINGS.nParams;
nV = SETTINGS.nValues;
nE = SETTINGS.nEdges;
N = 3;
pairs = {[1,2], [1,3], [2,3]};

switch CASE
    case 1
        name = 'SLICE';
        field = 'stimulus';  
    case 2
        name = 'STIM';
        field = 'slice';
end
      
for iEdge = 1:nE
    fromEdge = SETTINGS.edgeInds{iEdge}(1);
    toEdge = SETTINGS.edgeInds{iEdge}(2);
    pvals = zeros(nM, nL, N, nN, nV, N);
    qvals = zeros(nM, nL, N, nN, nV, N);
    
    for iMotif = 1:nM
        for iLogic = 1:nL
            D1 = load(['Analysis_' ALGORITHM '_M' num2str(iMotif) 'L' num2str(iLogic) 'S1']);
            D2 = load(['Analysis_' ALGORITHM '_M' num2str(iMotif) 'L' num2str(iLogic) 'S2']);
            D3 = load(['Analysis_' ALGORITHM '_M' num2str(iMotif) 'L' num2str(iLogic) 'S3']);
            D = {D1, D2, D3};

            for iS = 1:N
                for iNoise = 1:nN
                    fprintf('Creating CSV for MOTIF[%d] LOGIC[%d] NOISE[%d] [%d]\n', ...
                        iMotif, iLogic, iNoise, iS);
                    noise = SETTINGS.noiseNames{iNoise};

                    for iValue = 1:nV
                        value = SETTINGS.values{iValue};

                        for i = 1:N
                            switch CASE
                                case 1
                                    A1 = get_mat(D{iS}.(noise).(value).mean, ...
                                        fromEdge, toEdge, nP, pairs{i}(1));
                                    A2 = get_mat(D{iS}.(noise).(value).mean, ...
                                        fromEdge, toEdge, nP, pairs{i}(2));
                                case 2
                                    A1 = get_mat(D{pairs{i}(1)}.(noise).(value).mean, ...
                                        fromEdge, toEdge, nP, iS);
                                    A2 = get_mat(D{pairs{i}(2)}.(noise).(value).mean, ...
                                        fromEdge, toEdge, nP, iS);
                            end
                            
                            a1 = reshape(A1, 1, []);
                            a2 = reshape(A2, 1, []);

                            [~, p] = ttest2(a1, a2, 'Vartype', 'unequal');
                            pvals(iMotif, iLogic, iS, iNoise, iValue, i) = p;
                        end
                    end
                end
            end
        end
    end
        
    filename = sprintf('%s_%s_%s.csv', ALGORITHM, SETTINGS.edgeNames{iEdge}, name);
    fid = fopen(filename, 'wt');
    fprintf(fid, ['motif,logic,' field ',noise,', ...
        'p_IW_1_2,q_IW_1_2,p_IW_1_3,q_IW_1_3,p_IW_2_3,q_IW_2_3,', ...
        'p_NW_1_2,q_NW_1_2,p_NW_1_3,q_NW_1_3,p_NW_2_3,q_NW_2_3,', ...
        'p_ES_1_2,q_ES_1_2,p_ES_1_3,q_ES_1_3,p_ES_2_3,q_ES_2_3,', ...
        'p_ERS_1_2,q_ERS_1_2,p_ERS_1_3,q_ERS_1_3,p_ERS_2_3,q_ERS_2_3\n']);

    for iValue = 1:nV
        FDR = mafdr(reshape(pvals(:,:,:,:,iValue,:), 1, []), 'BHFDR', true);
        qvals(:,:,:,:,iValue,:) = reshape(FDR, nM, nL, N, nN, N);
    end

    for iMotif = 1:nM
        for iLogic = 1:nL
            for iS = 1:N
                for iNoise = 1:nN
                    fprintf(fid, '%d,%d,%d,%d', iMotif, iLogic, iS, iNoise);

                    for iValue = 1:nV
                        for i = 1:N
                            fprintf(fid, ',%f,%f', ...
                                pvals(iMotif, iLogic, iS, iNoise, iValue, i), ...
                                qvals(iMotif, iLogic, iS, iNoise, iValue, i));
                        end
                    end

                    fprintf(fid, '\n');
                end
            end
        end
    end

    fclose(fid);
end

end