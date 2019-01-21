% -------------------------------------------------------------------------
% LM_MERGE_STIMULUS consolidates ES and ERS across all true edges for the
% given algorithm across different stimulus conditions.
% -------------------------------------------------------------------------

function LM_merge_stimulus(ALGORITHM)

SETTINGS = LM_SETTINGS();
nL = SETTINGS.nLogics; % number of logics
nM = SETTINGS.nMotifs; % number of motifs
params = [1,5,9,13,17]; % subset of parameters

fid = fopen([ALGORITHM '_MERGED_STIMULUS_DELTA.csv'], 'wt');
fprintf(fid, 'CODE,NOI,S,DIR,FROM,TO,KA,KB,SIZE,TARGET,NP,NC,NL,IW,NW,ES,ERS\n');

fid2 = fopen([ALGORITHM '_MERGED_STIMULUS.csv'], 'wt');
fprintf(fid2, 'CODE,NOI,S,DIR,FROM,TO,KA,KB,SIZE,TARGET,NP,NC,NL,IW,NW,ES,ERS\n');

for iMotif = 1:nM
    A = SETTINGS.motifMatrices{iMotif};
    
    for iLogic = 1:nL
        % Go to next case if nan.
        if (iLogic == 1 && iMotif ~= 3)
            continue
        end

        [E, targets] = get_fan_edges(A, [1 2]);

        % Load data.
        code = ['M' num2str(iMotif) 'L' num2str(iLogic)];
        D1 = load(['Analysis_' ALGORITHM '_' code 'S2']);
        D{1} = D1.noise_000;
        D2 = load(['Analysis_' ALGORITHM '_' code 'S3']);
        D{2} = D2.noise_000;

        for iTarget = targets
            field = ['T' num2str(iTarget)];
            e = E.(field);
            [nEdges, ~, nStims] = size(e);

            if nStims < 2
                continue
            end

            checks = sum(e(:,6,:), 3);
            
            for iEdge = 1:nEdges
                if checks(iEdge) > 1
                    group = e(iEdge,:,:);
                    subgroup = group(:,:,group(:,6,:) == 1);

                    dir = subgroup(1,1,1);
                    noi = subgroup(1,3,1);
                    from = subgroup(1,4,1);
                    to = subgroup(1,5,1);
                    
                    np = subgroup(1,7,1);
                    nc = subgroup(1,8,1);
                    nl = subgroup(1,9,1);

                    iNOI = find(subgroup(:,2,:) == subgroup(:,3,:));
                    iUNOI = setdiff(1:size(subgroup,3), iNOI);
                    
                    noi_s = e(iEdge, 2, iNOI);
                    noi_nodes = get_subset(SETTINGS, code, noi_s + 1);
                    
                    for iA = params
                        for iB = params
                            fprintf(fid2, '%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n', ...
                                code(1:4), noi, noi_s, dir, from, to, iA, iB, ...
                                length(noi_nodes), iTarget, np, nc, nl, ...
                                D{noi_s}.IW.mean{iA,iB,1}(from, to), ...
                                D{noi_s}.NW.mean{iA,iB,1}(from, to), ...
                                D{noi_s}.ES.mean{iA,iB,1}(from, to), ...
                                D{noi_s}.ERS.mean{iA,iB,1}(from, to));
                        end
                    end
                                
                    for iunoi = iUNOI
                        unoi_s = e(iEdge, 2, iunoi);
                        
                        % Get size of overlapping nodes.
                        unoi_nodes = get_subset(SETTINGS, code, unoi_s + 1);
                        netsize = length(intersect(noi_nodes, unoi_nodes));

                        for iA = params
                            for iB = params
                                fprintf(fid, '%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n', ...
                                    code(1:4), noi, unoi_s, dir, from, to, iA, iB, ...
                                    netsize, iTarget, np, nc, nl, ...
                                    D{noi_s}.IW.mean{iA,iB,1}(from, to) - D{unoi_s}.IW.mean{iA,iB,1}(from, to), ...
                                    D{noi_s}.NW.mean{iA,iB,1}(from, to) - D{unoi_s}.NW.mean{iA,iB,1}(from, to), ...
                                    D{noi_s}.ES.mean{iA,iB,1}(from, to) - D{unoi_s}.ES.mean{iA,iB,1}(from, to), ...
                                    D{noi_s}.ERS.mean{iA,iB,1}(from, to) - D{unoi_s}.ERS.mean{iA,iB,1}(from, to));
                                
                                fprintf(fid2, '%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n', ...
                                    code(1:4), noi, unoi_s, dir, from, to, iA, iB, ...
                                    length(unoi_nodes), iTarget, np, nc, nl, ...
                                    D{unoi_s}.IW.mean{iA,iB,1}(from, to), ...
                                    D{unoi_s}.NW.mean{iA,iB,1}(from, to), ...
                                    D{unoi_s}.ES.mean{iA,iB,1}(from, to), ...
                                    D{unoi_s}.ERS.mean{iA,iB,1}(from, to));
                            end
                        end
                    end
                end
            end
        end
    end
end

fclose(fid);
fclose(fid2);

end

function subset = get_subset(SETTINGS, code, stim)

subsetCheck = sum(~cellfun('isempty', strfind(SETTINGS.subsetCodes, [code 'S' num2str(stim)])));

if subsetCheck == 1
    if stim == 2
        subset = [1 3 4 5];
    elseif stim == 3
        subset = [2 3 4 5];
    end
else
    subset = 1:5;
end

end