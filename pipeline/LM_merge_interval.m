% -------------------------------------------------------------------------
% LM_MERGE_INTERVAL consolidates ES and ERS across all true edges for the
% given algorithm across different intervals.
% -------------------------------------------------------------------------

function LM_merge_interval(ALGORITHM)

SETTINGS = LM_SETTINGS();
nL = SETTINGS.nLogics; % number of logics
nM = SETTINGS.nMotifs; % number of motifs
nT = SETTINGS.nTimeslices; % number of slices
params = [1,5,9,13,17]; % subset of parameters

fid = fopen([ALGORITHM '_MERGED_INTERVAL_DELTA.csv'], 'wt');
fprintf(fid, 'CODE,NOI,S,DIR,FROM,TO,INTERVAL,KA,KB,SIZE,TARGET,NP,NC,NL,IW,NW,ES,ERS\n');

fid2 = fopen([ALGORITHM '_MERGED_INTERVAL.csv'], 'wt');
fprintf(fid2, 'CODE,NOI,S,DIR,FROM,TO,INTERVAL,KA,KB,SIZE,TARGET,NP,NC,NL,IW,NW,ES,ERS\n');

% Time interval comparisons.
LHS = [1 1 2];
RHS = [2 3 3];
names = {'1-2', '1-3', '2-3'};

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
            [nEdges, ~, nStims] = size(E.(field));

            for iStim = 1:nStims
                for iEdge = 1:nEdges
                    e = E.(field)(iEdge,:,iStim);
                
                    if e(6) == 1
                        dir = e(1);
                        stim = e(2);
                        noi = e(3);
                        from = e(4);
                        to = e(5);
                        np = e(7);
                        nc = e(8);
                        nl = e(9);
                        
                        % Network size.
                        check = [code 'S' num2str(stim + 1)];
                        subsetCheck = sum(~cellfun('isempty', strfind(SETTINGS.subsetCodes, check)));
                        
                        if subsetCheck == 1
                            netsize = 4;
                        else
                            netsize = 5;
                        end

                        for iTime = 1:nT
                            lhs = LHS(iTime);
                            rhs = RHS(iTime);
                            
                            for iA = params
                                for iB = params
                                    fprintf(fid, '%s,%d,%d,%d,%d,%d,%s,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n', ...
                                        code(1:4), noi, stim, dir, ...
                                        from, to, names{iTime}, iA, iB, ...
                                        netsize, iTarget, np, nc, nl, ...
                                        D{stim}.IW.mean{iA,iB,lhs}(from,to) - D{stim}.IW.mean{iA,iB,rhs}(from,to), ...
                                        D{stim}.NW.mean{iA,iB,lhs}(from,to) - D{stim}.NW.mean{iA,iB,rhs}(from,to), ...
                                        D{stim}.ES.mean{iA,iB,lhs}(from,to) - D{stim}.ES.mean{iA,iB,rhs}(from,to), ...
                                        D{stim}.ERS.mean{iA,iB,lhs}(from,to) - D{stim}.ERS.mean{iA,iB,rhs}(from,to));
                                    
                                    fprintf(fid2, '%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n', ...
                                        code(1:4), noi, stim, dir, ...
                                        from, to, iTime, iA, iB, ...
                                        netsize, iTarget, np, nc, nl, ...
                                        D{stim}.IW.mean{iA,iB,iTime}(from,to), ...
                                        D{stim}.NW.mean{iA,iB,iTime}(from,to), ...
                                        D{stim}.ES.mean{iA,iB,iTime}(from,to), ...
                                        D{stim}.ERS.mean{iA,iB,iTime}(from,to));
                                end
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