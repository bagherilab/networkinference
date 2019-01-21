% -------------------------------------------------------------------------
% LM_GNW contains the entire pipeline for loading experimental data,
% generating nulls, running the inferences, and then calculating metrics.
% The data was generated using GeneNetWeaver yeast networks with 50 nodes.
% -------------------------------------------------------------------------

function LM_GNW(TASK, ALGORITHM, INDEX)

SETTINGS = LM_GNW_SETTINGS();
nR = SETTINGS.nReplicates;
nN = SETTINGS.nNulls;
nT = SETTINGS.nTimepoints;
nI = SETTINGS.nIntervals;
N = SETTINGS.nNodes;
NAME = SETTINGS.name;
ID = SETTINGS.id;

switch TASK
% PARSE GNW NETWORKS ------------------------------------------------------
%
% Load edge lists from GeneNetWeaver output files to remove stimulated
% nodes that (i) have an outdegree of 0 such that stimulating it has no
% effect on the rest of the network and (ii) all direct downstream
% neighbors have an outdegree of 0 such that we cannot compare performance
% of node near and far from the stimulated node.
% -------------------------------------------------------------------------
    case 0
        REPS = SETTINGS.replicates;
        d = load([SETTINGS.path 'ijsvecs' ID '.mat']);
        D = load([SETTINGS.path 'simulation_data' ID '.mat']);

        fid = fopen([SETTINGS.filePath 'Cases_' NAME '.txt'], 'wt');

        for iRep = 1:nR
            % Make full digraph.
            field = ['ijs_' num2str(N) '_' num2str(iRep)];
            from = d.(field)(:,1);
            to = d.(field)(:,2);
            dir = d.(field)(:,3);
            G = digraph(from,to);

            % Find cases where node has at least one output edge (can infer).
            outdeg = outdegree(G);
            inferable = find(outdeg ~= 0);
            
            fprintf(fid, 'NETWORK %s | inferable | ', REPS{iRep});
            fprintf(fid, 'N = %d | ', length(inferable));
            fprintf(fid, join(string(inferable), ', '));

            s.inferable{iRep} = inferable;

            % Find cases where there is a positive incoming edge (can stimulate).
            from_p = from(dir > 0);
            to_p = to(dir > 0);
            H = digraph(from_p, to_p);
            indeg = indegree(H);
            stimulatable = find(indeg ~= 0);

            fprintf(fid, '\nNETWORK %s | stimulatable | ', REPS{iRep});
            fprintf(fid, 'N = %d | ', length(stimulatable));
            fprintf(fid, join(string(stimulatable), ', '));

            s.stimulatable{iRep} = stimulatable;
            s.edges{iRep}.from = from_p;
            s.edges{iRep}.to = to_p;

            % Get nodes we can both infer and stimulate.
            inf_and_stim = intersect(stimulatable, inferable)';

            fprintf(fid, '\nNETWORK %s | inferable + stimulatable | ', REPS{iRep});
            fprintf(fid, 'N = %d | ', length(inf_and_stim));
            fprintf(fid, join(string(inf_and_stim), ', '));

            % Get reachable nodes for each node we can stimulate and infer.
            data = D.(['data_' num2str(N) '_' num2str(iRep)]);
            s.subSets{iRep} = struct();
            subnodes = [];

            for ias = inf_and_stim
                reachable = dfsearch(G, ias);
                field = ['S' num2str(ias)];

                % Remove nodes that have flat trajectory.
                sims = data(:,reachable,ias);
                flats = reachable(std(sims) < 1e-5);
                subsets = setdiff(reachable, flats);
                s.flats{iRep}.(field) = flats;

                if length(subsets) > 3 % at least three nodes
                    s.subSets{iRep}.(field) = subsets';
                    subnodes = [subnodes ias];
                end
            end

            s.subNodes{iRep} = subnodes;
            s.adjacency{iRep} = full(adjacency(G));

            fprintf(fid, '\nNETWORK %s | inferable + stimulatable + 4 | ', REPS{iRep});
            fprintf(fid, 'N = %d | ', length(subnodes));
            
            if ~isempty(subnodes)
                fprintf(fid, join(string(subnodes), ', '));
            end

            fprintf(fid, '\n\n');
        end

        fclose(fid);
        
        save([SETTINGS.filePath 'Networks_' NAME '.mat'], '-struct', 's');
        
% LOAD GNW DATA -----------------------------------------------------------
%
% Load simulation data from GeneNetWeaver output files into consolidated
% .mat files named as Simulation_FULL_GNW.mat.
% -------------------------------------------------------------------------
    case 1
        d = load([SETTINGS.path 'simulation_data' ID '.mat']);
        s.sim = cell(nR, 1);

        for iRep = 1:nR
            field = ['data_' num2str(N) '_' num2str(iRep)];
            s.sim{iRep} = d.(field);
        end

        save([SETTINGS.filePath 'Simulations_FULL_' NAME '.mat'], '-struct', 's');

% GENERATE NULLS ----------------------------------------------------------
%
% Load consolidated simulation data from case 0 and generate nulls by
% shuffling along the perturbation axis for a given node and time point.
% Save as Simulation_NULL_GNW.mat
% -------------------------------------------------------------------------
    case 2
        d = load([SETTINGS.filePath 'Simulations_FULL_' NAME '.mat']);

        for iRep = 1:nR
            trues = d.sim{iRep};
            nulls = cell(nN,1);

            % Flatten true data by timepoint.
            ftrues = reshape(trues, nT, []);

            for i = 1:nN
                rng(i); % set seed
                inds = randi(N*N, nT, N); % random indices
                for iT = 1:nT
                    nulls{i}(iT,:) = ftrues(iT,inds(iT, :));
                end
            end

            s.null{iRep} = nulls;
        end
        save([SETTINGS.filePath 'Simulations_NULL_' NAME '.mat'], '-struct', 's');
        
% RUN INFERENCE ON TRUE DATA ----------------------------------------------
% 
% Run selected inference algorithm on the subsets of the selected cases of
% data. Save as Results_(ALGORITHM)_GNW_TRUE.mat
% -------------------------------------------------------------------------
    case 3
        dir = start_parpool(SETTINGS, nR);
        weights = cell(nR,1);
        
        parfor iRep = 1:nR
            LOCAL = LM_GNW_SETTINGS();
            S = load([LOCAL.filePath 'Networks_' NAME '.mat']);
            D = load([LOCAL.filePath 'Simulations_FULL_' NAME '.mat']); D = D.sim;
            algorithm = get_algorithm(ALGORITHM, LOCAL);
            subnodes = S.subNodes{iRep};
            
            for iStim = subnodes
                field = ['S' num2str(iStim)];
                subset = S.subSets{iRep}.(field);
                d = D{iRep}(:,subset,iStim);

                if strcmp(ALGORITHM, 'GENIE3')
                    LOCAL.regulatorIDs = 1:length(subset);
                    algorithm = get_algorithm(ALGORITHM, LOCAL);
                end

                for iInterval = 1:nI
                    fprintf('Running [%s] on REP[%d] STIM[%d] INT[%d]\n', ...
                        ALGORITHM, iRep, iStim, iInterval);

                    if strcmp(ALGORITHM, 'BANJO')
                        subname = sprintf('R%dS%dI%d', iRep, iStim, iInterval);
                    else
                        subname = '';
                    end

                    if strcmp(ALGORITHM, 'BANJO')
                        javaaddpath(LOCAL.banjoJar);
                    end

                    interval = LOCAL.intervals{iInterval};
                    dI = d(interval,:)';
                    w = -Inf*ones(N);
                    w(subset,subset) = algorithm(dI, subname);
                    weights{iRep}.(field){iInterval} = w;
                end
            end
        end

        save([SETTINGS.filePath 'Results_' ALGORITHM '_' NAME '_TRUE.mat'], 'weights');
        stop_parpool(dir);

% RUN INFERENCE ON NULL DATA ----------------------------------------------
% 
% Run selected inference algorithm on the nulls generated across stimulus
% conditions. Save as Results_(ALGORITHM)_GNW_NULL.mat
% -------------------------------------------------------------------------
    case 4
        dir = start_parpool(SETTINGS, nR);
        weights = cell(nR, nN, nI);
        
        parfor iNull = 1:nN
            LOCAL = LM_GNW_SETTINGS();
            algorithm = get_algorithm(ALGORITHM, LOCAL);

            for iRep = 1:nR
                D = load([LOCAL.filePath 'Simulations_NULL_' NAME '.mat']); D = D.null;
                d = D{iRep}{iNull};

                for iInterval = 1:nI
                    fprintf('Running [%s] on REP[%d] NULL[%d] INT[%d]\n', ...
                        ALGORITHM, iRep, iNull, iInterval);

                    if strcmp(ALGORITHM, 'BANJO')
                        subname = sprintf('R%dN%dI%d', iRep, iNull, iInterval);
                    else
                        subname = '';
                    end

                    if strcmp(ALGORITHM, 'BANJO')
                        javaaddpath(LOCAL.banjoJar);
                    end

                    interval = LOCAL.intervals{iInterval};
                    dI = d(interval,:)';
                    w = algorithm(dI, subname);
                    weights{iRep, iNull, iInterval} = w;
                end
            end
        end

        save([SETTINGS.filePath 'Results_' ALGORITHM '_' NAME '_NULL.mat'], 'weights');
        stop_parpool(dir);

% CALCULATE METRICS -------------------------------------------------------
%
% Calculate metrics for each network, relevant stimulated node, and time
% interval. Save as Analysis_(ALGORITHM)_GNW.mat
% -------------------------------------------------------------------------
    case 5
        S = load([SETTINGS.filePath 'Networks_' NAME '.mat']);
        path = [SETTINGS.filePath 'Results_' ALGORITHM];
        parser = get_parser(ALGORITHM);
        
        % Initialize containers.
        s.IW = cell(nR,1);
        s.NW = cell(nR,1);
        s.ES = cell(nR,1);
        s.ERS = cell(nR,1);

        % Load data.
        trueData = load([path '_' NAME '_TRUE.mat']); trueData = trueData.weights;
        nullData = load([path '_' NAME '_NULL.mat']); nullData = nullData.weights;
        
        for iRep = 1:nR
            subnodes = S.subNodes{iRep};
            
            for iStim = subnodes
                field = ['S' num2str(iStim)];
                subset = S.subSets{iRep}.(field);
                
                for iInterval = 1:nI
                    dTrue = trueData{iRep}.(field){iInterval};
                    dNull = nullData(iRep,:,iInterval);

                    % Containers for metrics.
                    es = -Inf*ones(N); es(subset, subset) = 0;
                    ers = -Inf*ones(N); ers(subset, subset) = 0;

                    % Get true weights and ranking.
                    [trueWeights, trueRanks] = get_vects(dTrue, subset, parser);
                    
                    % Get null weights and rankings.
                    for i = 1:nN
                        [nullWeights, nullRanks] = get_vects(dNull{i}, subset, parser);
                        es(subset, subset) = es(subset, subset) + ...
                            (trueWeights > nullWeights) + 0.5*(trueWeights == nullWeights);
                        ers(subset, subset) = ers(subset, subset) + ...
                            (trueRanks < nullRanks) + 0.5*(trueRanks == nullRanks);
                    end
                    
                    s.ES{iRep}.(field){iInterval} = es/nN;
                    s.ERS{iRep}.(field){iInterval} = ers/nN;
                    
                    % Add in true weights.
                    trueMat = dTrue;
                    trueMat(subset, subset) = parser(trueMat(subset, subset));
                    s.IW{iRep}.(field){iInterval} = trueMat;

                    % Add in null weights.
                    nullMat = zeros(N,N,nN);
                    for i = 1:nN
                        nullMat(:,:,i) = parser(dNull{i});
                    end
                    s.NW{iRep}.(field){iInterval} = mean(nullMat, 3);
                end
            end
        end
        
        save([SETTINGS.filePath 'Analysis_' ALGORITHM '_' NAME '.mat'], '-struct', 's');

% CONSOLIDATE INTERVAL DATA -----------------------------------------------
% 
% Generate CSV summarizing true edge metrics under different intervals.
% Column headers: network, S (stimulated node), A (from node), B (to node),
% ES, ERS. Save as (ALGORITHM)_GNW_INTERVAL.csv
% -------------------------------------------------------------------------
    case 6
        REPS = SETTINGS.replicates;
        S = load([SETTINGS.filePath 'Networks_' NAME '.mat']);
        D = load([SETTINGS.filePath 'Analysis_' ALGORITHM '_' NAME '.mat']);
        
        fid = fopen([ALGORITHM '_' NAME '_INTERVAL_DELTA.csv'], 'wt');
        fprintf(fid, 'CODE,NOI,S,DIR,FROM,TO,INTERVAL,SIZE,TARGET,NP,NC,NL,IW,NW,ES,ERS\n');
        
        fid2 = fopen([ALGORITHM '_' NAME '_INTERVAL.csv'], 'wt');
        fprintf(fid2, 'CODE,NOI,S,DIR,FROM,TO,INTERVAL,SIZE,TARGET,NP,NC,NL,IW,NW,ES,ERS\n');
        
        % Time interval comparisons.
        LHS = [1 1 2];
        RHS = [2 3 3];
        names = {'1-2', '1-3', '2-3'};

        for iRep = 1:nR
            A = S.adjacency{iRep};
            subnodes = S.subNodes{iRep};

            if isempty(subnodes)
                continue
            end

            [E, targets] = get_fan_edges(A, subnodes);
                
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
                            sfield = ['S' num2str(stim)];
                            netsize = length(S.subSets{iRep}.(sfield));
                            
                            % Remove cases for flat trajectories.
                            flat = S.flats{iRep}.(sfield);
                            
                            if ~isempty(find(flat == to, 1))
                                continue
                            end
                            
                            if ~isempty(find(flat == from, 1))
                                continue
                            end

                            for iTime = 1:nI
                                lhs = LHS(iTime);
                                rhs = RHS(iTime);

                                fprintf(fid, '%s,%d,%d,%d,%d,%d,%s,%d,%d,%d,%d,%d,%f,%f,%f,%f\n', ...
                                    REPS{iRep}, noi, stim, dir, ...
                                    from, to, names{iTime}, ...
                                    netsize, iTarget, np, nc, nl, ...
                                    D.IW{iRep}.(sfield){lhs}(from,to) - D.IW{iRep}.(sfield){rhs}(from,to), ...
                                    D.NW{iRep}.(sfield){lhs}(from,to) - D.NW{iRep}.(sfield){rhs}(from,to), ...
                                    D.ES{iRep}.(sfield){lhs}(from,to) - D.ES{iRep}.(sfield){rhs}(from,to), ...
                                    D.ERS{iRep}.(sfield){lhs}(from,to) - D.ERS{iRep}.(sfield){rhs}(from,to));
                               
                                fprintf(fid2, '%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n', ...
                                    REPS{iRep}, noi, stim, dir, ...
                                    from, to, iTime, ...
                                    netsize, iTarget, np, nc, nl, ...
                                    D.IW{iRep}.(sfield){iTime}(from,to), ...
                                    D.NW{iRep}.(sfield){iTime}(from,to), ...
                                    D.ES{iRep}.(sfield){iTime}(from,to), ...
                                    D.ERS{iRep}.(sfield){iTime}(from,to));
                            end
                        end
                    end
                end
            end               
        end
        
        fclose(fid);
        fclose(fid2);
        
% CONSOLIDATE STIMULUS DATA -----------------------------------------------
%
% Generate CSV summarizing true edge metrics under different stimulus.
% Column headers: network, NOI (node of interest), S (other stimulated
% node), DIR (direction, 1 = from, -1 = to), FROM, TO, IW, NW, ES, ERS.
% Save as (ALGORITHM)_GNW_STIMULUS.csv
% -------------------------------------------------------------------------
    case 7
        REPS = SETTINGS.replicates;
        S = load([SETTINGS.filePath 'Networks_' NAME '.mat']);
        D = load([SETTINGS.filePath 'Analysis_' ALGORITHM '_' NAME '.mat']);
        
        fid = fopen([ALGORITHM '_' NAME '_STIMULUS_DELTA.csv'], 'wt');
        fprintf(fid, 'CODE,NOI,S,DIR,FROM,TO,SIZE,TARGET,NP,NC,NL,IW,NW,ES,ERS\n');
        
        fid2 = fopen([ALGORITHM '_' NAME '_STIMULUS.csv'], 'wt');
        fprintf(fid2, 'CODE,NOI,S,DIR,FROM,TO,SIZE,TARGET,NP,NC,NL,IW,NW,ES,ERS\n');

        for iRep = 1:nR
            A = S.adjacency{iRep};
            subnodes = S.subNodes{iRep};

            if isempty(subnodes)
                continue
            end

            [E, targets] = get_fan_edges(A, subnodes);
                
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
                        
                        if isempty(iNOI)
                            continue
                        end
                        
                        noi_s = subgroup(1, 2, iNOI);
                        noi_f = ['S' num2str(noi_s)];
                        
                        % Remove cases for flat trajectories.
                        noi_flat = S.flats{iRep}.(noi_f);
                        
                        if ~isempty(find(noi_flat == to, 1))
                            continue
                        end

                        if ~isempty(find(noi_flat == from, 1))
                            continue
                        end
                        
                        noi_nodes = S.subSets{iRep}.(noi_f);
                        
                        fprintf(fid2, '%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n', ...
                            REPS{iRep}, noi, noi_s, dir, from, to, ...
                            length(noi_nodes), iTarget, np, nc, nl, ...
                            D.IW{iRep}.(noi_f){1}(from,to), ...
                            D.NW{iRep}.(noi_f){1}(from,to), ...
                            D.ES{iRep}.(noi_f){1}(from,to), ...
                            D.ERS{iRep}.(noi_f){1}(from,to));

                        for iunoi = iUNOI
                            unoi_s = subgroup(1, 2, iunoi);
                            unoi_f = ['S' num2str(unoi_s)];
                            
                            % Remove cases for flat trajectories.
                            unoi_flat = S.flats{iRep}.(unoi_f);

                            if ~isempty(find(unoi_flat == to, 1))
                                continue
                            end
                            
                            if ~isempty(find(unoi_flat == from, 1))
                                continue
                            end
                            
                            % Get size of overlapping nodes.
                            unoi_nodes = S.subSets{iRep}.(unoi_f);
                            netsize = length(intersect(noi_nodes, unoi_nodes));

                            fprintf(fid, '%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n', ...
                                REPS{iRep}, noi, unoi_s, dir, from, to, ...
                                netsize, iTarget, np, nc, nl, ...
                                D.IW{iRep}.(noi_f){1}(from,to) - D.IW{iRep}.(unoi_f){1}(from,to), ...
                                D.NW{iRep}.(noi_f){1}(from,to) - D.NW{iRep}.(unoi_f){1}(from,to), ...
                                D.ES{iRep}.(noi_f){1}(from,to) - D.ES{iRep}.(unoi_f){1}(from,to), ...
                                D.ERS{iRep}.(noi_f){1}(from,to) - D.ERS{iRep}.(unoi_f){1}(from,to));

                            fprintf(fid2, '%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n', ...
                                REPS{iRep}, noi, unoi_s, dir, from, to, ...
                                length(unoi_nodes), iTarget, np, nc, nl, ...
                                D.IW{iRep}.(unoi_f){1}(from,to), ...
                                D.NW{iRep}.(unoi_f){1}(from,to), ...
                                D.ES{iRep}.(unoi_f){1}(from,to), ...
                                D.ERS{iRep}.(unoi_f){1}(from,to));
                        end
                    end
                end
            end               
        end
        
        fclose(fid);
        fclose(fid2);

% PLOT NETWORKS -----------------------------------------------------------
    case 8
        close all;
        S = load([SETTINGS.filePath 'Networks_' NAME '.mat']);
        A = S.adjacency{INDEX};
        G = digraph(A);
        subnodes = S.subNodes{INDEX};

        if isempty(subnodes)
            targets = [];
            for i = 1:N
                parents = A(:,i);
                if sum(parents) > 1
                    targets = [targets i];
                end
            end
        else
            [~, targets] = get_fan_edges(A, subnodes);
        end

        stim = S.stimulatable{INDEX};
        inf = S.inferable{INDEX};
        edges = S.edges{INDEX};

        subplot(1,2,1);
        h1 = plot(G, 'NodeColor', [0.6 0.6 0.6], 'EdgeColor', [0.1 0.1 0.7]);
        title('Stimulatable (g) / Inferable (r) / Both (b)');
        highlight(h1, edges.from, edges.to, 'EdgeColor', [0.7 0.1 0.1]);
        highlight(h1, stim, 'NodeColor', [0.1 0.7 0.1]);
        highlight(h1, inf, 'NodeColor', [0.7 0.1 0.1]);
        highlight(h1, intersect(inf, stim), 'NodeColor', [0.1 0.1 0.7]);

        subplot(1,2,2);
        h2 = plot(G, 'NodeColor', [0.6 0.6 0.6], 'EdgeColor', [0.1 0.1 0.7]);
        title('Stim + Inf + 4 (b) / Fan Targets (o) / Both (p) ');
        highlight(h2, edges.from, edges.to, 'EdgeColor', [0.7 0.1 0.1]);
        highlight(h2, subnodes, 'NodeColor', [0.1 0.1 0.9]);
        highlight(h2, targets, 'Marker', 'd', 'NodeColor', [1 0.5 0]);
        highlight(h2, intersect(subnodes, targets), 'NodeColor', [0.5 0 0.7]);

        h = gcf;
        h.Position = [100 100 1600 800];

% RUN VALIDATION ----------------------------------------------------------
%
% Generate validation .mat file of similarities between the inferred and
% the true networks.
% -------------------------------------------------------------------------
    case 9
        REPS = SETTINGS.replicates;
        S = load([SETTINGS.filePath 'Networks_' NAME '.mat']);
        D = load([SETTINGS.filePath 'Analysis_' ALGORITHM '_' NAME '.mat']);
        [TO, FROM] = meshgrid(1:N, 1:N);
        
        ids = [];
        TP = [];
        FP = [];
        TN = [];
        FN = [];

        fid = fopen([ALGORITHM '_' NAME '_ALLEDGES.csv'], 'wt');
        fprintf(fid, 'CASE,SUBCASE,FROM,TO,IW,ES,ERS\n');

        for iRep = 1:nR
            A = S.adjacency{iRep};
            subnodes = S.subNodes{iRep};

            if isempty(subnodes)
                continue
            end
            
            for iStim = subnodes
                field = ['S' num2str(iStim)];
                subset = S.subSets{iRep}.(field);
                
                IW = D.IW{iRep}.(field){1};
                ES = D.ES{iRep}.(field){1};
                ERS = D.ERS{iRep}.(field){1};

                iw = reshape(IW(subset,subset), 1, []);
                es = reshape(ES(subset,subset), 1, []);
                ers = reshape(ERS(subset,subset), 1, []);
                from = reshape(FROM(subset,subset), 1, []);
                to = reshape(TO(subset,subset), 1, []);
                
                edgelist = [from; to; iw; es; ers]';
                edges = edgelist(edgelist(:,1) ~= edgelist(:,2), :);

                if iRep == 1
                    for i = 1:length(edges)
                        fprintf(fid, '%s,%d,%d,%d,%f,%f,%f\n', ...
                            REPS{iRep}, iStim, edges(i,1), edges(i,2), edges(i,3), edges(i,4), edges(i,5));
                    end
                end

                m = get_measures(A, edges, subset);
                ids = [ids; iRep iStim];

                M = [m.N; m.IW; m.ES; m.IW_ERS; m.ES_ERS; m.ERS; m.IWes; m.IWers; m.IWiwers; m.IWesers]';

                TP = [TP; M(1,:)];
                FP = [FP; M(2,:)];
                TN = [TN; M(3,:)];
                FN = [FN; M(4,:)];
            end
        end
        
        output.ids = ids;
        output.TP = TP;
        output.FP = FP;
        output.TN = TN;
        output.FN = FN;
        
        save([SETTINGS.filePath 'Validation_' ALGORITHM '_' NAME '.mat'], '-struct', 'output');

% SAVE VALIDATION ---------------------------------------------------------
    case 10
        if INDEX == 0
            D = load([SETTINGS.filePath 'Validation_' ALGORITHM '_' NAME '.mat']);
        else
            D = load(['Validation_' ALGORITHM '.mat']); NAME = 'MERGED';
        end

        TP = D.TP;
        FP = D.FP;
        TN = D.TN;
        FN = D.FN;

        ACCURACY = (TP + TN)./(TP + TN + FP + FN);
        SENSITIVITY = TP./(TP + FN);
        SPECIFICITY = TN./(TN + FP);
        FALLOUT = FP./(FP + TN);
        PRECISION = TP./(TP + FP);

        % Verify that calculations don't have any division by zero errors.
        if sum(sum(isnan(SENSITIVITY))) > 0
            SENSITIVITY
        end

        if sum(sum(isnan(ACCURACY))) > 0
            ACCURACY
        end

        if sum(sum(isnan(SPECIFICITY))) > 0
            SPECIFICITY
        end

        if sum(sum(isnan(FALLOUT))) > 0
            FALLOUT
        end

        if sum(sum(isnan(PRECISION))) > 0
            nans = find(isnan(mean(PRECISION,2)));
            PRECISION(nans,:) = [];
        end

        D.SENSITIVITY = SENSITIVITY;
        D.ACCURACY = ACCURACY;
        D.SPECIFICITY = SPECIFICITY;
        D.PRECISION = PRECISION;

        fid = fopen([ALGORITHM '_' NAME '_VALIDATION.csv'], 'wt');
        fprintf(fid, 'CASE,LHS,RHS,N,MEAN,STD,MIN,MAX,PVAL,PVAL_LEFT,PVAL_RIGHT\n');

        scores = { 'TP', 'FP', 'TN', 'FN', 'ACCURACY', 'SENSITIVITY', 'SPECIFICITY', 'FALLOUT', 'PRECISION' };
        cases = { 'N', 'IW', 'ES', 'IWERS', 'ESERS', 'ERS', 'IWes', 'IWers', 'IWiwers', 'IWesers' };

        % Print out actual scores.
        for iScore = 1:length(scores)
            s = scores{iScore};

            for iCase = 1:length(cases)
                c = cases{iCase};

                a = D.(s)(:,iCase);
                a(isnan(a)) = [];

                fprintf(fid, '%s,%s,*,%d,%f,%f,%f,%f,*,*,*\n', ...
                    s, c, size(a,1), mean(a), std(a), min(a), max(a));
            end
        end

        % Print out deltas.
        for iScore = 1:length(scores)
            s = scores{iScore};

            for i = 1:length(cases)
                lhs = cases{i};

                for j = 1:length(cases)
                    rhs = cases{j};

                    a = D.(s)(:,i);
                    b = D.(s)(:,j);

                    pval = signrank(a, b);
                    pval_left = signrank(a, b, 'tail', 'left');
                    pval_right = signrank(a, b, 'tail', 'right');

                    fprintf(fid, '%s,%s,%s,%d,%f,%f,%f,%f,%f,%f,%f\n', ...
                        s, lhs, rhs, size(D.ids,1), ...
                        mean(a - b), std(a - b), min(a - b), max(a - b), ...
                        pval, pval_left, pval_right);
                end
            end
        end

        fclose(fid);
end

end

function [weights, ranks] = get_vects(A, subset, parser)
    n = length(subset);
    weight_vect = reshape(parser(A(subset, subset)), [], 1);
    [~, inds] = sort(weight_vect, 'descend');
    rank_vect(inds) = 1:(n*n);
    ranks = reshape(rank_vect, n, []);
    weights = reshape(weight_vect, n, []);
end