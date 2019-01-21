% -------------------------------------------------------------------------
% LM_MERGE_VALIDATION calculates similarity between the true and inferred
% networks for a given algorithm using various thresholding techniques.
% -------------------------------------------------------------------------

function output = LM_merge_validation(ALGORITHM)

SETTINGS = LM_SETTINGS();
nL = SETTINGS.nLogics; % number of logics
nM = SETTINGS.nMotifs; % number of motifs
nP = SETTINGS.nParams; % number of parameter values
N = 5; % number of nodes

[TO, FROM] = meshgrid(1:N, 1:N);

ids = [];
TP = [];
FP = [];
TN = [];
FN = [];

fid = fopen([ALGORITHM '_MERGED_ALLEDGES.csv'], 'wt');
fprintf(fid, 'CASE,SUBCASE,KA,KB,FROM,TO,IW,ES,ERS\n');

for iMotif = 1:nM
    A = SETTINGS.motifMatrices{iMotif};
    
    for iLogic = 1:nL
        for iStim = 2:3
            code = ['M' num2str(iMotif) 'L' num2str(iLogic) 'S' num2str(iStim)];
            D = load(['Analysis_' ALGORITHM '_' code '.mat']);
            D = D.noise_000;
            
            % Check for special cases.
            nanCodeCheck = sum(~cellfun('isempty', strfind(SETTINGS.nanCodes, code)));
            subsetCodeCheck = sum(~cellfun('isempty', strfind(SETTINGS.subsetCodes, code)));

            if nanCodeCheck == 1
                continue
            end

            % Adjust data for subset cases.
            if subsetCodeCheck == 1
                if iStim == 2
                    subset = [1 3 4 5];
                elseif iStim == 3
                    subset = [2 3 4 5];
                end
            else
                subset = 1:5;
            end

            for iA = 1:nP
                for iB = 1:nP
                    IW = D.IW.mean{iA,iB,1};
                    ES = D.ES.mean{iA,iB,1};
                    ERS = D.ERS.mean{iA,iB,1};
                    
                    iw = reshape(IW(subset,subset), 1, []);
                    es = reshape(ES(subset,subset), 1, []);
                    ers = reshape(ERS(subset,subset), 1, []);
                    from = reshape(FROM(subset,subset), 1, []);
                    to = reshape(TO(subset,subset), 1, []);
                    
                    edgelist = [from; to; iw; es; ers]';
                    edges = edgelist(edgelist(:,1) ~= edgelist(:,2), :);
                    
                    if iMotif == 6 && iLogic == 2
                        for i = 1:length(edges)
                            fprintf(fid, '%s,%d,%d,%d,%d,%d,%f,%f,%f\n', ...
                                code, iStim, iA, iB, edges(i,1), edges(i,2), edges(i,3), edges(i,4), edges(i,5));
                        end
                    end

                    m = get_measures(A, edges, subset);
                    ids = [ids; iMotif iLogic iStim iA iB];
                    
                    M = [m.N; m.IW; m.ES; m.IW_ERS; m.ES_ERS; m.ERS; m.IWes; m.IWers; m.IWiwers; m.IWesers]';
                    
                    TP = [TP; M(1,:)];
                    FP = [FP; M(2,:)];
                    TN = [TN; M(3,:)];
                    FN = [FN; M(4,:)];
                end
            end
        end
    end
end

output.ids = ids;
output.TP = TP;
output.FP = FP;
output.TN = TN;
output.FN = FN;

end