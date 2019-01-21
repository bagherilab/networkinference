% -------------------------------------------------------------------------
% LM_MERGE_INFERENCE prints out csvs of the merged inference results.
% -------------------------------------------------------------------------

function LM_merge_inference(ALGORITHM)

SETTINGS = LM_SETTINGS();
nV = SETTINGS.nValues;
N = 5;

[X, Y] = meshgrid(1:N, 1:N);
from = reshape(X, 1, []);
to = reshape(Y, 1, []);

filename = [ALGORITHM '_MERGED_INFERENCE.csv'];
fid = fopen(filename, 'wt');
fprintf(fid, 'case,nodeFrom,nodeTo,IW,NW,ES,ERS\n');
fclose(fid);

% CHAMELEON

mls = 'M6L1S2';
D = load(['Analysis_' ALGORITHM '_' mls]);
D = D.noise_000;

paramA = [  9,  9,  5,  5, 13, 13, 17, 17, 17];
paramB = [  9, 17,  9, 17, 13, 17,  9,  5, 13];
M = zeros(N, N, 9);
output = [ones(1,N*N); from; to; zeros(4,N*N)];

for iValue = 1:nV
    value = SETTINGS.values{iValue};
    
    for i = 1:9
        M(:,:,i) = D.(value).mean{paramA(i), paramB(i), 1};
    end
    
    m = std(M, 1, 3);
    output(3 + iValue, :) = reshape(m', 1, []);
end

dlmwrite(filename, output', 'delimiter', ',', '-append');

% REVERSE CHAMELEON

paramA = [5, 9, 13];
motifs = [2, 3, 6];
logics = [2, 3, 1];
M = zeros(N, N, 9);
output = [2*ones(1,N*N); from; to; zeros(4,N*N)];

for iValue = 1:nV
    value = SETTINGS.values{iValue};
    
    for i = 1:3
        for j = 1:3
            mls = ['M' num2str(motifs(i)) 'L' num2str(logics(j)) 'S2'];
            D = load(['Analysis_' ALGORITHM '_' mls]);
            D = D.noise_000;
            M(:,:,(i - 1)*3 + j) = D.(value).mean{paramA(j), 9, 1};
        end
    end
   
    m = std(M, 1, 3);
    output(3 + iValue, :) = reshape(m', 1, []);
end


dlmwrite(filename, output', 'delimiter', ',', '-append');

end