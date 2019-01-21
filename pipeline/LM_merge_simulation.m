% -------------------------------------------------------------------------
% LM_MERGE_SIMULATION prints out csvs of the merged simulations.
% -------------------------------------------------------------------------

function LM_merge_simulation()

SETTINGS = LM_SETTINGS();
N = 5;
T = length(SETTINGS.timePoints);

[X, Y] = meshgrid(1:N, 1:T);
node = reshape(X, 1, []);
time = reshape(Y, 1, []);

D = load('Simulations_FULL_S2');

filename = '_MERGED_SIMULATION.csv';
fid = fopen(filename, 'wt');
fprintf(fid, 'case,node,time,mean,std\n');
fclose(fid);

% CHAMELEON

iMotif = 6;
iLogic = 1;
paramA = [  9,  9,  5,  5, 13, 13, 17, 17, 17];
paramB = [  9, 17,  9, 17, 13, 17,  9,  5, 13];
M = zeros(N, T, 9);
output = [ones(1,N*T); node; time; zeros(2,N*T)];
 
for i = 1:9
    M(:,:,i) = D.noise_000{iMotif, iLogic}{paramA(i), paramB(i)};
end
   
m1 = mean(M, 3);
m2 = std(M, 1, 3);
output(4, :) = reshape(m1', 1, []);
output(5, :) = reshape(m2', 1, []);

dlmwrite(filename, output', 'delimiter', ',', '-append');

% REVERSE CHAMELEON

paramA = [5, 9, 13];
motifs = [2, 3, 6];
logics = [2, 3, 1];
M = zeros(N, T, 9);
output = [2*ones(1,N*T); node; time; zeros(2,N*T)];
 
for i = 1:3
    for j = 1:3
        M(:,:,(i - 1)*3 + j) = D.noise_000{motifs(i), logics(j)}{paramA(j), 9};
    end
end

m1 = mean(M, 3);
m2 = std(M, 1, 3);
output(4, :) = reshape(m1', 1, []);
output(5, :) = reshape(m2', 1, []);

dlmwrite(filename, output', 'delimiter', ',', '-append');

end