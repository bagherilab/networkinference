% -------------------------------------------------------------------------
% LM_MERGE_DIFFERENCE calculates average in group and out group city block
% distance for each motif and logic.
% -------------------------------------------------------------------------

function LM_merge_difference(ALGORITHM, FROM, TO, EDGE)

SETTINGS = LM_SETTINGS();
nM = SETTINGS.nMotifs;
nL = SETTINGS.nLogics;
nN = SETTINGS.nNoises;
nT = SETTINGS.nTimeslices;
nP = SETTINGS.nParams;

all_iw = zeros(nP, nP, nM, nL, nN, nT);
all_ers = zeros(nP, nP, nM, nL, nN, nT);

for iMotif = 1:nM
    for iLogic = 1:nL
        mls = ['M' num2str(iMotif) 'L' num2str(iLogic) 'S1'];
        D = load(['Analysis_' ALGORITHM '_' mls]);

        for iNoise = 1:nN
            noise = SETTINGS.noiseNames{iNoise};

            for iSlice = 1:nT
                for iParamA = 1:nP
                    for iParamB = 1:nP
                        all_iw(iParamA, iParamB, iMotif, iLogic, iNoise, iSlice) = ...
                            D.(noise).IW.mean{iParamA, iParamB, iSlice}(FROM, TO);
                        all_ers(iParamA, iParamB, iMotif, iLogic, iNoise, iSlice) = ...
                            D.(noise).ERS.mean{iParamA, iParamB, iSlice}(FROM, TO);
                    end
                end
            end
        end
    end
end

filename = sprintf('%s_%s_DIFF.csv', ALGORITHM, EDGE);
fid = fopen(filename, 'wt');
fprintf(fid, 'type,code,noise,slice,IW_in,ERS_in,IW_out,ERS_out\n');

diff_iw = zeros(nM*nL);
diff_ers = zeros(nM*nL);

for iNoise = 1:nN
    for iSlice = 1:nT
        IW = all_iw(:,:,:,:,iNoise,iSlice);
        ERS = all_ers(:,:,:,:,iNoise,iSlice);
        
        for iM1 = 1:nM
            for iM2 = 1:nM
                for iL1 = 1:nL
                    for iL2 = 1:nL
                        i1 = to_ind(iM1, iL1);
                        i2 = to_ind(iM2, iL2);
                        
                        diff_iw(i1, i2) = diff(IW(:,:,iM1,iL1), IW(:,:,iM2,iL2));
                        diff_ers(i1, i2) = diff(ERS(:,:,iM1,iL1), ERS(:,:,iM2,iL2));
                    end
                end
            end
        end
        
        for iM = 1:6
            in_inds = to_ind(iM, 1:6);
            n_in = length(in_inds);
            n_in_all = (n_in^2 - n_in)/2;
            D_in_iw = sum(sum(triu(diff_iw(in_inds,in_inds), 1)))/n_in_all;
            D_in_ers = sum(sum(triu(diff_ers(in_inds,in_inds), 1)))/n_in_all;
            
            out_inds = setdiff(1:36, in_inds);
            n_out = length(out_inds);
            n_out_all = (n_out^2 - n_out)/2;
            D_out_iw = sum(sum(triu(diff_iw(out_inds,out_inds), 1)))/n_out_all;
            D_out_ers = sum(sum(triu(diff_ers(out_inds,out_inds), 1)))/n_out_all;
            
            fprintf(fid, '%s,%d,%d,%d,%f,%f,%f,%f\n', 'M', iM, iNoise, iSlice, ...
                D_in_iw, D_in_ers, D_out_iw, D_out_ers);
        end
        
        
        for iL = 1:6
            in_inds = to_ind(1:6, iL);
            n_in = length(in_inds);
            n_in_all = (n_in^2 - n_in)/2;
            D_in_iw = sum(sum(triu(diff_iw(in_inds,in_inds), 1)))/n_in_all;
            D_in_ers = sum(sum(triu(diff_ers(in_inds,in_inds), 1)))/n_in_all;
            
            out_inds = setdiff(1:36, in_inds);
            n_out = length(out_inds);
            n_out_all = (n_out^2 - n_out)/2;
            D_out_iw = sum(sum(triu(diff_iw(out_inds,out_inds), 1)))/n_out_all;
            D_out_ers = sum(sum(triu(diff_ers(out_inds,out_inds), 1)))/n_out_all;
            
            fprintf(fid, '%s,%d,%d,%d,%f,%f,%f,%f\n', 'L', iL, iNoise, iSlice, ...
                D_in_iw, D_in_ers, D_out_iw, D_out_ers);
        end
    end
end

fclose(fid);

end

function i = to_ind(M,L)
i = (M - 1)*6 + L;
end

function [M, L] = from_ind(i)
M = ceil(i/6);
L = i - 6*(M - 1);
end

function d = diff(A, B) 
d = sum(sum(abs(A - B)));
end