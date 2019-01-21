% -------------------------------------------------------------------------
% ADD_NOISE adds gaussian noise with selected percent to matrix A.
% -------------------------------------------------------------------------
% INPUTS
%   A       = matrix of any size
%   percent = value between 0 and 1 indicating amount of noise
%
% We assume that standard gaussian has range of +/- 3 standard deviations
% so divide percent by three to scale noise to selected percent.
% -------------------------------------------------------------------------

function B = add_noise(A, percent)

noise = randn(size(A)); % generate random numbers from gaussian
scaled_noise = percent/3*noise; % scale noise to between +/- percent
B = A + A.*scaled_noise; % add noise to data

end