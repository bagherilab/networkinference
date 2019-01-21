% -------------------------------------------------------------------------
% GET_STIM_FUNCTION returns a function handle for the selected stimulus
% input shape as (1) step, (2) ramp up, (3) ramp down, (4) double step up,
% (5) double step down, or (6) double pulse
% -------------------------------------------------------------------------
% INPUTS
%   INPUT = input selection
%   tS = time at which step occurs
%   v = max value
%
% Returned function takes one input:
%	t = time
% -------------------------------------------------------------------------

function F = get_stim_function(INPUT, tS, v)

f = @(t) (1 - heaviside(t - tS));

switch INPUT
    case 1 % step
        F = @(t) v*f(t);
    case 2 % ramp up
        F = @(t) v*f(t).*(1/tS*t);
    case 3 % ramp down
        F = @(t) v*f(t).*(1 - 1/tS.*t);
    case 4 % step up
        F = @(t) v*f(t).*(1 + heaviside(t - tS/2))/2;
    case 5 % step down
        F = @(t) v*f(t).*((1 - heaviside(t - tS/2)) + heaviside(t - tS/2)/2);
    case 6 % double pulse
        F = @(t) v*f(t).*((1 - heaviside(t - tS/3)) + heaviside(t - 2*tS/3));
    case 7 % triple pulse
        F = @(t) v*f(t).*((1 - heaviside(t - tS/5)) + heaviside(t - 2*tS/5) + ...
                (1 - heaviside(t - 3*tS/5)) + heaviside(t - 4*tS/5) - 1);
end

end