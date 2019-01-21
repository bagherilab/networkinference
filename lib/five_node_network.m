% -------------------------------------------------------------------------
% FIVE_NODE_NETWORK evaluates the ODEs representing a five node network
% with a three-node upstream module containing a fan-in and a two-node
% downstream module emanating from the fan-in.
% -------------------------------------------------------------------------
% Inputs:
%   t       = time
%   y       = concentration
%   m       = 5x5 motif matrix where rows are inputs and columns are targets
%   s       = vector indicating which nodes are affected by stimulus
%   a       = vector of degradation rates
%   gateFun = function handle to gating activation
%   hillFun = function handle to hill activation
%   stimFun = function handle to step function for stimulation
% -------------------------------------------------------------------------

function dydt = five_node_network(t, y, m, s, a, gateFun, hillFun, stimFun)

dydt = zeros(5,1);

for i = 1:5
    if i == 3
        dydt(i) = gateFun(y(1), y(2)) - a(i)*y(i);
    else
        stim = stimFun(t)*s(i);
        dydt(i) = hillFun(y, m(:,i), stim) - a(i)*y(i);
    end
end

end