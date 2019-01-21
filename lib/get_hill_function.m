% -------------------------------------------------------------------------
% GET_HILL_FUNCTION returns a function handle to the single input
% activation hill equation.
% -------------------------------------------------------------------------
% INPUTS
%   k0 = affinity
%
% Returned function takes three inputs:
%	x = vector of concentrations of nodes
%	m = vector of node-to-node relationships
%   s = stimulation concentration
%
% Vector m is derived from the adjacency matric and represents which nodes
% do (1) or do not (0) activate the target node. Note that in this 
% system, we assume that there is only a single input to the target node
% other than possible external stimulation.
% -------------------------------------------------------------------------

function F = get_hill_function(k0)

F = @(x, m, s) ((k0*sum(x.*m) + s)/(1 + k0*sum(x.*m) + s));

end