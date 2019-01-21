% -------------------------------------------------------------------------
% GET_GATE_FUNCTION returns a function handle to the selected double input
% logic gate activation equation.
% -------------------------------------------------------------------------
% INPUTS
%   LOGIC = gate selection
%   ki    = affinity of node i
%   kj    = affinity of node j
%
% Returned functions take two inputs:
%	xi = concentration of node i
%   xj = concentration of node j
% -------------------------------------------------------------------------

function F = get_gate_function(LOGIC, ki, kj)

switch LOGIC
    case 'AND'
        F = @(xi, xj) ((ki*kj*xi*xj)/(1 + ki*xi + kj*xj + ki*kj*xi*xj));
    case 'OR'
        F = @(xi, xj) ((ki*xi + kj*xj + ki*kj*xi*xj)/(1 + ki*xi + kj*xj + ki*kj*xi*xj));
    case 'SUM'
        F = @(xi, xj) (((ki*xi)/(1 + ki*xi)) + ((kj*xj)/(1 + kj*xj)))/2;
    case 'NAND'
        F = @(xi, xj) 1 - ((ki*kj*xi*xj)/(1 + ki*xi + kj*xj + ki*kj*xi*xj));
    case 'NOR'
        F = @(xi, xj) 1 - ((ki*xi + kj*xj + ki*kj*xi*xj)/(1 + ki*xi + kj*xj + ki*kj*xi*xj));
    case 'SUB'
        F = @(xi, xj) 1 - (((ki*xi)/(1 + ki*xi)) + ((kj*xj)/(1 + kj*xj)))/2;
end

end