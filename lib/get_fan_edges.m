% -------------------------------------------------------------------------
% GET_FAN_EDGES gets edges to and from parents of a fan in node.
% -------------------------------------------------------------------------
% INPUTS
%   A       = adjancency matrix of any size
%   S       = vector of nodes that are stimulated
% OUTPUT
%   E       = cell array with field T# where # is the target of the fan in
%             with each entry being a matrix of edges
%   targets = vector of fan in target nodes
% -------------------------------------------------------------------------

function [E, targets] = get_fan_edges(A, S)

G = digraph(A);
N = length(A);
NODES = 1:N;
targets = [];
pars = [];
chils = [];
leafs = [];

% Find fan in target nodes.
for i = 1:N
    parents = A(:,i);

    if sum(parents) > 1
        children = A(i,:);
        c = find(children == 1);
        leaf = 0;
        
        for j = 1:length(c)
            cofc = A(c(j), :);
            
            if sum(cofc) == 0
                leaf = leaf + 1;
            end
        end

        targets = [targets i];
        pars = [pars sum(parents)];
        chils = [chils length(c)];
        leafs = [leafs leaf];
    end
end

% For each target node, find edges from parents.
E = struct();

for iT = 1:length(targets)
    all_edges = [];
    t = targets(iT);
    par = pars(iT);
    chil = chils(iT);
    leaf = leafs(iT);
    
    % Parents of fan in.
    parents = NODES(A(:,t) == 1);
    
    % Stimulated parents.
    stim = intersect(S, parents);
    
    % Edges to or from a parent node.
    for s = stim
        edges = [];
        
        for p = parents
            edges_from = NODES(A(p,:) == 1);
            edges_to = NODES(A(:,p) == 1);
            
            for ef = edges_from
                reachable = dfsearch(G, s);
                edges = [edges; 0 s p p ef all(ismember([p ef], reachable)) par chil leaf];
            end
            
            for et = edges_to
                reachable = dfsearch(G, s);
                edges = [edges; 1 s p et p all(ismember([et p], reachable)) par chil leaf];
            end
        end
        
        all_edges = cat(3, all_edges, edges);
    end
    
    E.(['T' num2str(t)]) = all_edges;
end

end