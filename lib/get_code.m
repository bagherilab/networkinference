% -------------------------------------------------------------------------
% GET_CODE return the string code and indicies for a given case.
%
% CASE 1
%   Code in the form M#L#S# and indices for the motif structure,
%   logic gate, and stimulus type given an INDEX ranging from 1 to
%   nM*nL*nS where nM is the number of motifs, nL is the number of logic
%   gates and nS is the number of stimulus types.
%
% CASE 2
%   Code in the form N#S#A# and indicies for the noise percent, stimulus
%   type, and parameter A value given an INDEX ranging from 1 to nN*nS*nP
%   where nN is the number of noise types, nS is the number of stimulus
%   types, and nP is the number of parameters.
% -------------------------------------------------------------------------

function [code, iA, iB, iC] = get_code(INDEX, CASE)

SETTINGS = LM_SETTINGS();

switch CASE
    case 1
        nL = SETTINGS.nLogics; % number of logic gates
        nS = SETTINGS.nStims; % number of stimulus types
        iMotif = ceil(INDEX/nL/nS); % index of motif matrix
        iLogic = ceil((INDEX - (iMotif - 1)*nL*nS)/nS); % index of logic gate
        iStim = INDEX - (iMotif - 1)*nL*nS - (iLogic - 1)*nS; % index of stimulus type
        code = ['M' num2str(iMotif) 'L' num2str(iLogic) 'S' num2str(iStim)];
        iA = iMotif; iB = iLogic; iC = iStim;
    case 2
        nS = SETTINGS.nStims; % number of stimulus types
        nP = SETTINGS.nParams; % number of parameters
        iNoise = ceil(INDEX/nS/nP); % index of noise type
        iStim = ceil((INDEX - (iNoise - 1)*nS*nP)/nP); % index of stimulus type
        iParamA = INDEX - (iNoise - 1)*nS*nP - (iStim - 1)*nP; % index of parameter A
        code = ['N' num2str(iNoise) 'S' num2str(iStim) 'A' num2str(iParamA)];
        iA = iNoise; iB = iStim; iC = iParamA;
end

end