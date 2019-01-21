% -------------------------------------------------------------------------
% WRAP_BANJO_JAR acts a Matlab wrapper for the Java based BANJO algorithm.
% It converts the input data matrix into a tab-separated file, passes the
% file and a custom settings file (copied from a template), and runs the
% BANJO jar. The report file is parsed to get the adjacency matrix and then
% all three files are deleted.
% -------------------------------------------------------------------------

function weights = wrap_banjo_jar(data, code, template, path)

tic

% Number of variables.
n = size(data,2);

% File names.
dataFile = [path code '.data.txt'];
settingsFile = [path code '.settings.txt'];
reportFile = [path code '.report.txt'];

% Convert matrix data to required BANJO data format.
dlmwrite(dataFile, data, '\t');

% Create custom settings file using template.
fTemplate = fileread(template);
fTemplate = strrep(fTemplate, '[PATH]', path(1:end - 1));
fTemplate = strrep(fTemplate, '[CODE]', code);
fTemplate = strrep(fTemplate, '[NUMVARS]', num2str(n));
fTemplate = strrep(fTemplate, '[PARENT]', num2str(min(10, n - 1)));
fSettings = fopen(settingsFile,'w');
fprintf(fSettings, fTemplate);
fclose(fSettings);

% Load and run BANJO.
import edu.duke.cs.banjo.application.Banjo
Banjo.main(['settingsFile=' settingsFile]);

% Load and parse results file.
fReport = fileread(reportFile);
bestNetwork = regexp(fReport, 'Best network overall\n[-]+[\nA-z\s:0-9\-\.,]+[-]+', 'match');
network = regexp(bestNetwork, '[0-9]*[\s]{3}1:[\s]{3}[0-9\s]*\n', 'match');
network = network{1};

% Compile results into matrix.
weights = zeros(n);
for i = 1:n
    split = strsplit(network{i},'\n');
    split = strsplit(split{1}, ' '); % split along whitespace
    
    if str2double(split(1)) ~= (i - 1)
        error(['Missing row in report for ' code]);
    end
    
    iChild = str2double(split{1}) + 1; % get index of child
    nParents = str2double(split{3});
    for j = 1:nParents
        iParent = str2double(split{3 + j}) + 1; % get index of parent
        weights(iParent, iChild) = 1; % add to weight matrix
    end
end

delete(dataFile);
delete(settingsFile);
delete(reportFile);

toc

end