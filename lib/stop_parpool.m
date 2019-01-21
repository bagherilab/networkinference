% -------------------------------------------------------------------------
% STOP_PARPOOL stops a parallel pool and remove the temporary directory
% containing project files.
% -------------------------------------------------------------------------

function stop_parpool(dir)

% Delete parallel pool object.
p = gcp('nocreate');
delete(p);

% Remove associated temporary file directory.
rmdir(dir,'s');

end