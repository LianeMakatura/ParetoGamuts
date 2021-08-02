% use this to make sure the desired Python distribution is called via
% matlab

function fixEnvPath(desiredPythonDist)
    paths = strsplit(getenv('PATH'), ':');

    if ~any(contains(paths, desiredPythonDist)) %desired path not there
        setenv('PATH', [strcat(desiredPythonDist, ':'), getenv('PATH')]);
    end
end
