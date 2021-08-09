function data = readFromJSON(fileName)
    str = fileread(fileName); % dedicated for reading files as text
    data = jsondecode(str); % Using the jsondecode function to parse JSON from string
end

