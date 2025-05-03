function Path = getPath(str)
    % Get the parent folder of the current script's directory
    [parentFolder, ~, ~] = fileparts(pwd);
    
    % Construct the path to the folder named str within the parent folder
    Path = fullfile(parentFolder, str);

    % Ensure the str folder exists, create it if it doesn't
    if ~exist(Path, 'dir')
        mkdir(Path);
    end
end