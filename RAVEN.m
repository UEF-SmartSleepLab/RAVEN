function RAVEN(Fs,Fs_cap,filePath, filename, channels_of_interest,reference_signals, channels_of_interest_cap ,reference_signals_cap, microevents, save_path, mode, format, Thresholds)

%% (0) Ensure that all the necessary information is defined

fprintf('Analysis started.');
fprintf('\n');

if iscell(microevents) && any(ismember(microevents, {'spindle', 'kcomplex', 'delta_wave'})) && isempty(channels_of_interest{1})
    error('Please, ensure that you defined the channels of interest.');
end

if iscell(microevents) && any(ismember(microevents, {'cap'})) && isempty(channels_of_interest_cap{1})
    error('Please, ensure that you defined the channels of interest for CAP analysis.');
end

if isempty(save_path)
    error('Please, define the "Save Path"');
end


if iscell(microevents) && any(ismember(microevents, {'spindle', 'kcomplex', 'delta_wave', 'cap'})) && isempty(reference_signals{1}) && isempty(reference_signals_cap{1})
    reference = 'Yes';
elseif iscell(microevents) && any(ismember(microevents, {'spindle', 'kcomplex', 'delta_wave', 'cap'})) && ~isempty(reference_signals{1}) && ~isempty(reference_signals_cap{1})
    reference = 'No';
elseif iscell(microevents) && any(ismember(microevents, {'spindle', 'kcomplex', 'delta_wave'})) && isempty(reference_signals{1}) && isempty(channels_of_interest_cap{1})
    reference = 'Yes';
elseif iscell(microevents) && any(ismember(microevents, {'spindle', 'kcomplex', 'delta_wave'})) && ~isempty(reference_signals{1}) && isempty(channels_of_interest_cap{1})
    reference = 'No';
elseif iscell(microevents) && any(ismember(microevents, {'cap'})) && isempty(reference_signals_cap{1}) && isempty(channels_of_interest{1})
    reference = 'Yes';
elseif iscell(microevents) && any(ismember(microevents, {'cap'})) && ~isempty(reference_signals_cap{1}) && isempty(channels_of_interest{1})
    reference = 'No';
end



%% (1) Create folder(s) at save path for data(s)

if strcmp(mode, 'Group Analysis')
    % Get time of current run to avoid overriting existing analyse
    dt = datetime('now');

    % Extract components
    hourValue = hour(dt);
    minuteValue = minute(dt);
    secondValue = round(second(dt));

    % Combine into a formatted string and save in a variable
    timeVar = sprintf('%02d:%02d:%02d', hourValue, minuteValue, secondValue);
    timeVar = strrep(timeVar,':','_');

    % Find all files from given directory
    if strcmp(format, 'EDF')
        files = dir(fullfile(filePath, '*.edf'));
        % File names including .edf
        filenames = {files.name}';
        % Variable for file names without .edf
        fileNames = cell(length(filenames),1);

        % Initialize variables to store file names
        for i = 1 : length(filenames)
            [~, fileName, ~] = fileparts(filenames{i});

            % Check if the first character is a digit
            if isstrprop(fileName(1), 'digit')
                % Remove special characters but keep the initial digit
                validName = matlab.lang.makeValidName(fileName);
                fileNames{i} = validName(2:end);
            else
                fileNames{i} = matlab.lang.makeValidName(fileName); % Make valid name
            end
        end

        % Add the desired text "_EEG_analysis" to the name and create
        newFolderNames = cell(length(fileNames),1);
        for i = 1 : length(fileNames)
            newFolderNames{i} = [fileNames{i}, '_EEG_analysis_', timeVar];
        end

        % Initialise variable to save the saveFolderPaths
        saveFolderPaths = cell(length(newFolderNames),1);
        % Create a new folder path using save path and the modified name
        for i = 1 : length(newFolderNames)
            saveFolderPath = fullfile(save_path, newFolderNames{i});
            saveFolderPaths{i} = saveFolderPath;
            % Create the folder if it doesn't exist
            if ~exist(saveFolderPath, 'dir')
                mkdir(saveFolderPath);
            end
        end
    elseif strcmp(format,'SLF')
        allItems = dir(filePath);
        % Filter out the folders
        folders = allItems([allItems.isdir]);
        % Remove the '.' and '..' entries
        folders = folders(~ismember({folders.name}, {'.', '..'}));
        % Pick up the file names / ID
        filenames = {folders.name}';
        % Add the desired text "_EEG_analysis" to the name and create
        newFolderNames = cell(length(filenames),1);
        for i = 1 : length(filenames)
            newFolderNames{i} = [filenames{i}, '_EEG_analysis_', timeVar];
        end

        % Initialise variable to save the saveFolderPaths
        saveFolderPaths = cell(length(newFolderNames),1);

        % Create a new folder path using save path and the modified name
        for i = 1 : length(newFolderNames)
            saveFolderPath = fullfile(save_path, newFolderNames{i});
            saveFolderPaths{i} = saveFolderPath;
            % Create the folder if it doesn't exist
            if ~exist(saveFolderPath, 'dir')
                mkdir(saveFolderPath);
            end
        end
    end
elseif strcmp(mode, 'Individual Focus')

     % Extract the base filename without extension
     [~, fileName, ~] = fileparts(filename);
     % Check if the first character is a digit
     if isstrprop(fileName(1), 'digit')
         % Remove special characters but keep the initial digit
         validName = matlab.lang.makeValidName(fileName);
         fileName = validName(2:end);
     else
         fileName = matlab.lang.makeValidName(fileName); % Make valid name
     end

     % Get time of current run to avoid overwriting existing analyse
     dt = datetime('now');
     % Extract components
     hourValue = hour(dt);
     minuteValue = minute(dt);
     secondValue = round(second(dt));
     % Combine into a formatted string and save in a variable
     timeVar = sprintf('%02d:%02d:%02d', hourValue, minuteValue, secondValue);
     timeVar = strrep(timeVar,':','_');

     % Add the desired text "_EEG_analysis" to the name
     newFolderName = [fileName, '_EEG_analysis_', timeVar];
     % Create the new folder path using save_path and the modified filename
     saveFolderPath = fullfile(save_path, newFolderName);

     % Create the folder if it doesn't exist
     if ~exist(saveFolderPath, 'dir')
         mkdir(saveFolderPath);
     end
end



%% Analysis
% Find the signals of interest from given file. Filter and downsample
% the signals for analysis.

if strcmp(mode, 'Individual Focus')
    filenames = {filename};
    saveFolderPaths = {saveFolderPath};
end

% Create a waitbar
hWaitbar = waitbar(0,'Analysing patients...');

    for i = 1 : length(filenames)
        % Update waitbar for current patient
        waitbar(i / length(filenames), hWaitbar, sprintf('Analysing patient %d out of %d', i, length(filenames)));

        filename = filenames{i};
        [~,fileName,~] = fileparts(filename);
        % Check if the first character is a digit
        if isstrprop(fileName(1), 'digit')
            % Remove special characters but keep the initial digit
            validName = matlab.lang.makeValidName(fileName);
            fileName = validName(2:end);
        else
            fileName = matlab.lang.makeValidName(fileName); % Make valid name
        end
        saveFolderPath = saveFolderPaths{i};


        % (2) Signal formatting
        % Find the signals of interest from given file. Filter and
        % downsample the signals for analysis.
        if strcmp(format, 'EDF')
            mat_signal_formatted = edf_mat_converter(filePath,filename);
        elseif strcmp(format, 'SLF')
            mat_signal_formatted = slf_mat_converter(filePath,filename,channels_of_interest, channels_of_interest_cap);
        end

        % Pick up the raw EEG signals according to the channels of interest
        % for artifact identification
        % Initialize a struct to hold the extracted signals and sampling frequencies
        rawEEGSignals = struct();

        if isempty(channels_of_interest{1})
            channels_of_interest = channels_of_interest_cap;
        end
    
        % Find the index for each channel and extract data
        for k = 1:length(channels_of_interest)
            target_channel = channels_of_interest{k};
    
            % Find the index of the target channel
            channel_names = mat_signal_formatted(1, :); % Extract the channel names row
            channel_index = find(strcmp(channel_names, target_channel));
    
            % Check if the channel was found
            if isempty(channel_index)
                warning('Channel %s not found in the signal variable.', target_channel);
                % Assign empty values if the channel is not found
                rawEEGSignals.(target_channel).signal = [];
            else
                % Extract the signal and sampling frequency
                extracted_signal = mat_signal_formatted{2, channel_index};
                originalFs = mat_signal_formatted{3, channel_index};
                extracted_signal = resample(extracted_signal, Fs, originalFs);
    
                % Assign the extracted signal and sampling frequency to the struct
                rawEEGSignals.(target_channel).signal = extracted_signal;
            end
        end



        if isempty(channels_of_interest{1})
            disp('Only CAP analysis is performed.');
        else
            mat_signal_filtered = EEG_formatting(Fs,mat_signal_formatted,channels_of_interest,reference_signals,reference);
        end



        % (3) K-complex detection for channels of interest
        if iscell(microevents) && any(strcmp('kcomplex', microevents))
            Detect_Kcomplex(Fs, mat_signal_filtered, fileName, saveFolderPath,Thresholds);
        else
            fprintf('K-complex detection was not included in the analysis');
            fprintf('\n');
        end



        % (4) Spindle detection for channels of interest
        if iscell(microevents) && any(strcmp('spindle', microevents))
            Detect_Spindles(Fs, mat_signal_filtered, fileName, saveFolderPath,Thresholds);
        else
            fprintf('Spindle detection was not included in the analysis');
            fprintf('\n');
        end



        % (5) Delta waves detection (slow wave sleep) for channels of interest
        if iscell(microevents) && any(strcmp('delta_wave', microevents))
            Detect_delta_waves(Fs, mat_signal_filtered, fileName, saveFolderPath, rawEEGSignals,Thresholds);
        else
            fprintf('Delta wave detection was not included in the analysis');
            fprintf('\n');
        end



        % (6) CAP detection for channels of interest
        if iscell(microevents) && any(strcmp('cap', microevents))
            fprintf('Signal formatting started for CAP detection');
            fprintf('\n');
            mat_signal_filtered = EEG_formatting(Fs_cap,mat_signal_formatted,channels_of_interest_cap,reference_signals_cap,reference);
            Detect_CAP(Fs_cap, mat_signal_filtered, fileName, saveFolderPath,Thresholds);
        else
            fprintf('CAP detection was not included in the analysis');
            fprintf('\n');
        end
    end

    close(hWaitbar);
    fprintf('EEG microstructure analysis completed for all patients.');
    fprintf('\n');

end








%%-------------------------------------------------------------------------
% ALL FUNCTIONS

%%-------------------------------------------------------------------------
% EDF CONVERSION TO MAT
%%-------------------------------------------------------------------------



function mat_signal = edf_mat_converter(Path,filename)

    fprintf('\n');
    fprintf('Picking up signals from .edf file');
    fprintf('\n');

    % Get a list of all files in the folder
    edfFiles = dir(fullfile(Path, '*.edf'));

    % Initialize variable to store the full path of the selected file
    edf_file = '';

    % Loop through the list of files to find the matching filename
    for k = 1:length(edfFiles)
        if strcmp(edfFiles(k).name, filename)
            edf_file = fullfile(Path, edfFiles(k).name);
            break; % Exit the loop once the file is found
        end
    end

    % Check if the file was found
    if isempty(edf_file)
        error('File with the name "%s" was not found in the directory "%s".', filename, Path);
    end

    [~, baseFileName, ~] = fileparts(edf_file);

    % Replace all space with "_" in baseFileName for eval function to work
    % properly
    % Check if the first character is a digit
    if isstrprop(baseFileName(1), 'digit')
        % Remove special characters but keep the initial digit
        validName = matlab.lang.makeValidName(baseFileName);
        baseFileName = validName(2:end);
    else
        baseFileName = matlab.lang.makeValidName(baseFileName); % Make valid name
    end


    signalData = edfread(edf_file);

    times = signalData.("Record Time");
    stepSize = diff(times);

    % Calculate the average time difference
    averageStepSize = seconds(mean(stepSize));

    if averageStepSize > 1
        % Function for changing the sampling frequency of the data to
        % be samples/ 1 second
        signalData = resample_data(signalData, averageStepSize);
    end


    % Make the signals to be samples/second and save the data in mat
    % files with sampling frequencies
    signalData_processed = concatenate_data_table(signalData);
    mat_signal = signalData_processed;
    varName = ['signalData_' baseFileName];

    % Create a variable with the dynamic name
    eval([varName ' = mat_signal;']);

end











function resampled_sigalData = resample_data(data, step)
% Pick up the column labels
columnLabels = data.Properties.VariableNames;

% Count the amount of signals (columns)
numSignals = width(data);

% Count the amount of time steps
numSteps = height(data);

newData = cell(numSteps*step, numSignals);

for i = 1 : numSignals
    % Pick up the current columns signal
    currentSignal = data{:,i};

    % Each cell has to be now divided by the step
    % Create a new signal column to which include the divided data
    newSignal = cell(numSteps*step, 1);

    pos = 1;

    for j = 1 : numSteps
        cellData = currentSignal{j};

        % Calculate the size of each cell
        cellSize = ceil(length(cellData) / step);

        % Determine the starting indices for each cell
        startIndex = [1, cumsum(repmat(cellSize, 1, step-1)) + 1];
        % Ensure the last cell extends to the end of the array
        startIndex(end) = min(startIndex(end), length(cellData));

        % Create the cell array
        cellArray = mat2cell(cellData, diff([startIndex, length(cellData)+1]), 1);

         % Append the cellArray to newSignal
        for k = 1 : step
            newSignal{pos} = cellArray{k};
            pos = pos + 1;
        end

    end

    newData(:,i) = newSignal;


end

% Create a timetable from newData
newTimetable = timetable();


for i = 1:length(columnLabels)
    % Get the current column label
    label = columnLabels{i};

    % Add the data from newData(:, i) to the timetable
    newTimetable.(label) = newData(:, i);
end

% Define the time vector (seconds from 1 to the end of newData)
timeVector = (1:size(newData, 1))';
newTimetable.Time = seconds(timeVector);

resampled_sigalData = newTimetable;

end











function concatSignals = concatenate_data_table(signalData)

    % Get fieldnames of signalData
    variables = signalData.Properties.VariableNames;

    % Initialize concatenated signals cell array
    concatSignals = cell(3, numel(variables));  % Preallocate for efficiency

    % Loop through each variable (column) in signalData
    for i = 1:numel(variables)

        if iscell(signalData.(variables{i}))
            % Extract current signals
            signals_cell = signalData.(variables{i});

            % Fs calculation
            samplingFreq = length(signals_cell{1});

            % Concatenate signals into a single vector
            concatenated_signal = vertcat(signals_cell{:});

            % Store the variable name in the first row
            concatSignals{1, i} = variables{i};

            % Store the concatenated signal in the second row
            concatSignals{2, i} = concatenated_signal;

            % Store the Fs to the third row
            concatSignals{3,i} = samplingFreq;
        else
            concatSignals{1, i} = variables{i};
            concatSignals{2, i} = signalData.(variables{i});
            concatSignals{3,i} = [];

        end
    end

end






%%-------------------------------------------------------------------------
% SLF TO MAT CONVERSION
%%-------------------------------------------------------------------------

function mat_signal = slf_mat_converter(Path,fileName,channels_micro, channels_cap)

    % Combine all the channels of interest.
    channels = union(channels_micro, channels_cap);

    % Remove empty cells from channels.
    non_empty_indices = ~cellfun('isempty', channels);
    channels = channels(non_empty_indices);

    % Create variable including all the signal names, signals, and sampling
    % frequencys
    signalData = cell(3,numel(channels));

    % Get a list of all items in the directory
    items = dir(Path);

    % Filter to find directories with the specified name
    matchingFolder = items([items.isdir] & strcmp({items.name}, fileName));

    % Check if the folder exists
    if ~isempty(matchingFolder)
        % Get the full path of the matching folder
        specificFolderPath = fullfile(Path, matchingFolder(1).name);

        % List all subfolders inside the specific folder
        subItems = dir(specificFolderPath);

        % Filter out '.' and '..' and ensure they are directories
        subFolders = subItems([subItems.isdir] & ~ismember({subItems.name}, {'.', '..'}));

        % Get the names of the subfolders
        subFolderNames = {subFolders.name};

        % Find subfolders matching the desired channels
        matchingSubFolders = subFolderNames(contains(subFolderNames, channels));

        % Identify missing channels
        missingChannels = setdiff(channels, matchingSubFolders);

        if ~isempty(missingChannels)
            fprintf('The following channels are incorrect: %s\n', strjoin(missingChannels, ', '));
            error('Channels that found in directory under patient ID: %s\n', strjoin(subFolderNames, ', '));
        end


         % Loop through matching subfolders
        for i = 1:length(matchingSubFolders)
            % Get the current subfolder name
            subFolderName = matchingSubFolders{i};

            % Build the full path to the current subfolder
            currentSubFolderPath = fullfile(specificFolderPath, subFolderName);

            % Get the list of all .npy files from current dir. There should
            % be only one .npy file for each channel.
            npyFiles = dir(fullfile(currentSubFolderPath,'*.npy'));
            jsonFiles = dir(fullfile(currentSubFolderPath,'*.json'));

            if length(npyFiles) ~= 1
                error('The %s folder does not contain .npy file OR there are more than 1 .npy file.', subFolderName);
            end

            if length(jsonFiles) ~= 1
                error('The %s folder does not contain .json file for attributes OR there are more than 1 .json file.', subFolderName)
            end

            % Get the name of .json file for attributes
            attriName = jsonFiles.name;

            % Build the full path to the "data.npy" file
            dataFilePath = fullfile(currentSubFolderPath, npyFiles.name);
            attributeFilePath = fullfile(currentSubFolderPath, attriName);
            attributesRead = fileread(attributeFilePath);

            % Check if the data file exists
            if exist(dataFilePath, 'file')
                % Load the .npy file (requires npy-matlab library)
                data = readNPY(dataFilePath);

                % Get the correct channel
                channel = channels{i};

                % Get the attributes
                attributes = jsondecode(attributesRead);

                % Extract the sampling rate for signals
                samplingRate = attributes.sampling_rate;

                % Create a variable name dynamically based on the folder name
                variableName = ['signal_' subFolderName];

                % Assign the data to the dynamically created variable
                eval([variableName ' = data;']);

                % Append the information to signalData
                signalData{1,i} = channel;
                signalData{2,i} = data;
                signalData{3,i} = samplingRate;

                % Display success message
                fprintf('Successfully loaded data from "%s".\n', subFolderName);
            else
                fprintf('File "%s" not found in folder "%s".\n', npyFiles.name, subFolderName);
            end
        end
    else
        fprintf('No folder with the name "%s" found in "%s".\n', fileName, Path);
    end

    mat_signal = signalData;

end


function data = readNPY(filename)
% Function to read NPY files into matlab.

[shape, dataType, fortranOrder, littleEndian, totalHeaderLength, ~] = readNPYheader(filename);

if littleEndian
    fid = fopen(filename, 'r', 'l');
else
    fid = fopen(filename, 'r', 'b');
end

try

    [~] = fread(fid, totalHeaderLength, 'uint8');

    % read the data
    data = fread(fid, prod(shape), [dataType '=>' dataType]);

    if length(shape)>1 && ~fortranOrder
        data = reshape(data, shape(end:-1:1));
        data = permute(data, [length(shape):-1:1]);
    elseif length(shape)>1
        data = reshape(data, shape);
    end

    fclose(fid);

catch me
    fclose(fid);
    rethrow(me);
end
end

function [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename)

fid = fopen(filename);

% verify that the file exists
if (fid == -1)
    if ~isempty(dir(filename))
        error('Permission denied: %s', filename);
    else
        error('File not found: %s', filename);
    end
end

try
    
    dtypesMatlab = {'uint8','uint16','uint32','uint64','int8','int16','int32','int64','single','double', 'logical'};
    dtypesNPY = {'u1', 'u2', 'u4', 'u8', 'i1', 'i2', 'i4', 'i8', 'f4', 'f8', 'b1'};
    
    
    magicString = fread(fid, [1 6], 'uint8=>uint8');
    
    if ~all(magicString == [147,78,85,77,80,89])
        error('readNPY:NotNUMPYFile', 'Error: This file does not appear to be NUMPY format based on the header.');
    end
    
    majorVersion = fread(fid, [1 1], 'uint8=>uint8');
    minorVersion = fread(fid, [1 1], 'uint8=>uint8');
    
    npyVersion = [majorVersion minorVersion];
    
    headerLength = fread(fid, [1 1], 'uint16=>uint16');
    
    totalHeaderLength = 10+headerLength;
    
    arrayFormat = fread(fid, [1 headerLength], 'char=>char');
    
    
    r = regexp(arrayFormat, '''descr''\s*:\s*''(.*?)''', 'tokens');
    if isempty(r)
        error('Couldn''t parse array format: "%s"', arrayFormat);
    end
    dtNPY = r{1}{1};    
    
    littleEndian = ~strcmp(dtNPY(1), '>');
    
    dataType = dtypesMatlab{strcmp(dtNPY(2:3), dtypesNPY)};
        
    r = regexp(arrayFormat, '''fortran_order''\s*:\s*(\w+)', 'tokens');
    fortranOrder = strcmp(r{1}{1}, 'True');
    
    r = regexp(arrayFormat, '''shape''\s*:\s*\((.*?)\)', 'tokens');
    shapeStr = r{1}{1}; 
    arrayShape = str2num(shapeStr(shapeStr~='L'));

    
    fclose(fid);
    
catch me
    fclose(fid);
    rethrow(me);
end
end
%%------------------------------------------------------------------------
% FORMATTING THE SIGNALS
%%-------------------------------------------------------------------------




function mat_signal_filtered = EEG_formatting(fs_new, signals, channels_of_interest, reference_signals,reference)



    % Extract the original sampling frequency for channels of interest
    original_Fs = [];

    for i = 1:length(channels_of_interest)
        % Find the colum index of the signal names
        idx = find(strcmp(signals(1,:), channels_of_interest{i}));

        % If the signal name exist in the cell array
        if ~isempty(idx)
            % Append the corresponding sampling frequency to original_Fs
            original_Fs = [original_Fs, signals{3, idx}];
        else
            singnalsInEdf = signals(1,:);
            % Display a message if the signal name is not found
            fprintf('Signal name %s is incorrect.\n', channels_of_interest{i});
            fprintf('Signal names found in edf: %s\n', strjoin(singnalsInEdf, ', '));
            error('Check the input signal names. If the signal names are correct, please ensure there is no space after the commas.');
        end
    end

    % Check that the sampling frequencies are identical.
    if isscalar(unique(original_Fs))
        disp('All sampling frequencies are identical.');
        original_Fs = original_Fs(1);
    else
        error('Sampling frequencies are not identical. Channels of interest must have identical sampling frequencies');
    end



    %Design filters
    hpFilt = designfilt('highpassiir','FilterOrder',3, ...
        'PassbandFrequency',0.1,'PassbandRipple',0.001, ...
        'SampleRate',original_Fs,'DesignMethod','cheby1');
    lpFilt = designfilt('lowpassiir','FilterOrder',20, ...
        'PassbandFrequency',30,'PassbandRipple',0.001, ...
        'SampleRate',original_Fs,'DesignMethod','cheby1');


    % Initialize a struct to hold the extracted signals and sampling frequencies
    signals_struct = struct();

    % Find the index for each channel and extract data
    for i = 1:length(channels_of_interest)
        target_channel = channels_of_interest{i};

        % Find the index of the target channel
        channel_names = signals(1, :); % Extract the channel names row
        channel_index = find(strcmp(channel_names, target_channel));

        % Check if the channel was found
        if isempty(channel_index)
            warning('Channel %s not found in the signal variable.', target_channel);
            % Assign empty values if the channel is not found
            signals_struct.(target_channel).signal = [];
            signals_struct.(target_channel).sampling_frequency = [];
        else
            % Extract the signal and sampling frequency
            extracted_signal = signals{2, channel_index};
            extracted_sampling_frequency = signals{3, channel_index};

            % Assign the extracted signal and sampling frequency to the struct
            signals_struct.(target_channel).signal = extracted_signal;
            signals_struct.(target_channel).sampling_frequency = extracted_sampling_frequency;
        end
    end


    sampling_frequency = signals_struct.(cell2mat(channels_of_interest(1))).sampling_frequency;

    mat_signal_filtered = table();

        if strcmp(reference, 'No') && fs_new ~= sampling_frequency

            fprintf('\n');
            fprintf('subtracting the reference channels, filtering, and resampling signals');
            fprintf('\n');

            % Process reference signals
            for j = 1:length(reference_signals)
                % Split the reference signal description (e.g., 'C3-M2')
                ref_signal_desc = reference_signals{j};
                parts = strsplit(ref_signal_desc, '-');
                if length(parts) ~= 2
                    warning('Invalid reference signal format: %s', ref_signal_desc);
                    continue;
                end

                % Get the signals to be subtracted
                signal1 = parts{1};
                signal2 = parts{2};

                if isfield(signals_struct, signal1) && isfield(signals_struct, signal2)
                    signal1_data = signals_struct.(signal1).signal;
                    signal2_data = signals_struct.(signal2).signal;
                    reference_signal_data = signal1_data - signal2_data;

                    % Filter and resample the reference signal
                    current_filtered_signal = filter_and_resample_EEG(reference_signal_data, lpFilt, hpFilt, fs_new, sampling_frequency);

                    % Add the filtered reference signal to the table with the corresponding name
                    mat_signal_filtered.(signal1) = current_filtered_signal;
                else
                    error('Signals for reference %s are not available.', ref_signal_desc);
                end
            end

        else

            fprintf('\n');
            fprintf('Filtering and downsampling signals');
            fprintf('\n');

            for j = 1 : numel(fieldnames(signals_struct))
                % Get the current signal and its name
                current_signal = signals_struct.(cell2mat(channels_of_interest(j))).signal;
                current_signal_name = cell2mat(channels_of_interest(j));

                % Filter and resample the current signal
                current_filtered_signal = filter_and_resample_EEG(current_signal,lpFilt,hpFilt,fs_new,sampling_frequency);

                % Add the current filtered signal to the table with the corresponding name
                mat_signal_filtered.(current_signal_name) = current_filtered_signal;
            end

        end


end

function [resampled_channel] = filter_and_resample_EEG(channel,lpFilt,hpFilt,fs_new,sampling_frequency)
if isa(channel, 'single')
    channel = double(channel);
end
% High pass filter
channel = filtfilt(hpFilt,channel);
% Low pass filter
channel = filtfilt(lpFilt,channel);
%resample
resampled_channel = resample(channel,fs_new, sampling_frequency);
end









%%-------------------------------------------------------------------------
% K-COMPLEX DETECTION
%%-------------------------------------------------------------------------




function kComplexes = Detect_Kcomplex(Fs, mat_signal_filtered, filename, save_path, Thresholds)




channels_to_analyse = mat_signal_filtered.Properties.VariableNames;


% Default thresholds
thrs_power = Thresholds.Power_Kcomplex;           % Delta band
thrs_rms = Thresholds.RMS_Kcomplpex;              % Filtered signal (base frequency)
thrs_hilbert = Thresholds.Hilbert_Kcomplex;          % Delta band
thrs_cwt = Thresholds.CWT_Kcomplex;              % Filtered signal (base frequency)
thrs_wave_length = Thresholds.WaveLength_Kcomplex;    % Delta band (in seconds)

% Alternative annotation name
aasmEvent_name = 'UNSURE';


% Define column names
columnNames = channels_to_analyse;

% Create an empty cell array with specified column names
dataTableKcomplex = cell(1, numel(columnNames)); % Start with one row

% Assign column names to the first row of the cell array
dataTableKcomplex(1, :) = columnNames;

event = 'kComplex';



% Analyse the channels of interest for K-complexes

for channel = 1 : length(channels_to_analyse)

    chan = cell2mat(channels_to_analyse(channel));
    EEGbf = mat_signal_filtered.(chan);


    fprintf('\n');
    fprintf('Detecting K-complexes for channel %s', chan);
    fprintf('\n');


    % Window size and window step for caluculating powers in signal
    WinL = round(0.3 * Fs);  % Window length (0.3 seconds)
    WinS = floor(0.1 * Fs);  % Step size in samples (0.1 seconds)


    % Filter the data for delta band

    f_low = 0.1; % Lower cutoff frequency
    f_high = 5.5; % Upper cutoff frequency
    order = 2; % Filter order

    % Design the bandpass filter using the Butterworth filter for delta band
    [b, a] = butter(order, [f_low, f_high]/(Fs/2), 'bandpass');

    % Apply the bandpass filter to the signal
    EEG_delta = filter(b, a, EEGbf);

    delta_power = zeros(size(EEG_delta));

    padL = floor((0.1/2)*Fs); % Zero padding length for the window

    fprintf('Calculating delta frequency range power. This takes few seconds ');
    fprintf('\n');

    for i = 1:WinS:length(EEG_delta) - WinL + 1
        window = EEG_delta(i:i+WinL-1);
        if i > 1
            window = padarray(window,padL,'pre');
        end
        delta_power(i:i+WinL-1) = log10(sum(window.^2) / WinL);
    end






    WinL = 30*Fs;
    overlap = round(0.5*WinL);

    % Calculate the total number of windows
    total_windows = floor((length(EEGbf) - WinL) / (WinL - overlap)) + 1;

    % Initialize RMS values array for delta band
    envelope_values_delta = zeros(size(EEG_delta));

    % Loop through the signal and calculate Hilbert transfomr
    % for each window using envelope
    for i = 1:total_windows
        % Calculate the start and end indices of the current window
        start_idx = (i-1) * (WinL - overlap) + 1;
        end_idx = start_idx + WinL - 1;

        % Ensure the indices are within the signal length
        if end_idx > length(EEG_delta)
            end_idx = length(EEG_delta);
        end

        % Extract the windowed segment from the delta band signal
        window_segment = EEG_delta(start_idx:end_idx);

        % Calculate the envelope for the current window
        envelope = abs(hilbert(window_segment));

        % Assign the envelope values to the corresponding positions in
        % signal
        envelope_values_delta(start_idx:end_idx) = envelope;
    end




    WinL = 1*Fs;
    overlap = round(0.5*WinL);

    % Calculate the total number of windows
    total_windows = floor((length(EEGbf) - WinL) / (WinL - overlap)) + 1;

    % Initialize RMS values array for delta band
    rms_values_delta = zeros(size(EEG_delta));
    rms_bf = zeros(size(EEGbf));

    % Loop through the signal and calculate RMS for each window
    for i = 1:total_windows
        % Calculate the start and end indices of the current window
        start_idx = (i-1) * (WinL - overlap) + 1;
        end_idx = start_idx + WinL - 1;

        % Ensure the indices are within the signal length
        if end_idx > length(EEG_delta)
            end_idx = length(EEG_delta);
        end

        % Extract the windowed segment
        window_segment = EEG_delta(start_idx:end_idx);
        window_segment_bf = EEGbf(start_idx:end_idx);

        % Calculate the RMS value for the current window
        rms_value = sqrt(mean(window_segment .^ 2));
        rms_value_bf = sqrt(mean(window_segment_bf .^ 2));

        % Assign the RMS value to the corresponding positions
        rms_values_delta(start_idx:end_idx) = rms_value;
        rms_bf(start_idx:end_idx) = rms_value_bf;
    end










    % Detect kcomplex candidates with point system using thresholds. If all
    % thresholds are exceeded Kcomplex is most likely detected. kcomplex should
    % be at least 0.5 s in length.



    % Update the WinL and make the window overlap to be 50%
    WinL = round(30*Fs);
    overlap = round(0.5*WinL);

    % Calculate the total number of windows
    total_windows = floor((length(EEGbf) - WinL) / (WinL - overlap)) + 1;

    kComplexes = cell(total_windows, 1);

    h = waitbar(0, 'Detecting Kcomplexes over night...', 'Name', 'kcomlex detection');

    for i = 1:total_windows
        % Calculate start and end indices for the current window
        startIndex = 1 + (i - 1) * (WinL - overlap);
        endIndex = startIndex + WinL - 1;

        % Ensure endIndex does not exceed the length of the signal
        endIndex = min(endIndex, length(EEGbf));

        power_window = delta_power(startIndex:endIndex);

        % If there are no values over thresholds in both rms window and
        % envelope window, skip the whole iteration. In this scenario
        % the events will not have the needed points at the end in any
        % way. In other words, the result will be false true or false.

        rms_window = rms_values_delta(startIndex:endIndex);
        envelope_window = envelope_values_delta(startIndex:endIndex);

        if ~any(rms_window >= thrs_rms) || ~any(envelope_window >= thrs_hilbert)
            continue;
        end




        % CALCULATING THE CWT FROM STANDARD SIGNAL
        [cwtResult, frequencies] = cwt(EEGbf(startIndex:endIndex), 'amor', Fs);

        % Find regions where CWT values exceed the threshold
        highFreqIndices = frequencies >= 0.5 & frequencies <= 6;
        cwtHighFreq = abs(cwtResult(highFreqIndices, :));
        % Create masks for threshold exceeding regions
        mask = cwtHighFreq > thrs_cwt;



        % Find indices over the threshold
        indices_over_threshold = find(power_window >= thrs_power);

        % Check if there are any candidates
        if isempty(indices_over_threshold)
            % If no consecutive indices, break out of the loop and move to the next iteration
            continue;
        end



        % Find consecutive indices
        consecutive_indices = find_consecutive_kcomplexes(indices_over_threshold);
        % Check if there are any consecutive runs within the desired length range
        if isempty(consecutive_indices)
            % If no consecutive indices, break out of the loop and move to the next iteration
            continue;
        end



        % Map indices from power_window to kComplex_power for each consecutive run
        filtered_consecutive_indices = cell(size(consecutive_indices));

        for j = 1:length(consecutive_indices)
            % Map indices from power_window
            real_indices = consecutive_indices{j} + startIndex - 1;

            first_index = real_indices(1,1);

            % Calculate starting point in seconds
            starting_time_seconds = first_index / Fs;

            % Store filtered indices and starting time as a struct or tuple
            filtered_consecutive_indices{j} = struct('indices', real_indices, 'start_sec', starting_time_seconds);
        end





        event_indices = cell(size(filtered_consecutive_indices));

        % Pick up the indices to speed up the process
        for j = 1:length(filtered_consecutive_indices)
            event_indices{j} = filtered_consecutive_indices{j,1}.indices;
        end






        % Initialize points for the current window
        points_per_consecutive_part = zeros(size(event_indices));

        % CANDIDATE EVALUATION


        % Iterate over each consecutive run and check values for rms
        for j = 1:size(event_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = event_indices{j}(1);
            end_index = event_indices{j}(end);

            % Subset rms_values array between start_index and end_index
            subset_rms_values = rms_values_delta(start_index:end_index);
            subset_rms_bf = rms_bf(start_index:end_index);

            % Process the result based on the rms values
            if any(subset_rms_values >= thrs_rms)
                points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
            end

            % Process the result based on the rms values
            if any(subset_rms_bf >= 180)
                points_per_consecutive_part(j) = points_per_consecutive_part(j) - 1;
            end
         end







         % Iterate over each consecutive run and check values in Hilbert
         % trasform of kcomplex filtered signal
         for j = 1:size(filtered_consecutive_indices, 1)
             % Get start and end indices of the current consecutive run
             start_index = event_indices{j}(1);
             end_index = event_indices{j}(end);

             subset_envelope = envelope_values_delta(start_index:end_index);

             % Check if any value in envelope is over the threshold
             if any(subset_envelope >= thrs_hilbert)
                 points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
             end
         end



         % Iterate over each consecutive run and check the amplitude
         % differrence
         for j = 1:size(filtered_consecutive_indices, 1)
             % Get start and end indices of the current consecutive run
             start_index = event_indices{j}(1);
             end_index = event_indices{j}(end);

             signal_part = EEGbf(start_index:end_index);
             % Calculate the peak-to-peak value
             peak_to_peak = max(signal_part) - min(signal_part);
             threshold = 70;

             % Check if any value in envelope is over the threshold
             if peak_to_peak >= threshold
                 points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
             end
         end



        % Wave length check. In this part the code checks the wave
        % length of the signal in detected area.
        for j = 1:size(filtered_consecutive_indices, 1)
             % Get start and end indices of the current consecutive run
             start_index = event_indices{j}(1);
             end_index = event_indices{j}(end);

             % Extract the corresponding part from EEG_sigma
             delta_consecutive_part = EEG_delta(start_index:end_index,1);

             [~, locs] = findpeaks(delta_consecutive_part);
             if length(locs) > 1
                 period_samples = diff(locs); % Differences between consecutive peaks
                 periods = period_samples / Fs; % Convert to time periods

                 % Count the number of periods greater than the threshold
                 num_long_periods = sum(periods >= thrs_wave_length);

                 % Process the result based on the count of long periods
                 if num_long_periods <= 2 && num_long_periods ~= 0
                     points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
                 end
             end
        end





         %Iterate over each consecutive run and check values in CWT
         %coefficients of kcomplex frequency filtered signal
         for j = 1:size(filtered_consecutive_indices, 1)
             % Get start and end indices of the current consecutive run.
             % Map the indices to find the area in current CWT transform
             % for epoch.
             start_index = event_indices{j}(1) - startIndex + 1;
             end_index = event_indices{j}(end) - startIndex + 1;

             timeIndices = start_index:end_index;

             % Pick up the values from cwt results for candidate area
             candidate_areas = bwconncomp(mask(:,timeIndices));

             % Initialize a variable for the number of areas in the specified time range
             numAreasInRange = candidate_areas.NumObjects;

             % Check if averageOfValuesOverOne exceeds the threshold
             if numAreasInRange > 0
                 points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
             end
         end



         % Iterate over each consecutive run and check the zerocrossing
         for j = 1:size(filtered_consecutive_indices, 1)
             % Get start and end indices of the current consecutive run
             start_index = event_indices{j}(1);
             end_index = event_indices{j}(end);

             signal_part_zero_crossing = EEG_delta(start_index:end_index);

             zc = diff(sign(signal_part_zero_crossing)) ~= 0;
             zc_count = sum(zc);

             % Check the zero crossing value
             if zc_count >= 2
                 points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
             end

             % Subtract one point if zc_count is too high. This
             % indicates slow wave sleep or artifact.
             if zc_count >= 14
                 % Process the result based on the values in envelope
                 points_per_consecutive_part(j) = points_per_consecutive_part(j) - 1;
             end
         end







         % POINT CALCULATION

         % Check if any consecutive part has accumulated six points
         if any(points_per_consecutive_part == 6)
             % Find indices of consecutive parts with six points
             indices_five_points = find(points_per_consecutive_part == 6);

             % Create a cell array to store Kcomplex indices for the current window
             kComplex_candidates = cell(length(indices_five_points), 1);

             % Iterate over indices and store corresponding segments in current_kcomplex_indices
             for k = 1:length(indices_five_points)
                 current_consecutive_part = filtered_consecutive_indices{indices_five_points(k)};
                 kComplex_indices = current_consecutive_part.indices(1,1):current_consecutive_part.indices(end,1);
                 current_kcomplex_starting_time = current_consecutive_part.start_sec;
                 current_kcomplex_duration = length(kComplex_indices) / Fs;
                 kComplex_candidates{k} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', chan, 'sampling_frequency', Fs, ...
                                                    'start_sec', current_kcomplex_starting_time, 'duration', current_kcomplex_duration, ...
                                                    'indices', kComplex_indices);
             end

             % Store Kcomplex indices for the current window in
             % kcomplex_candidates
             kComplexes{i} = kComplex_candidates;
         end

         % Update waitbar
         waitbar(i/total_windows, h);



     end

     close(h);

     % Remove empty cells from kcomplex (indices)
     non_empty_indices = ~cellfun('isempty', kComplexes);
     kComplexes = kComplexes(non_empty_indices);





     % Connect those event areas that has overlap with each other due to
     % the overlap in analyzing windows

     for i = 1:size(kComplexes, 1) - 1

         current_indices = [];
         next_indices = [];

         for j = 1:numel(kComplexes{i})
             current_cell = kComplexes{i, 1}{j, 1}.indices;
             current_indices{j} = current_cell;
         end

         for k = 1:numel(kComplexes{i+1})
             next_cell = kComplexes{i+1, 1}{k, 1}.indices;
             next_indices{k} = next_cell;
         end

         for n = 1:numel(current_indices)

             % Check if the current cell is not empty
             if ~isempty(current_indices)
                     current_set = current_indices{n}; % Extract the set from the single set cell

                     % Check if the next cell has sets
                     if ~isempty(next_indices)  % Check if the next cell has at least one set
                         % Find the indices of sets in the next cell that include any index from the current set
                         idx_included = find(cellfun(@(x) any(ismember(current_set, x)), next_indices));

                         if ~isempty(idx_included)
                             % Combine the indices in cells that have overlapping indices
                             combined_indices = unique([current_set, next_indices{idx_included}]);

                             % Update the current cell with combined indices
                             kComplexes{i, 1}{j, 1}.indices = combined_indices;

                             % Remove all included sets from the next cell
                             kComplexes{i + 1}(idx_included) = [];
                             next_indices(idx_included) = [];
                         end
                     end
             else
                 % Continue the loop if the current cell is empty
                 continue;
             end
         end
     end

     % Remove empty cells from kComplexes
     non_empty_indices = ~cellfun('isempty', kComplexes);
     kComplexes = kComplexes(non_empty_indices);

     % Initialize a new cell array to store the resulting structs
     resulting_structs = {};

     % Loop through each cell of kComplexes
     for i = 1:numel(kComplexes)
         % Check if the current cell contains structs
         if iscell(kComplexes{i})
             % If there are structs, extract them individually
             for j = 1:numel(kComplexes{i})
                 resulting_structs{end+1} = kComplexes{i}{j};
             end
         else
             % If there's only one struct, extract it
             resulting_structs{end+1} = kComplexes{i};
         end
     end

     kComplexes = resulting_structs';




     % Connect those detected events that are closer than 0.5 second
     % from each other

     difference = 0.5 * Fs; % in seconds
     connectedEvents = cell(length(kComplexes),1);

     i = 1;

     while i <= length(kComplexes) - 1
         current_event_indices = kComplexes{i,1}.indices;
         current_event_last = kComplexes{i,1}.indices(end);
         current_event_start = kComplexes{i, 1}.start_sec;
         current_event_dur = kComplexes{i, 1}.duration;

         next_event_indices = kComplexes{i+1,1}.indices;
         next_event_first = kComplexes{i+1,1}.indices(1);

         if (next_event_first - current_event_last) <= difference
             combined_indices = current_event_indices(1):next_event_indices(end);
             combined_duration = length(combined_indices)/Fs;
             connectedEvents{i} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event,  'input_channel', chan, 'sampling_frequency', Fs, ...
                                         'indices', combined_indices, 'start_sec', current_event_start, ...
                                         'duration', combined_duration);
             i = i + 2;
         else
             connectedEvents{i} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', chan, 'sampling_frequency', Fs, ...
                                         'indices', current_event_indices, 'start_sec', current_event_start, ...
                                         'duration', current_event_dur);
             i = i +1;
         end
     end

     % Remove empty cells from connectedEvents
     non_empty_indices = ~cellfun('isempty', connectedEvents);
     connectedEvents = connectedEvents(non_empty_indices);


     % Check the indices and remove if there are none
     for i = 1:length(connectedEvents)
         check = connectedEvents{i,1}.indices;

         if isempty(check)
             connectedEvents{i,1} = [];
         end
     end

     % Remove empty cells from connectedEvents
     non_empty_indices = ~cellfun('isempty', connectedEvents);
     connectedEvents = connectedEvents(non_empty_indices);

     kComplexes = connectedEvents;


     dataTableKcomplex{2,channel} = kComplexes;
end

    % Save
    save(fullfile(save_path,[filename '_' 'kcomplexes']), 'dataTableKcomplex');
    fprintf('K-complex analysis completed.\n');
end









function filtered_consecutive_indices = find_consecutive_kcomplexes(indices)
    % This function identifies and filters consecutive runs of indices

    % Initialize variables
    consecutive_indices = {};
    current_run = [];

    % Sampling frequency (Fs) in Hz
    Fs = 512;

    % Initialize the minimum and maximum length of consecutive runs in data points
    min_length = 0.49 * Fs;
    max_length = 2.5 * Fs;

    % Iterate through the input array of indices
    for i = 1:length(indices)-1
        if indices(i+1) == indices(i) + 1
            % Indices are consecutive
            current_run = [current_run; indices(i)];  % Store as column vector
        else
            % End of consecutive run
            current_run = [current_run; indices(i)];  % Store as column vector
            % Check if the consecutive run length is within the specified range
            if length(current_run) >= min_length && length(current_run) <= max_length
                % Add valid consecutive run to the cell array
                consecutive_indices{end+1} = current_run;
            end
            current_run = [];  % Reset the current run
        end
    end

    % Add the last index if it's part of a consecutive run
    current_run = [current_run; indices(end)];  % Store as column vector
    % Check the length of the last consecutive run
    if length(current_run) >= min_length && length(current_run) <= max_length
        % Add valid consecutive run to the cell array
        consecutive_indices{end+1} = current_run;
    end

    % Filter consecutive indices based on length
    filtered_consecutive_indices = consecutive_indices';
end










%%-------------------------------------------------------------------------
% SPINDLE DETECTION
%%-------------------------------------------------------------------------





function spindles = Detect_Spindles(Fs, mat_signal_filtered, filename, save_path, Thresholds)

% Default thresholds
thrs_sigma = Thresholds.Power_spindle;       % Sigma band power
thrs_relative = Thresholds.RelativePower_spindle;    % Relative power (sigma band - delta band)
thrs_hilbert = Thresholds.Hilbert_spindle;       % Sigma band
thrs_cwt = Thresholds.CWT_spindle;          % Sigma band
thrs_rms_sigma = Thresholds.RMS_spindle;     % Sigma band
thrs_P2P = Thresholds.P2P_spindle;          % Sigma band

% Alternative annotation name
aasmEvent_name = 'UNSURE';

channels_to_analyse = mat_signal_filtered.Properties.VariableNames;

% Create a table where to save all detected spindles from different
% channels

% Define column names
columnNames = channels_to_analyse;

% Create an empty cell array with specified column names
dataTableSpindle = cell(1, numel(columnNames)); % Start with one row

% Assign column names to the first row of the cell array
dataTableSpindle(1, :) = columnNames;

event = 'spindle';

% Clip the signals to be in range -100,100
clipped_signals = clip(mat_signal_filtered, -100, 100);


for channel = 1 : length(channels_to_analyse)

    chan = cell2mat(channels_to_analyse(channel));
    EEGbf = clipped_signals.(chan);

    fprintf('\n');
    fprintf('Detecting spindles for channel %s', chan);
    fprintf('\n');

    % Window size and window step for caluculating powers in signal
    WinL = round(0.3 * Fs);  % Length of a window (0.3 seconds)
    WinS = floor(0.1 * Fs);  % Step size in samples (0.1 seconds)

    % Filter the data for sigma band and broad band excluding the delta band
    % using the butter function

    f_low = 11; % Lower cutoff frequency
    f_high = 16; % Upper cutoff frequency
    order = 2; % Filter order

    % Design the bandpass filter using the Butterworth filter for Sigma band
    [b, a] = butter(order, [f_low, f_high]/(Fs/2), 'bandpass');

    % Apply the bandpass filter to the signal
    EEG_sigma = filter(b, a, EEGbf);

    % Design the bandpass filter using Butterworth filter for Delta band
    f_low_d = 4.5; % Lower cutoff frequency
    f_high_d = 30; % Upper cutoff frequency

    [b_d, a_d] = butter(order, [f_low_d, f_high_d]/(Fs/2), 'bandpass');

    % Apply the bandpass filter to the signal
    EEG_delta = filter(b_d, a_d, EEGbf);

    % Calculate the absolut powers of sigma and delta filtered signals and
    % calculate the ratio for relative sigma power
    absolute_sigma_power = zeros(size(EEG_sigma));

    padL = floor((0.1/2)*Fs); % Zero padding length for the window

    fprintf('Filtering signal for sigma band. This takes few seconds ');
    fprintf('\n');

    for i = 1:WinS:length(EEG_sigma) - WinL + 1
        window = EEG_sigma(i:i+WinL-1);
        window = padarray(window,padL,'pre');
        absolute_sigma_power(i:i+WinL-1) = log10(sum(window.^2) / WinL);
    end

    absolute_delta_power = zeros(size(EEG_delta));

    fprintf('Excluding delta band from signal. This takes few seconds ');
    fprintf('\n');

    for i = 1:WinS:length(EEG_delta) - WinL + 1
        window = EEG_delta(i:i+WinL-1);
        if i > 1
            window = padarray(window,padL,'pre');
        end
        absolute_delta_power(i:i+WinL-1) = log10(sum(window.^2) / WinL);
    end

    ratio_power = absolute_sigma_power./absolute_delta_power;



    % Window size and window step for caluculating amplitude changes
    taoL = floor(30 * Fs);  % Window (tao) length in samples (60 seconds)
    tao_overlap = round(0.5*taoL);
    tao0L = floor(2 * Fs);  % Window (tao_0) length in samples (2 seconds)
    tao0_overlap = round(0.5*tao0L);

    total_windows = floor((length(EEG_sigma) - taoL) / (taoL - tao_overlap)) + 1;
    Peak2Peak_sigma = zeros(size(EEG_sigma));


    for i = 1:total_windows
        % Calculate start and end indices for the current window
        startIndex = 1 + (i - 1) * (taoL - tao_overlap);
        endIndex = startIndex + taoL - 1;
        % Ensure endIndex does not exceed the length of the signal
        endIndex = min(endIndex, length(EEG_delta));

        current_window = EEG_sigma(startIndex:endIndex);

        tao0_windows = floor((length(current_window) - tao0L) / (tao0L - tao0_overlap)) + 1;

        for j = 1:tao0_windows
            startI = 1 + (j -1) * (tao0L - tao0_overlap);
            endI = startI + tao0L - 1;
            endI = min(endI, length(current_window));
            tao0_window = current_window(startI:endI);
            peak2peak_value = max(tao0_window) - min(tao0_window);
            Peak2Peak_index = startIndex + startI - 1 : startIndex + endI - 1;
            Peak2Peak_sigma(Peak2Peak_index) = peak2peak_value;
        end
    end


    % Calculate the first IMF 
    [imf_s,~,~] = emd(EEG_sigma);
    first_imf = imf_s(:,1);


    % Detect spindle candidates with point system using thresholds. If all
    % thresholds are exceeded spindle is most likely detected. Detected
    % spindles are at least 0.5 second in length.


    % Update the WinL and make the window overlap to be 50%
    WinL = 30*Fs;
    WinL_rms = round(Fs * 0.3);
    overlap = round(0.5*WinL);

    % Calculate the total number of windows
    total_windows = floor((length(EEGbf) - WinL) / (WinL - overlap)) + 1;

    spindles = cell(total_windows, 1);

    h = waitbar(0, 'Detecting spindles over night...', 'Name', 'Spindle detection');

    for i = 1:total_windows

        % Calculate start and end indices for the current window
        startIndex = 1 + (i - 1) * (WinL - overlap);
        endIndex = startIndex + WinL - 1;

        % Ensure endIndex does not exceed the length of the signal
        endIndex = min(endIndex, length(EEGbf));

        % Calculate the rms values for sigma band and original signal
        rms_values_sigma = zeros(1,length(startIndex:endIndex));
        rms_values_bf = zeros(1,length(startIndex:endIndex));

        for k = 1:WinS:length(startIndex:endIndex)-WinL_rms+1
            window_start = startIndex + k - 1;
            window_end = min(window_start + WinL_rms - 1, endIndex);
            window_segment_sigma = EEG_sigma(window_start:window_end);
            window_segment_bf = EEGbf(window_start:window_end);
            window_segment_sigma = padarray(window_segment_sigma, padL, 0, 'pre');
            window_segment_bf = padarray(window_segment_bf, padL, 0, 'pre');
            window_rms_sigma = sqrt(mean(window_segment_sigma.^2));
            window_rms_bf = sqrt(mean(window_segment_bf.^2));
            rms_values_sigma(k:k+WinL_rms-1) = window_rms_sigma;
            rms_values_bf(k:k+WinL_rms-1) = window_rms_bf;
        end

        % Extract a windows from different signal
        power_window = absolute_sigma_power(startIndex:endIndex);
        ratio_window = ratio_power(startIndex:endIndex);

        % Find indices over the threshold
        indices_over_threshold = find(power_window >= thrs_sigma);

        % Check if there are any consecutive runs within the desired length range
        if isempty(indices_over_threshold)
            % If no consecutive indices, break out of the loop and move to the next iteration
            continue;
        end


        % Find consecutive indices
        consecutive_indices = find_consecutive_spindles(indices_over_threshold);
        % Check if there are any consecutive runs within the desired length range
        if isempty(consecutive_indices)
            % If no consecutive indices, break out of the loop and move to the next iteration
            continue;
        end


        % Map indices to correspond the original signals indices
        % for each consecutive run
        filtered_consecutive_indices = cell(size(consecutive_indices));

        for j = 1:length(consecutive_indices)
            % Map indices from power_window to absolute_sigma_power
            absolute_sigma_power_indices = consecutive_indices{j} + startIndex - 1;
            start_time_seconds = absolute_sigma_power_indices(1,1) / Fs;
            % Store filtered indices and starting time as a struct or tuple
            filtered_consecutive_indices{j} = struct('indices', absolute_sigma_power_indices, 'start_sec', start_time_seconds);
        end


        % CALCULATING THE CWT FROM FIRST MODE
        [cwtResult, frequencies] = cwt(first_imf(startIndex:endIndex), 'amor', Fs);

        % Find regions where CWT values exceed the threshold
        highFreqIndices = frequencies >= 10 & frequencies <= 17;
        cwtHighFreq = abs(cwtResult(highFreqIndices, :));
        % Create masks for threshold exceeding regions
        mask = cwtHighFreq > thrs_cwt;


        % Initialize points for the current window
        points_per_consecutive_part = zeros(size(consecutive_indices));

        % Iterate over each consecutive run and check the amplitude
        % changes
        for j = 1:size(consecutive_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = consecutive_indices{j}(1);
            end_index = consecutive_indices{j}(end);
            % Extract the corresponding part amplitude cahnge values
            Peak2Peak_part = Peak2Peak_sigma(start_index:end_index);
            % Process the result based on peak to peak values
            if any(Peak2Peak_part >= thrs_P2P)
                % Increment the counter for the current consecutive part
                points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
            end
            % Process the result based on peak to peak values. subtract
            % one point if the peak to peak value is too high. This
            % indicates the breaking of a spindle.
            if any(Peak2Peak_part >= 100)
                % Increment the counter for the current consecutive part
                points_per_consecutive_part(j) = points_per_consecutive_part(j) - 1;
            end
        end



        % Iterate over each consecutive run and check values for rms.
        % This is done to the original signal to filter out the
        % artefacts
        for j = 1:size(consecutive_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = consecutive_indices{j}(1);
            end_index = consecutive_indices{j}(end);
            % Subset rms_values array between start_index and end_index
            subset_rms_values = rms_values_bf(start_index:end_index);
            % Find indices where the RMS values exceed the threshold
            indices_over_threshold_rms = find(subset_rms_values >= 35, 1);
            % Process the result based on the rms values
            if ~isempty(indices_over_threshold_rms)
                % Increment the counter for the current consecutive part
                points_per_consecutive_part(j) = points_per_consecutive_part(j) - 1;
            end
        end




        % Iterate over each consecutive run and check values in ratio_window
        for j = 1:size(consecutive_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = consecutive_indices{j}(1);
            end_index = consecutive_indices{j}(end);
            % Extract the corresponding part from ratio_window
            ratio_consecutive_part = abs(ratio_window(start_index:end_index,1));
            % Process the result based on the values in ratio_window
            if any(ratio_consecutive_part >= thrs_relative)
                % Increment the counter for the current consecutive part
                points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
            end
        end



        % Iterate over each consecutive run and check rms values in
        % sigma band
        for j = 1:size(consecutive_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = consecutive_indices{j}(1);
            end_index = consecutive_indices{j}(end);
            % Extract the corresponding part from ratio_window
            max_rms_sigma = max(rms_values_sigma(start_index:end_index));
            % Process the result based on the values in ratio_window
            if any(max_rms_sigma >= thrs_rms_sigma)
                % Increment the counter for the current consecutive part
                points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
            end
        end




        % Iterate over each consecutive run and check values in Hilbert
        % trasform of sigma filtered signal
        for j = 1:size(consecutive_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = consecutive_indices{j}(1);
            end_index = consecutive_indices{j}(end);
            % Extract the corresponding part from EEG_sigma
            sigma_consecutive_part = EEG_sigma(start_index:end_index);
            % Compute the envelope of Hilbert transform for sigma_consecutive
            % part
            envelope = abs(hilbert(sigma_consecutive_part));
            % Process the result based on the values in envelope
            if any(envelope >= thrs_hilbert)
                % Process the result based on the values in envelope
                points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
            end
        end




        % Iterate over each consecutive run and check if there are
        % areas that exceeds the CWT threshold.
        for j = 1:size(consecutive_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = consecutive_indices{j}(1);
            end_index = consecutive_indices{j}(end);
            timeIndices = start_index:end_index;
            % Pick up the values from cwt results for candidate area
            candidate_areas = bwconncomp(mask(:,timeIndices));
            % Initialize a variable for the number of areas in the specified time range
            numAreasInRange = candidate_areas.NumObjects;
            % Check if values in candidate area exceeds the cwt
            % threshold
            if numAreasInRange > 0
                % Process the result based on the average value
                points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
            end
        end





        % Check if any consecutive part has accumulated five points
        if any(points_per_consecutive_part == 5)
            % Find indices of consecutive parts with five points
            indices_three_points = find(points_per_consecutive_part == 5);
            % Create a cell array to store spindle indices for the current window
            current_spindle_indices = cell(length(indices_three_points), 1);
            % Iterate over indices and store corresponding segments in current_spindle_indices
            for k = 1:length(indices_three_points)
                current_consecutive_part_indices = filtered_consecutive_indices{indices_three_points(k)};
                current_spindle_indices{k} = current_consecutive_part_indices.indices(1,1):current_consecutive_part_indices.indices(end,1);
                current_spindle_starting_time = filtered_consecutive_indices{indices_three_points(k)};
                current_spindle_starting_time = current_spindle_starting_time.start_sec;
                current_spindle_duration = length(current_spindle_indices{k,1}) / Fs;
                current_spindle_amplitude = max(EEGbf(current_spindle_indices{k,1})) - min(EEGbf(current_spindle_indices{k,1}));
                current_spindle_indices{k} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', columnNames{channel}, 'sampling_frequency', Fs, ...
                                                        'indices', current_spindle_indices{k}, 'start_sec', current_spindle_starting_time, ...
                                                        'duration', current_spindle_duration, 'amplitude', current_spindle_amplitude);
            end
            % Store spindle indices for the current window in all_spindle_indices
            spindles{i} = current_spindle_indices;
        end

        % Update waitbar
        waitbar(i/total_windows, h);

    end
    close(h);

    % Remove empty cells from all_spindle_indices
    non_empty_indices = ~cellfun('isempty', spindles);
    spindles = spindles(non_empty_indices);


    % This part connects those spindle candidates that share same
    % indices due to the overlap of analyzing window
    for i = 1:size(spindles, 1) - 1
        current_indices = [];
        next_indices = [];

        for j = 1:numel(spindles{i})
            current_cell = spindles{i, 1}{j, 1}.indices;
            current_indices{j} = current_cell;
        end

        for k = 1:numel(spindles{i+1})
            next_cell = spindles{i+1, 1}{k, 1}.indices;
            next_indices{k} = next_cell;
        end

        for n = 1:numel(current_indices)
            % Check if the current cell is not empty
            if ~isempty(current_indices)
                    current_set = current_indices{n}; % Extract the set from the single set cell
                    % Check if the next cell has sets
                    if ~isempty(next_indices) %&& size(next_indices, 1) >= 1 % Check if the next cell has at least one set
                        % Find the indices of identical sets in the next cell
                        idx_identical = find(cellfun(@(x) isequal(current_set, x), next_indices));
                        if ~isempty(idx_identical)
                            % Remove all identical sets from the next cell
                            spindles{i + 1}(idx_identical) = [];
                            next_indices(idx_identical) = [];
                        end
                    end
            else
                % Continue the loop if the current cell is empty
                continue;
            end
        end
    end

    % Remove empty cells from all_spindle_indices
    non_empty_indices = ~cellfun('isempty', spindles);
    spindles = spindles(non_empty_indices);

    % Initialize a new cell array to store the resulting structs
    resulting_structs = {};

    % Loop through each cell of spindles
    for i = 1:numel(spindles)
        % Check if the current cell contains structs
        if iscell(spindles{i})
            % If there are structs, extract them individually
            for j = 1:numel(spindles{i})
                resulting_structs{end+1} = spindles{i}{j};
            end
        else
            % If there's only one struct, extract it
            resulting_structs{end+1} = spindles{i};
        end
    end
    spindles = resulting_structs';



    % Connect those spindle candidates that are close to each other
    difference = 0.5 * Fs; % in seconds
    connectedEvents = cell(length(spindles),1);

    i = 1;
    while i <= length(spindles) - 1
        current_event_indices = spindles{i,1}.indices;
        current_event_last = spindles{i,1}.indices(end);
        current_event_start = spindles{i, 1}.start_sec;
        current_event_dur = spindles{i, 1}.duration;
        current_event_amp = spindles{i, 1}.amplitude;
        next_event_indices = spindles{i+1,1}.indices;
        next_event_first = spindles{i+1,1}.indices(1);

        if (next_event_first - current_event_last) <= difference
            combined_indices = current_event_indices(1):next_event_indices(end);
            combined_duration = length(combined_indices)/Fs;
            new_amplitude = max(EEGbf(combined_indices)) - min(EEGbf(combined_indices));
            connectedEvents{i} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', columnNames{channel}, 'sampling_frequency', Fs, ...
                                        'indices', combined_indices, 'start_sec', current_event_start, ...
                                        'duration', combined_duration, 'amplitude', new_amplitude);
            i = i + 2;
        else
            connectedEvents{i} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', columnNames{channel}, 'sampling_frequency', Fs, ...
                                        'indices', current_event_indices, 'start_sec', current_event_start, ...
                                        'duration', current_event_dur, 'amplitude', current_event_amp);
            i = i +1;
        end
    end

    % Remove empty cells from connectedEvents
    non_empty_indices = ~cellfun('isempty', connectedEvents);
    connectedEvents = connectedEvents(non_empty_indices);


    % Check the indices and remove if there are none
    for i = 1:length(connectedEvents)
        check = connectedEvents{i,1}.indices;
        if isempty(check)
            connectedEvents{i,1} = [];
        end
    end

    % Remove empty cells from connectedEvents
    non_empty_indices = ~cellfun('isempty', connectedEvents);
    connectedEvents = connectedEvents(non_empty_indices);
    spindles = connectedEvents;

    dataTableSpindle{2,channel} = spindles;
end
    % Save
    save(fullfile(save_path,[filename '_' 'spindles']), 'dataTableSpindle');
    fprintf('Spindle analysis completed.\n');

end




function filtered_consecutive_indices = find_consecutive_spindles(indices)
    % This function identifies and filters consecutive runs of indices

    % Initialize variables
    consecutive_indices = {};
    current_run = [];

    % Sampling frequency (Fs) in Hz
    Fs = 512;

    % Initialize the minimum and maximum length of consecutive runs in data points
    min_length = 0.5 * Fs;
    max_length = 10 * Fs;

    % Iterate through the input array of indices
    for i = 1:length(indices)-1
        if indices(i+1) == indices(i) + 1
            % Indices are consecutive
            current_run = [current_run; indices(i)];  % Store as column vector
        else
            % End of consecutive run
            current_run = [current_run; indices(i)];  % Store as column vector
            % Check if the consecutive run length is within the specified range
            if length(current_run) >= min_length && length(current_run) <= max_length
                % Add valid consecutive run to the cell array
                consecutive_indices{end+1} = current_run;
            end
            current_run = [];  % Reset the current run
        end
    end

    % Add the last index if it's part of a consecutive run
    current_run = [current_run; indices(end)];  % Store as column vector
    % Check the length of the last consecutive run
    if length(current_run) >= min_length && length(current_run) <= max_length
        % Add valid consecutive run to the cell array
        consecutive_indices{end+1} = current_run;
    end

    % Filter consecutive indices based on length
    filtered_consecutive_indices = consecutive_indices';
end


function clipped_array = clip(array, lower_bound, upper_bound)
    % This function clips the signals from all the channels
    clipped_array = min(max(array, lower_bound), upper_bound);
end












%%-------------------------------------------------------------------------
% DELTA WAVE / SWS DETECTION
%%-------------------------------------------------------------------------





function delta_waves = Detect_delta_waves(Fs, mat_signal_filtered, filename, save_path, rawEEGSignals, Thresholds)



channels_to_analyse = mat_signal_filtered.Properties.VariableNames;

% Default thresholds
thrs_power = Thresholds.Power_SWS;           % Delta band 
thrs_hilbert = Thresholds.Hilbert_SWS;         % Delta band 
thrs_rms = Thresholds.RMS_SWS;              % Delta band
thrs_wave_length = Thresholds.WaveLength_SWS;    % Delta band 
thrs_cwt = Thresholds.CWT_SWS;              % Delta band 

% Alternative annotation name
aasmEvent_name = 'UNSURE';


% Create a table where to save all detected slow waves from different
% channels

% Define column names
columnNames = channels_to_analyse;
% Create an empty cell array with specified column names
dataTableSWS = cell(1, numel(columnNames)); % Start with one row
% Assign column names to the first row of the cell array
dataTableSWS(1, :) = columnNames;

event = 'Delta_waves';


for channel = 1 : length(channels_to_analyse)
    chan = cell2mat(channels_to_analyse(channel));
    EEGbf = mat_signal_filtered.(chan);

    rawSignal = rawEEGSignals.(chan).signal;

    fprintf('\n');
    fprintf('Detecting delta waves for channel %s', chan);
    fprintf('\n');

    % Window size and window step for caluculating powers in signal
    WinL = round(2 * Fs);  % Window length in samples
    WinS = floor(1 * Fs);  % Step size in samples

    % Filter the data for delta band
    f_low = 0.1; % Lower cutoff frequency
    f_high = 4.5; % Upper cutoff frequency
    order = 2; % Filter order

    % Design the bandpass filter using the Butterworth filter for Delta band
    [b, a] = butter(order, [f_low, f_high]/(Fs/2), 'bandpass');

    % Apply the bandpass filter to the signal
    EEG_delta = filter(b, a, EEGbf);
    absolute_delta_power = zeros(size(EEG_delta));
    padL = floor((0.1/2)*Fs); % Zero padding length for the window

    fprintf('Calculating delta power. This takes few seconds ');
    fprintf('\n');

    for i = 1:WinS:length(EEG_delta) - WinL + 1
        window = EEG_delta(i:i+WinL-1);
        if i > 1
            window = padarray(window,padL,'pre');
        end
        absolute_delta_power(i:i+WinL-1) = log10(sum(window.^2) / WinL);
    end

    % Calculate the imf from delta filtered signal using EMD function
    [imf_d,~,~] = emd(EEG_delta);
    % Choose the second imf to use
    second_imf = imf_d(:,2);

    % Update the window lengths for RMS
    WinL = 1*Fs;
    overlap = round(0.5*WinL);

    % Calculate the total number of windows
    total_windows = floor((length(EEGbf) - WinL) / (WinL - overlap)) + 1;

    % Initialize RMS values array for delta band
    rms_values_delta = zeros(size(EEG_delta));
    rms_raw = zeros(size(EEGbf));

    % Loop through the signal and calculate RMS for each window
    for i = 1:total_windows
        % Calculate the start and end indices of the current window
        start_idx = (i-1) * (WinL - overlap) + 1;
        end_idx = start_idx + WinL - 1;

        % Ensure the indices are within the signal length
        if end_idx > length(EEG_delta)
            end_idx = length(EEG_delta);
        end

        % Extract the windowed segments
        window_segment = EEG_delta(start_idx:end_idx);
        window_segment_raw = rawSignal(start_idx:end_idx);

        % Calculate the RMS value for the current window
        rms_value = sqrt(mean(window_segment .^ 2));
        rms_value_raw = sqrt(mean(window_segment_raw .^ 2));

        % Assign the RMS value to the corresponding positions
        rms_values_delta(start_idx:end_idx) = rms_value;
        rms_raw(start_idx:end_idx) = rms_value_raw;
    end


    % Update the WinL and make the window overlap to be 50%
    WinL = 30*Fs;
    overlap = round(0.5*WinL);

    % Calculate the total number of windows
    total_windows = floor((length(EEGbf) - WinL) / (WinL - overlap)) + 1;

    % Initialize the score table for possible candidates for later
    % evaluation.
    delta_waves = cell(total_windows, 1);

    h = waitbar(0, 'Detecting SWS over night...', 'Name', 'SWS detection');


    for i = 1:total_windows
        % Calculate start and end indices for the current window
        startIndex = 1 + (i - 1) * (WinL - overlap);
        endIndex = startIndex + WinL - 1;
        % Ensure endIndex does not exceed the length of the signal
        endIndex = min(endIndex, length(EEGbf));

        % Extract the corresponding part from EEG_delta
        delta_consecutive_part = EEG_delta(startIndex:endIndex);
        % Compute the envelope of Hilbert transform for delta_consecutive
        % part
        envelope = abs(hilbert(delta_consecutive_part));

        % CALCULATING THE CWT FROM SECOND MODE
        [cwtResult, frequencies] = cwt(second_imf(startIndex:endIndex), 'amor', Fs);

        % Find regions where CWT values exceed the threshold
        highFreqIndices = frequencies >= 0.1 & frequencies <= 30;
        cwtHighFreq = abs(cwtResult(highFreqIndices, :));
        % Create masks for threshold exceeding regions
        mask = cwtHighFreq > thrs_cwt;
        % Find connected components in the masks
        CC = bwconncomp(mask);

        % Initialize a cell array to store the column indices for each connected component
        componentIndices = cell(CC.NumObjects, 1);

        % Minimum length for each component
        minLength = 2 * Fs;

        % Iterate over each connected component
        for j = 1:CC.NumObjects
            % Get the linear indices of the current component
            linearIndices = CC.PixelIdxList{j};
            % Convert linear indices to row and column subscripts
            [~, cols] = ind2sub(size(mask), linearIndices);
            % Store the unique column indices
            componentIndices{j} = unique(cols);
        end

        % Merge overlapping components
        componentIndices = mergeOverlappingComponents(componentIndices);

        % Remove components that are shorter than the specified minimum length
        validComponents = cellfun(@(x) length(x) >= minLength, componentIndices);
        componentIndices = componentIndices(validComponents)';

        % Map indices from power_window to EEGbf for each component
        filtered_consecutive_indices = cell(size(componentIndices));

        for j = 1:length(componentIndices)
            real_indices = componentIndices{j} + startIndex - 1;
            filtered_consecutive_indices{j} = real_indices;
        end


        % Combine components that are within 0.5 * Fs distance from each other
        maxDistance = 0.5 * Fs;
        filtered_consecutive_indices = combineCloseComponents(filtered_consecutive_indices, maxDistance)';


        for j = 1:length(filtered_consecutive_indices)

            start_time_seconds = filtered_consecutive_indices{j,1}(1) / Fs;
            % Store filtered indices and starting time as a struct or tuple
            filtered_consecutive_indices{j} = struct('indices', filtered_consecutive_indices{j,1}, 'start_sec', start_time_seconds);
        end


        if isempty(filtered_consecutive_indices)
            continue;
        end

        event_indices = cell(size(filtered_consecutive_indices));

        % Pick up the indices to speed up the process
        for j = 1:length(filtered_consecutive_indices)
            event_indices{j} = filtered_consecutive_indices{j,1}.indices;
        end

        % Initialize points for the current window
        points_per_consecutive_part = zeros(size(event_indices));


        % Check if the power exceeds the threshold in candidate area
        for j = 1:size(event_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = event_indices{j}(1);
            end_index = event_indices{j}(end);
            % Extract the corresponding part from EEG_sigma
            power_part = absolute_delta_power(start_index:end_index,1);
            % Process the result based on the wave length
            if any(power_part >= thrs_power)
               points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
            end
        end



        % Wave length check. In this part the code checks the wave
        % length of the signal in detected area.
        for j = 1:size(filtered_consecutive_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = event_indices{j}(1);
            end_index = event_indices{j}(end);
            % Extract the corresponding part from EEG_sigma
            delta_consecutive_part = EEG_delta(start_index:end_index,1);
            [~, locs] = findpeaks(delta_consecutive_part);
            if length(locs) > 1
                period_samples = diff(locs); % Differences between consecutive peaks
                periods = period_samples / Fs; % Convert to time periods
                % Count the number of periods greater than the threshold
                num_long_periods = sum(periods >= thrs_wave_length);
                % Process the result based on the count of long periods
                if num_long_periods >= 3
                    points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
                end
            end
        end




        for j = 1:size(event_indices, 1)
            start_index = event_indices{j}(1);
            end_index = event_indices{j}(end);
            % Extract the corresponding part from EEG_sigma
            rms_part_delta = rms_values_delta(start_index:end_index);
            rms_part_raw = rms_raw(start_index:end_index);
            if max(rms_part_delta) >= thrs_rms
                % Process the result based on the values in envelope
                points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
            end
            if max(rms_part_raw) >= 50
                % Process the result based on the values in envelope
                points_per_consecutive_part(j) = points_per_consecutive_part(j) - 1;
            end
        end


        for j = 1:size(event_indices, 1)
            start_index = event_indices{j}(1);
            end_index = event_indices{j}(end);
            % Extract the corresponding part from EEG
            rawSignalPart = rawSignal(start_index:end_index);
            FT_EEG = fft(rawSignalPart);
            L = length(rawSignalPart);

            f = Fs/L*(0:L-1);
            % Filter out frequencies
            filtered_FT_EEG = FT_EEG(f >= 60 & f <= 90);
            fft_magnitude = abs(filtered_FT_EEG);
            mean_fft_magnitude = mean(fft_magnitude);

            if mean_fft_magnitude >= 100
                points_per_consecutive_part(j) = points_per_consecutive_part(j) - 1;
            end 

        end



        for j = 1:size(event_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = event_indices{j}(1) - startIndex + 1;
            end_index = event_indices{j}(end) - startIndex + 1;
            envelope_area = envelope(start_index:end_index);
            % Process the result based on the values in envelope
            if any(envelope_area >= thrs_hilbert)
                % Process the result based on the values in envelope
                points_per_consecutive_part(j) = points_per_consecutive_part(j) + 1;
            end
        end



        % Check if any consecutive part has accumulated four points
        if any(points_per_consecutive_part == 4)
            % Find indices of consecutive parts with four points
            indices_four_points = find(points_per_consecutive_part == 4);
            % Create a cell array to store sws indices for the current window
            current_sws_indices = cell(length(indices_four_points), 1);
            % Iterate over indices and store corresponding segments in current_sws_indices
            for k = 1:length(indices_four_points)
                current_consecutive_part_indices = filtered_consecutive_indices{indices_four_points(k)};
                current_sws_indices{k} = current_consecutive_part_indices.indices(1,1):current_consecutive_part_indices.indices(end,1);
                current_sws_starting_time = filtered_consecutive_indices{indices_four_points(k)};
                current_sws_starting_time = current_sws_starting_time.start_sec;
                current_sws_duration = length(current_sws_indices{k,1}) / Fs;
                current_sws_amplitude = max(EEGbf(current_sws_indices{k,1})) - min(EEGbf(current_sws_indices{k,1}));
                current_sws_EEGbf = EEGbf(current_sws_indices{k,1});
                current_sws_indices{k} = struct('indices', current_sws_indices{k}, 'start_sec', current_sws_starting_time, ...
                                                        'duration', current_sws_duration, 'amplitude', current_sws_amplitude, ...
                                                        'filtered_signal', current_sws_EEGbf);
            end
            % Store sws indices for the current window in all_sws_indices
            delta_waves{i} = current_sws_indices;
        end

        % Update waitbar
        waitbar(i/total_windows, h);
    end

    close(h);

    % Remove empty cells from all_sws_indices
    non_empty_indices = ~cellfun('isempty', delta_waves);
    delta_waves = delta_waves(non_empty_indices);


    for i = 1:size(delta_waves, 1) - 1
        current_indices = [];
        next_indices = [];
        for j = 1:numel(delta_waves{i})
            current_cell = delta_waves{i, 1}{j, 1}.indices;
            current_indices{j} = current_cell;
        end

        for k = 1:numel(delta_waves{i+1})
            next_cell = delta_waves{i+1, 1}{k, 1}.indices;
            next_indices{k} = next_cell;
        end

        for n = 1:numel(current_indices)
            % Check if the current cell is not empty
            if ~isempty(current_indices)
                    current_set = current_indices{n}; % Extract the set from the single set cell
                    % Check if the next cell has sets
                    if ~isempty(next_indices) %&& size(next_indices, 1) >= 1 % Check if the next cell has at least one set
                        % Find the indices of identical sets in the next cell
                        idx_identical = find(cellfun(@(x) isequal(current_set, x), next_indices));
                        if ~isempty(idx_identical)
                            % Remove all identical sets from the next cell
                            delta_waves{i + 1}(idx_identical) = [];
                            next_indices(idx_identical) = [];
                        end
                    end
            else
                % Continue the loop if the current cell is empty
                continue;
            end
        end
    end

    % Remove empty cells from all_sws_indices
    non_empty_indices = ~cellfun('isempty', delta_waves);
    delta_waves = delta_waves(non_empty_indices);

    % Initialize a new cell array to store the resulting structs
    resulting_structs = {};

    % Loop through each cell of sws
    for i = 1:numel(delta_waves)
        % Check if the current cell contains structs
        if iscell(delta_waves{i})
            % If there are structs, extract them individually
            for j = 1:numel(delta_waves{i})
                resulting_structs{end+1} = delta_waves{i}{j};
            end
        else
            % If there's only one struct, extract it
            resulting_structs{end+1} = delta_waves{i};
        end
    end

    delta_waves = resulting_structs';


    % Connect those slow wave sequences that are close to each other
    difference = 2 * Fs; % in seconds
    connectedEvents = cell(length(delta_waves),1);

    i = 1;
    while i <= length(delta_waves) - 1
        current_event_indices = delta_waves{i,1}.indices;
        current_event_last = delta_waves{i,1}.indices(end);
        current_event_start = delta_waves{i, 1}.start_sec;
        current_event_dur = delta_waves{i, 1}.duration;
        current_event_amp = delta_waves{i, 1}.amplitude;
        current_event_filt = delta_waves{i, 1}.filtered_signal;
        next_event_indices = delta_waves{i+1,1}.indices;
        next_event_first = delta_waves{i+1,1}.indices(1);

        if (next_event_first - current_event_last) <= difference
            combined_indices = current_event_indices(1):next_event_indices(end);
            combined_duration = length(combined_indices)/Fs;
            new_amplitude = max(EEGbf(combined_indices)) - min(EEGbf(combined_indices));
            connectedEvents{i} = struct('indices', combined_indices, 'start_sec', current_event_start, ...
                                        'duration', combined_duration, 'amplitude', new_amplitude, ...
                                        'filtered_signal', EEGbf(combined_indices));
            i = i + 2;
        else
            connectedEvents{i} = struct('indices', current_event_indices, 'start_sec', current_event_start, ...
                                        'duration', current_event_dur, 'amplitude', current_event_amp, ...
                                        'filtered_signal', current_event_filt);
            i = i +1;
         end
     end

    % Remove empty cells from connectedEvents
    non_empty_indices = ~cellfun('isempty', connectedEvents);
    connectedEvents = connectedEvents(non_empty_indices);

    delta_waves = connectedEvents;


    dataTableSWS{2,channel} = delta_waves;




    % Remove structures that are empty due merging
    events = dataTableSWS{2,channel};
    for i = 1 : numel(dataTableSWS{2,channel})
        current_SWS_indices = events{i,1}.indices;
        if isempty(current_SWS_indices)
            dataTableSWS{2,channel}{i,1} = [];
        end
    end

    % Remove empty cells from dataTableSWS
    non_empty_indices = ~cellfun('isempty', dataTableSWS{2,channel});
    dataTableSWS{2,channel} = dataTableSWS{2,channel}(non_empty_indices);



    % Make a check for connected events to remove duplicates or events that
    % have few indices difference
    difference = 1 * Fs; % in seconds
    events = dataTableSWS{2,channel};
    connectedEvents = cell(length(events),1);

    i = 1;
    while i <= length(events) - 1
        current_event_indices = events{i,1}.indices;
        current_event_last = events{i,1}.indices(end);
        current_event_start = events{i, 1}.start_sec;
        current_event_dur = events{i, 1}.duration;
        current_event_amp = events{i, 1}.amplitude;
        next_event_indices = events{i+1,1}.indices;
        next_event_first = events{i+1,1}.indices(1);

        if (next_event_first - current_event_last) <= difference
            combined_indices = current_event_indices(1):next_event_indices(end);
            combined_duration = length(combined_indices)/Fs;
            new_amplitude = max(EEGbf(combined_indices)) - min(EEGbf(combined_indices));
            connectedEvents{i} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', columnNames{channel}, 'sampling_frequency', Fs, ...
                                        'indices', combined_indices, 'start_sec', current_event_start, ...
                                        'duration', combined_duration, 'amplitude', new_amplitude);
            i = i + 2;
        else
           connectedEvents{i} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', columnNames{channel}, 'sampling_frequency', Fs, ...
                                       'indices', current_event_indices, 'start_sec', current_event_start, ...
                                       'duration', current_event_dur, 'amplitude', current_event_amp);
            i = i +1;
        end
    end

    % Remove empty cells from connectedEvents
    non_empty_indices = ~cellfun('isempty', connectedEvents);
    connectedEvents = connectedEvents(non_empty_indices);

    events = connectedEvents;
    dataTableSWS{2,channel} = events;
end

    delta_waves = dataTableSWS;


    % Save
    save(fullfile(save_path,[filename '_' 'delta_waves']), 'dataTableSWS');
    fprintf('Delta-wave analysis completed.\n');

end







% Function to merge overlapping components
function mergedComponents = mergeOverlappingComponents(components)
    mergedComponents = {};
    n = length(components);
    merged = false(1, n);

    for i = 1:n
        if merged(i)
            continue;
        end

        % Start with the current component
        currentComponent = components{i};
        startIdx = min(currentComponent);
        endIdx = max(currentComponent);

        % Check for overlaps with subsequent components
        for j = i+1:n
            if merged(j)
                continue;
            end

            % Check if there is overlap
            if ~isempty(intersect(currentComponent, components{j}))
                % Merge the components
                currentComponent = union(currentComponent, components{j});
                startIdx = min(startIdx, min(components{j}));
                endIdx = max(endIdx, max(components{j}));
                merged(j) = true;
            end
        end

        % Store the merged component
        mergedComponents{end+1} = (startIdx:endIdx)';
        merged(i) = true;
    end

     for i = 1:n
        if ~merged(i)
            mergedComponents{end+1} = components{i};
        end
    end
end



% Function to combine close components
function combinedComponents = combineCloseComponents(components, maxDistance)
    % Sort components by their starting index
    startIdxs = cellfun(@(x) min(x), components);
    [~, sortedIdx] = sort(startIdxs);
    components = components(sortedIdx);

    n = length(components);
    combinedComponents = {};
    merged = false(1, n);

    for i = 1:n
        if merged(i)
            continue;  % Skip already combined components
        end

        % Start with the current component
        currentComponent = components{i};
        startIdx = min(currentComponent);
        endIdx = max(currentComponent);

        % Combine components that are within maxDistance from each other
        for j = i+1:n
            if merged(j)
                continue;  % Skip already combined components
            end

            % Check if the distance between endIdx and startIdx of the next component is within maxDistance
            if startIdxs(j) - endIdx <= maxDistance
                % Combine the components
                currentComponent = union(currentComponent, components{j});
                startIdx = min(startIdx, min(components{j}));
                endIdx = max(endIdx, max(components{j}));
                merged(j) = true;  % Mark as combined
            end
        end

        % Store the combined component
        combinedComponents{end+1} = (startIdx:endIdx)';
    end
end








%%-------------------------------------------------------------------------
% CAP DETECTION
%%-------------------------------------------------------------------------







function cap = Detect_CAP(Fs, mat_signal_filtered, filename, save_path, Thresholds)

channels_to_analyse = mat_signal_filtered.Properties.VariableNames;

% Default thresholds
thrs_power = Thresholds.Power_CAP;         % Filtered signal (base frequency)
thrs_descriptor = Thresholds.Desc_CAP;  % For all band frequencies
thrs_cwt = Thresholds.CWT_CAP;          % Filtered signal (base frequency)
amp_diff = Thresholds.AmpDiff_CAP;         % Filtered signal (base frequency)

% Alternative annotation name
aasmEvent_name = 'UNSURE';

signals = mat_signal_filtered.Variables;
signalNam = channels_to_analyse;

bands = {'Delta', 'Theta', 'Alpha', 'Sigma', 'Beta'};
band_freq = [0.5, 4; 4, 8; 8, 12; 12, 16; 16, 30]; % Frequency ranges for each band

% Initialize cell array to store filtered signals
filtered_signals = cell(numel(signalNam), numel(bands));

for i = 1 : numel(signalNam)
    signal = signals(:,i);
    % Filter the signals for each frequency band
    for j = 1:numel(bands)
        % Design the bandpass filter using the Butterworth filter
        [b, a] = butter(2, band_freq(j, :)/(Fs/2), 'bandpass');
        % Apply the bandpass filter to the signal
        filtered_signals{i, j} = filter(b, a, signal);
    end
end


EEGbf = signals(:,1);

WinL = 2*Fs;
overlap = round(0.5*WinL);

% Calculate the total number of windows
total_windows = floor((length(EEGbf) - WinL) / (WinL - overlap)) + 1;

% Initialize RMS values array for unfiltered signal
rms_bf = zeros(length(EEGbf),1);

% Loop through the signal and calculate RMS for each window. Point is to
% find large changes in rms valiues to spot artefacts
for i = 1:total_windows
    % Calculate the start and end indices of the current window
    start_idx = (i-1) * (WinL - overlap) + 1;
    end_idx = start_idx + WinL - 1;
    % Ensure the indices are within the signal length
    if end_idx > length(EEGbf)
        end_idx = length(EEGbf);
    end
    % Extract the windowed segment from the delta band signal
    window_segment_bf = EEGbf(start_idx:end_idx,1);
    % Calculate the RMS value for the current window
    rms_value_bf = sqrt(mean(window_segment_bf .^ 2));
    % Assign the RMS value to the corresponding positions in rms_values_delta
    rms_bf(start_idx:end_idx) = rms_value_bf;
end




% Create a table with named columns
Signal_info = table('Size', [1, numel(signalNam) + 1], ...
                    'VariableTypes', [{'cell'}, repmat({'cell'}, 1, numel(signalNam))], ...
                    'VariableNames', [{'Type'}, signalNam]);

% Assign the signal band names to the first row and channel names to the second row
Signal_info{1, 'Type'} = {'Signal power'};


powerWinL = floor(3*Fs);
overlap = round(0.5*powerWinL);

for i = 1:numel(signalNam)
        channel = signalNam{i};
        signal = signals(:,i);

        total_windows = floor((length(signal) - powerWinL) / (powerWinL - overlap)) + 1;

        % Initialize arrays to hold power and contribution count
        power = zeros(size(signal));
        contribution_count = zeros(size(signal));

        h = waitbar(0, sprintf('Calculating power for %s...', channel), 'Name', 'CAP detection');

        for j = 1:total_windows
            % Calculate start and end indices for the current window
            startIndex = 1 + (j - 1) * (powerWinL - overlap);
            endIndex = startIndex + powerWinL - 1;
            % Ensure endIndex does not exceed the length of the signal
            endIndex = min(endIndex, length(signal));
            % Extract the current window
            current_window = signal(startIndex:endIndex);
            % Calculate power as the sum of squared magnitudes of DFT coefficients
            window_power = log10(sum(current_window.^2) / powerWinL);
            % Accumulate power and track contributions
            indices = startIndex:endIndex;
            power(indices) = power(indices) + window_power;
            contribution_count(indices) = contribution_count(indices) + 1;

            waitbar(j / total_windows, h, sprintf('Calculating power for %s: %d%%', channel, round(j / total_windows * 100)));
        end
        close(h);

        % Normalize the power by the contribution count
        power = power ./ max(contribution_count, 1);

        Signal_info.(channel){1} = power;

end







    descriptors = cell(numel(signalNam),numel(bands));


    for n = 1 : numel(signalNam)
        fprintf('Calculating EMD for delta band: %s\n', signalNam{n});
        fprintf('\n');
        [imf_delta,~,~] = emd(filtered_signals{n,1});

        fprintf('Calculating EMD for theta band: %s\n', signalNam{n});
        fprintf('\n');
        [imf_theta,~,~] = emd(filtered_signals{n,2});

        fprintf('Calculating EMD for alpha band: %s\n', signalNam{n});
        fprintf('\n');
        [imf_alpha,~,~] = emd(filtered_signals{n,3});

        fprintf('Calculating EMD for sigma band: %s\n', signalNam{n});
        fprintf('\n');
        [imf_sigma,~,~] = emd(filtered_signals{n,4});

        fprintf('Calculating EMD for beta band: %s\n', signalNam{n});
        fprintf('\n');
        [imf_beta,~,~] = emd(filtered_signals{n,5});

        second_imf_delta = imf_delta(:,2);
        second_imf_theta = imf_theta(:,2);
        second_imf_alpha = imf_alpha(:,2);
        second_imf_sigma = imf_sigma(:,2);
        second_imf_beta = imf_beta(:,2);

        imfs = {second_imf_delta, second_imf_theta, second_imf_alpha, second_imf_sigma, second_imf_beta};


        % Calculate the amplitude changes in each IMF
        mean_amplitudes_short = cell(1,numel(imfs));
        mean_amplitudes_long = cell(1,numel(imfs));
        WinL_short = 2 * Fs;
        WinL_long = 60 * Fs;
        overlap_short = 0.5*WinL_short;
        overlap_long = 0.5*WinL_long;

        for j = 1:numel(imfs)
            imf_signal = imfs{j};
            mean_amplitudes_short{j} = zeros(length(imf_signal),1);
            mean_amplitudes_long{j} = zeros(length(imf_signal),1);
            short_average = zeros(length(imf_signal),1);
            long_average = zeros(length(imf_signal),1);

            % Total windows for longer window average
            total_windows = floor((length(imf_signal) - WinL_long) / (WinL_long - overlap_long)) + 1;

            h = waitbar(0, ['Calculating short window averages for ', bands{j}, ' band IMF.', newline, ...
                'Signal: ', signalNam{n}], 'Name', 'CAP detection');

            for k = 1:total_windows
                startIndex = 1 + (k - 1) * (WinL_long - overlap_long);
                endIndex = startIndex + WinL_long - 1;

                % Ensure endIndex does not exceed the length of the signal
                endIndex = min(endIndex, length(imf_signal));

                segment = imf_signal(startIndex:endIndex);
                mean_segment = mean(abs(segment));
                long_average(startIndex:endIndex) = mean_segment;

                waitbar(k / total_windows, h);
            end
            close(h);

            mean_amplitudes_long{j} = long_average;

            % Total windows for shorter window average
            total_windows = floor((length(imf_signal) - WinL_short) / (WinL_short - overlap_short)) + 1;

            h = waitbar(0, ['Calculating long window averages for ', bands{j}, ' band IMF.', newline, ...
                'Signal: ', signalNam{n}], 'Name', 'CAP detection');

            for k = 1:total_windows
                startIndex = 1 + (k - 1) * (WinL_short - overlap_short);
                endIndex = startIndex + WinL_short - 1;

                % Ensure endIndex does not exceed the length of the signal
                endIndex = min(endIndex, length(imf_signal));

                segment = imf_signal(startIndex:endIndex);
                mean_segment = mean(abs(segment));
                short_average(startIndex:endIndex) = mean_segment;

                waitbar(k / total_windows, h);
            end
            close(h);

            mean_amplitudes_short{j} = short_average;
        end





        % Window length and overlap
        WinL = 2 * Fs; % 2 seconds window length
        overlap = 0.5 * WinL; % 50% overlap

        for j = 1:numel(imfs)
            imf_signal = imfs{j}; % Current IMF signal
            mean_long = mean_amplitudes_long{j}; % Long window average
            mean_short = mean_amplitudes_short{j}; % Short window average

            descriptors{n, j} = zeros(size(imf_signal)); % Initialize descriptor for this IMF

            h = waitbar(0, ['Calculating descriptors for ', bands{j}, ' band', newline, ...
                'Signal: ', signalNam{n}], 'Name', 'CAP detection');

            total_windows = floor((length(imf_signal) - WinL) / (WinL - overlap)) + 1;

            for k = 1:total_windows
                startIndex = 1 + (k - 1) * (WinL - overlap);
                endIndex = startIndex + WinL - 1;

                % Ensure endIndex does not exceed the length of the signal
                endIndex = min(endIndex, length(imf_signal));

                % Short and long average for this window
                S = mean_short(startIndex:endIndex);
                L = mean_long(startIndex:endIndex);

                % Descriptor calculation
                descriptor = (S - L) ./ L;

                % Assign the descriptor value to the entire window
                descriptors{n,j}(startIndex:endIndex) = descriptor;

                waitbar(k / total_windows, h);
            end
            close(h);
        end
    end





    % Parameters
    epochL = 60 * Fs; % Length of epoch in samples
    epoch_overlap = round(0.5 * epochL); % Overlap in samples
    epochs = floor((length(EEGbf) - epochL) / (epochL - epoch_overlap)) + 1;

    merge_distance = 2 * Fs;
    min_region_length = 2 * Fs;

    A_phase_candidates = cell(0,numel(signalNam));

    h = waitbar(0, 'Finding A phase candidates...', 'Name', 'CAP detection');

    for i = 1:epochs
        startIndex = 1 + (i - 1) * (epochL - epoch_overlap);
        endIndex = startIndex + epochL - 1;

        % Ensure endIndex does not exceed the length of the signal
        endIndex = min(endIndex, length(EEGbf));

        all_binary_images = zeros(length(descriptors), epochL);


        % Iterate over each CWT result to plot masks
        for j = 1:length(descriptors)
            descriptor_C3 = descriptors{1,j};
            descriptor_C4 = descriptors{2,j};

            % Find regions where descriptors exceed the threshold
            exceed_C3 = descriptor_C3(startIndex:endIndex) >= thrs_descriptor;
            exceed_C4 = descriptor_C4(startIndex:endIndex) >= thrs_descriptor;

            % Find the intersection of exceedance regions
            combined_exceed = exceed_C3 & exceed_C4;

            % Find start and end indices of contiguous regions
            start_regions = find(diff([0; combined_exceed; 0]) == 1);
            end_regions = find(diff([0; combined_exceed; 0]) == -1);

            % Ensure regions are aligned
            if length(start_regions) > length(end_regions)
                end_regions = [end_regions; length(combined_exceed)];
            elseif length(end_regions) > length(start_regions)
                start_regions = [1; start_regions];
            end

            if isempty(start_regions) || isempty(end_regions)
                merged_binary_image = zeros(size(combined_exceed));
            else
                merged_start = [];
                merged_end = [];
                current_start = start_regions(1);
                current_end = end_regions(1);

                for k = 2:length(start_regions)
                    if start_regions(k) - current_end <= merge_distance
                        current_end = end_regions(k);
                    else
                        merged_start = [merged_start; current_start];
                        merged_end = [merged_end; current_end];
                        current_start = start_regions(k);
                        current_end = end_regions(k);
                    end
                end
                merged_start = [merged_start; current_start];
                merged_end = [merged_end; current_end];

                % Remove regions shorter than min_region_length
                valid_regions = (merged_end - merged_start + 1) >= min_region_length;
                merged_start = merged_start(valid_regions);
                merged_end = merged_end(valid_regions);

                % Clip end indices to stay within the bounds of epochL
                merged_end = min(merged_end, epochL);

                % Create a binary image for remaining valid regions
                merged_binary_image = zeros(1, epochL);
                for k = 1:length(merged_start)
                    merged_binary_image(merged_start(k):merged_end(k)) = 1;
                end

                all_binary_images(j,:) = merged_binary_image;
            end
        end




        % Combine all binary images into a single binary image
        combined_binary_image = any(all_binary_images, 1);

        % Find start and end indices of contiguous regions in the combined binary image
        start_regions = find(diff([0, combined_binary_image, 0]) == 1);
        end_regions = find(diff([0, combined_binary_image, 0]) == -1);

        % Ensure regions are aligned
        if length(start_regions) > length(end_regions)
            end_regions = [end_regions, length(combined_binary_image)];
        elseif length(end_regions) > length(start_regions)
            start_regions = [1, start_regions];
        end




        for k = 1:numel(start_regions)
            region_start = start_regions(k);
            region_end = end_regions(k);

            % Clip end indices to stay within the bounds of epochL
            region_end = min(region_end, epochL);

            % Extract the region from all_binary_images
            region_data = all_binary_images(:, region_start:region_end);
            region_sum = sum(region_data, 1);

            % Check if the sum of columns is at least 3 at least one point
            % and if any row is not entirely zero
            if any(region_sum >= 3)
                % Check if first two rows are empty
                first_two_empty = all(all(region_data(1:2, :) == 0, 2));
                % Check if last three rows are empty
                last_three_empty = all(all(region_data(end-2:end, :) == 0, 2));

                % Determine subtype based on the conditions
                if first_two_empty && ~last_three_empty
                    subtype = 'A1';
                elseif ~first_two_empty && last_three_empty
                    subtype = 'A3';
                else
                    subtype = 'A2';
                end

                % Add valid region and subtype to the lists
                valid_region = (region_start:region_end)'; % Convert to column vector
                % Map indices to real signal indices
                valid_region = valid_region + startIndex - 1;

                A_phase_candidates{end+1, 1} = valid_region;
                A_phase_candidates{end, 2} = subtype; % Add subtype
            end
        end




        % Update the waitbar
        waitbar(i / epochs, h, sprintf('Finding A phase candidates...'), 'Name', 'CAP detection');
    end

    close(h);






    % Descriptors show the local amplitude change that points out from the
    % background EEG. Next thing is to find the indices for each cahnnel where the
    % indices exceeds the thresholds. Later compare these areas between the
    % filtered channels to find the CAP A phase candidates for each epoch
    % that are 60 seconds long. Both channels must exceed the
    % threshold for A phase candidate to be accepted.


    % Connect those A phases that share common indices due to the window
    % overlap in candidate detection

    i = 1;
    while i < length(A_phase_candidates)
        current_A = A_phase_candidates{i,1};
        current_subtype = A_phase_candidates{i,2};
        next_A = A_phase_candidates{i+1,1};

        if any(ismember(current_A, next_A))
            % Combine the unique elements of the two rows
            combined_A = unique([current_A; next_A]);

            % Store the combined result in the current row
            A_phase_candidates{i,1} = combined_A;
            A_phase_candidates{i,2} = current_subtype;

            % Remove the next row, since it has been merged
            A_phase_candidates(i+1,:) = [];
        else
            i = i + 1;
        end
    end




    points_per_candidate = zeros(size(A_phase_candidates(:,1)));

    first_signal_power = Signal_info.(signalNam{1}){1,1};
    second_signal_power = Signal_info.(signalNam{2}){1,1};


    h = waitbar(0, sprintf('Evaluating candidates 0%%'), 'Name', 'CAP detection');

    for i = 1:size(A_phase_candidates, 1)

        % Pick up the indices for current A phase candidate
        indices = A_phase_candidates{i,1};
        phase_start = indices(1);
        phase_end = indices(end);


        % Preallocate space for storing number of areas and other results
        numAreas = zeros(1, 2);  % For storing number of areas for both signals
        isEqualLength = zeros(1, 2);  % For checking the equal length condition

        % Check the candidates with CWT
        for j = 1:numel(signalNam)
            % Extract the current signal
            current_signal = signals(:,j);

            % Perform CWT
            [cwtResult, frequencies] = cwt(current_signal(phase_start:phase_end), 'amor', Fs);

            highFreqIndices = frequencies >= 0.5 & frequencies <= 30;

            cwtHighFreq = abs(cwtResult(highFreqIndices, :));

            % Create masks for threshold exceeding regions
            mask = cwtHighFreq >= thrs_cwt;

            % Find connected components in the masks
            CC = bwconncomp(mask);

            % Get the number of areas
            numAreas(j) = CC.NumObjects;

            % Get the plot length (width)
            area_length = size(mask, 2);  % Length of the plot

            % Get the properties of the connected components
            props = regionprops(CC, 'BoundingBox');

             % Extract the widths (lengths) of the bounding boxes
            widths = arrayfun(@(x) x.BoundingBox(3), props);

            % Check if any bounding box width equals the plot length
            isEqualLength(j) = any(widths == area_length);
        end

        % Now compare results for both signals
        if numAreas(1) >= 2 && numAreas(2) >= 2 || isEqualLength(1) >= 1 && numAreas(2) >= 2 || isEqualLength(2) >= 1 && numAreas(1) >= 2
            points_per_candidate(i) = points_per_candidate(i) + 1;
        end




        % Check if the emd threshold is exceeded in both channels
        rms_values = rms_bf(phase_start:phase_end);

        if any(rms_values > 90)
            points_per_candidate(i) = points_per_candidate(i) - 1;
        end




        % Check if the amplitude difference in the candidate area exceeds
        % the threshold in both channels

        % Preallocate space for storing the amplitude differences
        amplitude_differences = zeros(1,2);

        for j = 1 : numel(signalNam)
            % Extract the current signal
            current_signal = signals(:,j);

            % Calculate the amplitude difference for the current candidate
            amplitude_difference = max(current_signal(phase_start:phase_end)) - min(current_signal(phase_start:phase_end));

            % Add the amplitude difference to the preallocated variable
            amplitude_differences(j) = amplitude_difference;
        end

        % Now compare results for both signals
        if amplitude_differences(1) > amp_diff && amplitude_differences(2) > amp_diff
            points_per_candidate(i) = points_per_candidate(i) + 1;
        end





        % Check if the power threshold is exceeded in both channels

        % Get the power values for current candidate from both signals
        first_candidate_power = first_signal_power(phase_start:phase_end);
        second_candidate_power = second_signal_power(phase_start:phase_end);

        % Now compare results for both signals
        if any(first_candidate_power >= thrs_power) && any(second_candidate_power >= thrs_power)
            points_per_candidate(i) = points_per_candidate(i) + 1;
        end


        waitbar(i / size(A_phase_candidates,1), h, sprintf('Evaluating candidates %d%%', round(i / size(A_phase_candidates,1) * 100)), 'Name', 'CAP detection');


    end

    close(h);




    % Find all the A phase candidates that has 3 points and add them  to
    % the A_phases
    A_phases = {};
    for i = 1:numel(points_per_candidate)
        if points_per_candidate(i) == 3
            A_phases{end+1,1} = A_phase_candidates{i,1};
            A_phases{end, 2} = A_phase_candidates{i,2};
        end
    end



    % Extract the first element from each cell in the first column
    first_elements = cellfun(@(x) x(1), A_phases(:,1));

    % Sort the rows of A_phases based on the extracted first elements
    [~, sorted_idx] = sort(first_elements);

    % Apply the sorting index to A_phases
    A_phases = A_phases(sorted_idx, :);



    difference = 2*Fs;


    j = 1;

    while j <= length(A_phases) - 1
        current_array = A_phases{j,1};
        current_array_end = A_phases{j,1}(end);
        current_subtype = A_phases{j,2};
        next_array = A_phases{j+1,1};
        next_array_start = A_phases{j+1,1}(1);

        if (next_array_start - current_array_end) <= difference
            combined_arrays = vertcat(current_array, next_array);
            A_phases{j,2} = current_subtype;
            A_phases{j,1} = (min(combined_arrays):max(combined_arrays))';
            A_phases{j+1,1} = [];
            A_phases{j+1,2} = [];

            j = j + 2;

        else
           A_phases{j,2} = current_subtype;
           A_phases{j,1} = current_array;
           j = j + 1;

        end
    end
    % Remove empty cells
    non_empty_indices = ~cellfun('isempty', A_phases(:,1));
    A_phases = A_phases(non_empty_indices, :);



   % Calculate the CAP cycles. If the next A phase is under 60 s apart
    % from previous A phase, it is included to the cycle. Mark the spaces
    % as B phase. Minimum cycle is: A B

    cycle_limit = 60 * Fs; % Define the threshold

    cycles = {}; % Initialize the cell array to hold all cycles

    currentCycle = {};


    for i = 1:length(A_phases) - 1
        % Get the current and next cells
        currentCell = A_phases{i,1};
        nextCell = A_phases{i+1,1};

        % Get the last element of the current cell and the first element of the next cell
        lastElementCurrent = currentCell(end);
        firstElementNext = nextCell(1);


        if (firstElementNext - lastElementCurrent) > 0 && ...
            (firstElementNext - lastElementCurrent) <= cycle_limit

            B_phase =(lastElementCurrent+1:firstElementNext-1)';

            currentCycle(1,end+1) = {currentCell};
            currentCycle{2,end} = 'A';

            currentCycle{1,end+1} = B_phase;
            currentCycle{2,end} = 'B';

        else
            currentCycle(1,end+1) = {currentCell};
            currentCycle{2,end} = 'A';

            cycles{end+1} = currentCycle';
            currentCycle = {};

        end

    end



    % Remove all those cells from cycles that are too short to be a cycle.
    % In other words, remove the isolated A-phases that are non-CAP.
    for i = length(cycles):-1:1
        current_cycle = cycles{i};

        if length(current_cycle) < 3
            % Remove the element from the cell array
            cycles(i) = [];
        end
    end

    cycles = cycles';



    % now the remaining cycles represent the CAP sequences. Create a vaiable sequences to
    % distinguish the individual cycles and sequences.

    sequences = cycles;

    % Pick up the cycles from sequences for saving
    cycles = {};

    for i = length(sequences):-1:1
        current_sequence = sequences{i};
        
        for k = 1:length(current_sequence)/2
            A = (k-1)*2 + 1;
            B = A + 1;
            cycles{end+1} = current_sequence(A:B,:); 
        end
    end
    
    cycles = cycles';


    % Remove those sequences that do not meet the criteria for CAP sequence

    for i = length(sequences):-1:1
        current_cycle = sequences{i};

        if length(current_cycle) < 5
            % Remove the element from the cell array
            sequences(i) = [];
        end
    end


    % Remove the last A-phase from each sequence that is counted as non-CAP

    for i = length(sequences):-1:1
        current_sequence = sequences{i};
        current_sequence(end,:) = [];

        sequences{i} = current_sequence;
    end




    

    % Create a data table for CAPs
    event = 'A_phase';

    for i = 1:length(A_phases)
        current_event = A_phases{i,1};

        subtype = A_phases{i,2};
        starting_time = current_event(1)/Fs;
        duration = length(current_event)/Fs;


        current_struct = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'sampling_frequency', Fs, ...
                                'subtype', subtype, 'start_sec', starting_time, 'duration', duration, ...
                                'indices', current_event);

        A_phases{i} = current_struct;

    end

    event = 'sequence';

    for i = 1:length(sequences)
        current_event = sequences(i);

        starting_time = current_event{1,1}{1,1}(1)/Fs;
        ending_time = current_event{1,1}{end,1}(end)/Fs;
        duration = ending_time - starting_time;

        current_struct = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'sampling_frequency', Fs, ...
                                'start_sec', starting_time, 'duration', duration, ...
                                'Sequence', sequences(i));

        sequences{i} = current_struct;

    end


    event = 'cycle';

    for i = 1:length(cycles)
        current_event = cycles(i);

        starting_time = current_event{1,1}{1,1}(1)/Fs;
        ending_time = current_event{1,1}{2,1}(end)/Fs;
        duration = ending_time - starting_time;

        current_struct = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'sampling_frequency', Fs, ...
                                'start_sec', starting_time, 'duration', duration, ...
                                'cycle', cycles(i));

        cycles{i} = current_struct;

    end





    cap = struct('A_phases', {A_phases},'Sequences', {sequences}, 'Cycles', {cycles});

    % Save
    save(fullfile(save_path,[filename '_' 'CAP']), 'cap');
    fprintf('CAP analysis completed.\n');

end

