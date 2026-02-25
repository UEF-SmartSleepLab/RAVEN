function RAVEN_desktop(app, Fs, Fs_cap, filePath, filename, channels_to_load, channels_of_interest, channels_of_interest_cap, microevents, save_path, mode, format, Thresholds)



%%-------------------------------------------------------------------------
% ------------------------- FOLDERS FOR SAVING ----------------------------
%%-------------------------------------------------------------------------

if strcmp(mode, 'Group Analysis')
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
            newFolderNames{i} = [fileNames{i}, '_EEG_analysis_'];
        end

        % Initialise variable to save the saveFolderPaths
        saveFolderPaths = cell(length(newFolderNames),1);
        % Create a new folder path using save path and the modified name
        for i = 1 : length(newFolderNames)
            saveFolderPath = fullfile(save_path, newFolderNames{i});
            saveFolderPaths{i} = saveFolderPath;
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
            newFolderNames{i} = [filenames{i}, '_EEG_analysis_'];
        end

        % Initialise variable to save the saveFolderPaths
        saveFolderPaths = cell(length(newFolderNames),1);
        % Create a new folder path using save path and the modified name
        for i = 1 : length(newFolderNames)
            saveFolderPath = fullfile(save_path, newFolderNames{i});
            saveFolderPaths{i} = saveFolderPath;
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

     % Add the desired text "_EEG_analysis" to the name
     newFolderName = [fileName, '_EEG_analysis_'];
     % Create the new folder path using save_path and the modified filename
     saveFolderPath = fullfile(save_path, newFolderName);
end



%% 


%%-------------------------------------------------------------------------
% ------------------------------ ANALYSIS ---------------------------------
%%-------------------------------------------------------------------------



% Find the signals of interest from given file. Filter and downsample
% the signals for analysis.

if strcmp(mode, 'Individual Focus')
    filenames = {filename};
    saveFolderPaths = {saveFolderPath};
end

failedFiles = {};

% Create a waitbar
d = uiprogressdlg(app.UIFigure, ...
    'Title','RAVEN Analysis', ...
    'Message','Starting...', ...
    'Cancelable','on');
    
    for i = 1 : length(filenames)
        try 
            % Clear THIS variable in case there was an error during last
            % file
            if exist('THIS', 'var')
                clear THIS
            end
            % Update waitbar for current patient
            d.Value = (i-1)/length(filenames);
            d.Message = sprintf('Analysing patient %d of %d', i, length(filenames));
            drawnow
    
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
    
            % Get time of current run to avoid overriting existing analyse
            dt = datetime('now');
        
            % Extract components
            hourValue = hour(dt);
            minuteValue = minute(dt);
            secondValue = round(second(dt));
        
            % Combine into a formatted string and save in a variable
            timeVar = sprintf('%02d:%02d:%02d', hourValue, minuteValue, secondValue);
            timeVar = strrep(timeVar,':','_');
    
            saveFolderPath = [saveFolderPath, timeVar];
    
            % Create the folder if it doesn't exist
            if ~exist(saveFolderPath, 'dir')
                mkdir(saveFolderPath);
            end
    
    
            % (2) Signal formatting
            % Find the signals of interest from given file. Filter and
            % downsample the signals for analysis.
            if strcmp(format, 'EDF')
                loadedRawSignals = edf_mat_converter(filePath,filename,channels_to_load, d);
            elseif strcmp(format, 'SLF')
                loadedRawSignals = slf_mat_converter(filePath,filename,channels_to_load, d);
            end
    
            % Pick up the raw EEG signals according to the channels of interest
            % for artifact identification
            % Initialize a struct to hold the extracted signals and sampling frequencies
            rawEEGSignals = loadedRawSignals;



%%-------------------------------------------------------------------------
% ------------------------ DETECTION ALGORITHMS----------------------------
%%-------------------------------------------------------------------------
    
    
            if isempty(channels_of_interest)
                % disp('Only CAP analysis is performed.');
            else
                try
                    [mat_signal_filtered, rawEEGSignals] = EEG_formatting(Fs,loadedRawSignals,channels_of_interest,d);
                catch THIS
                    failedFiles{end+1,1} = filename;
                    failedFiles{end, 2} = THIS.message;
                end
            end

            if d.CancelRequested
                if isvalid(d), close(d); end
                return
            end
    
    
    
            % (3) K-complex detection 
            if iscell(microevents) && any(strcmp('kcomplexes', microevents)) && ~exist('THIS','var')
                Detect_Kcomplex(app, Fs, mat_signal_filtered, fileName, saveFolderPath, rawEEGSignals, Thresholds, d);
            end


            if d.CancelRequested
                if isvalid(d), close(d); end
                return
            end
    
    
            % (4) Spindle detection 
            if iscell(microevents) && any(strcmp('spindles', microevents)) && ~exist('THIS','var')
                Detect_Spindles(app, Fs, mat_signal_filtered, fileName, saveFolderPath, rawEEGSignals, Thresholds, d);
            end

            if d.CancelRequested
                if isvalid(d), close(d); end
                return
            end
    
    
            % (5) Delta waves detection (slow wave sleep) 
            if iscell(microevents) && any(strcmp('delta_waves', microevents)) && ~exist('THIS','var')
                Detect_delta_waves(app, Fs, mat_signal_filtered, fileName, saveFolderPath, rawEEGSignals,Thresholds, d);
            end

            if d.CancelRequested
                if isvalid(d), close(d); end
                return
            end
    
    
    
            % (6) CAP detection 
            if iscell(microevents) && any(strcmp('cap', microevents))

                try
                    mat_signal_filtered = EEG_formatting(Fs_cap,loadedRawSignals,channels_of_interest_cap, d);
                catch THIS
                    failedFiles{end+1,1} = filename;
                    failedFiles{end, 2} = THIS.message;
                end
                if ~exist('THIS', 'var')
                    Detect_CAP(app, Fs_cap, mat_signal_filtered, fileName, saveFolderPath,Thresholds, d);
                end
            end

            if d.CancelRequested
                if isvalid(d), close(d); end
                return
            end

            d.Value = i/length(filenames);

        catch ME
            failedFiles{end+1,1} = filename;
            failedFiles{end, 2} = ME.message;
        end
    end

    if ~isempty(failedFiles)
        % Create a string from both columns for comparison
        rowsAsStrings = cellfun(@(x, y) [x '||' y], failedFiles(:,1), failedFiles(:,2), 'UniformOutput', false);
        
        % Find unique rows
        [~, uniqueIdx] = unique(rowsAsStrings, 'stable');
        
        % Extract unique rows
        failedFiles = failedFiles(uniqueIdx, :);

        logFilePath = fullfile(save_path, 'failed_files_log.txt');

        fid = fopen(logFilePath, 'w');

        if fid == -1
            warning('Could not open file to write failed filenames.');
        else
            for i = 1:size(failedFiles, 1)
                fprintf(fid, '%s\t%s\n', failedFiles{i, 1}, failedFiles{i, 2});
            end
            fclose(fid);
        end
    end

    close(d);

    app.UIFigure.Pointer = 'arrow';
    drawnow;

    % Remove duplicate rows (you already did this above, but safe)
    if ~isempty(failedFiles)
        rowsAsStrings = cellfun(@(x,y) [x '||' y], failedFiles(:,1), failedFiles(:,2), 'UniformOutput', false);
        [~, uniqueIdx] = unique(rowsAsStrings, 'stable');
        failedFiles = failedFiles(uniqueIdx, :);
    end

    % ---------------------------------------------------------
    % 1) CASE: ALL FILES OK
    % ---------------------------------------------------------
    if isempty(failedFiles)

        uialert(app.UIFigure, ...
            'Analysis completed successfully. All subjects were processed.', ...
            'RAVEN Analysis Finished', ...
            'Icon','success');    
        return
    end


    % ---------------------------------------------------------
    % 2) CASE: SOME FILES FAILED
    % ---------------------------------------------------------

    nFailed = size(failedFiles,1);
    total   = length(filenames);

    summaryMsg = sprintf([ ...
        'Analysis completed with warnings.\n\n' ...
        '%d of %d subjects could not be processed.\n' ...
        'Failed subjects are skipped in results.\n\n' ...
        'You can review the list now or check the log file later.\n\n' ...
        'Do you want to see the failed files?' ], nFailed, total);

    choice = uiconfirm(app.UIFigure, summaryMsg, ...
        'RAVEN Analysis Completed with Errors', ...
        'Options', {'Show Failed Files','OK'}, ...
        'DefaultOption', 1);


    % ---------------------------------------------------------
    % 3) SHOW DETAILED LIST IF USER WANTS
    % ---------------------------------------------------------
    if strcmp(choice,'Show Failed Files')

        % Build formatted message
        detailLines = cell(nFailed,1);
        for k = 1:nFailed
            detailLines{k} = sprintf('%s\n   → %s', ...
                failedFiles{k,1}, failedFiles{k,2});
        end

        detailMsg = strjoin(detailLines, sprintf('\n\n'));

        uialert(app.UIFigure, ...
            sprintf('The following subjects failed during analysis:\n\n%s', detailMsg), ...
            'Failed Subjects', ...
            'Icon','warning');
    end
end




%%-------------------------------------------------------------------------
% ---------------------------- ALL FUNCTIONS ------------------------------
%%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------
% ------------------------ EDF CONVERSION TO MAT --------------------------
%%-------------------------------------------------------------------------



function loadedRawSignals = edf_mat_converter(Path,filename,channelsToLoad, mainDialog)


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

    % % Check if the file was found
    % if isempty(edf_file)
    %     error('File with the name "%s" was not found in the directory "%s".', filename, Path);
    % end

    % Load the .edf file 
    [~, signals, ~, ch_info] = edf_fileread(edf_file,'signal',false,channelsToLoad);

    % Check main dialog cancel
    if isvalid(mainDialog) && mainDialog.CancelRequested
        if isvalid(mainDialog), close(mainDialog); end
        error('User requested cancel. Exiting analysis.');  % propagate cancellation
    end

    fsVec = ch_info.fs(:);

    loadedRawSignals = struct();

    for c = 1:numel(channelsToLoad)
        chanName = channelsToLoad{c};
        % Ensure valid field name
        fieldName = matlab.lang.makeValidName(chanName);
    
        % Assign both signal data and sampling frequency
        loadedRawSignals.(fieldName) = struct( ...
            'data', signals{c}, ...   % or signals(:,c) if numeric matrix
            'fs',   fsVec(c) ...
        );
    end

end





% EDF reader written by Akseli Leino and Sami Nikkonen
%
% HOW TO USE
% [header, signals, events, ch_info] = readedf_file(filename, filetype, printopt, targetSignals)
% 
% OUTPUTS
% header = full edf header as struct
% signals = cell array of stored signals (empty in event-only files)
% events = cell array of scored events (empty in signal-only files)
% ch_info = most relevant info from header for targetSignals
% 
% INPUTS
% filename = filename of the edf file
% filetype = signal/event/both, whether the EDF-file is a signal-only file, event-only file, or contains both signals and annotations (defaults to both)
% printopt = true/false, toggles printing info to command window (defaults to false)
% targetSignals = names of the target signals (if not specified, all signals will be read)

function [header, signals, events, ch_info] = edf_fileread(filename, filetype, printopt, targetSignals)

fname=filename;

% Handle omitted input variables
if ~exist('targetSignals', 'var')
    targetSignals = {};
elseif isempty(targetSignals)
    targetSignals = {};
end

if ~exist('printopt', 'var')
    printopt=false;
elseif isempty(printopt)
    printopt=false;
elseif ~islogical(printopt)
    error('Invalid printopt. Must be true/false.')
end

if ~exist('filetype', 'var')
    filetype='both';
elseif isempty(filetype)
    filetype='both';
end

if strcmp(filetype,'signal')
    ftype_sig=true;
    ftype_ev=false;
elseif strcmp(filetype,'event')
    ftype_sig=false;
    ftype_ev=true;
elseif strcmp(filetype,'both')
    ftype_sig=true;
    ftype_ev=true;
else
    error('Invalid filetype. Accepted types: signal,event,both')
end

% Open file
[fid,msg] = fopen(fname, 'r', 'l', 'ISO-8859-1');
if fid == -1
    error(msg)
end

% Header will be always read first 
[ch_info, header, targetSignals] = get_header(fid, fname, targetSignals, printopt);

% Declare signals and events first as empty so they wont be missing if not read
events=[];
signals=[];

% Read either signals or events or both depending on file type
if ftype_sig
    [signals, annotationsStart] = get_signals(fid, header, targetSignals, printopt);   
end

% If signal not read, i.e. the file is annotations only file, set
% annotation start as empty
if ~exist('annotationsStart', 'var')
    annotationsStart=[];
end
if ftype_ev
    events = get_events(fid,annotationsStart);
end

end

%%
% The function reads the header bytes from the edf
% Input:
% fid: The file id of the edf file
% targetSignals: cell array of selected signals e.g. {'C4', 'E2'}
% printopt = toggles printing info to command window
%
% Output:
% ch_info: struct containing main signal information
% header: all header information
% targetSignals: array of the indices of the selected signals

function [ch_info, header, targetSignals] = get_header(fid, fname, targetSignals, printopt)
% Get header info not dependent on signal/record/sample size
header.ver        = str2double(char(fread(fid,8)'));
header.patientID  = fread(fid,80,'*char')';
header.recordID   = fread(fid,80,'*char')';
header.startdate  = fread(fid,8,'*char')';
header.starttime  = fread(fid,8,'*char')';
header.bytes      = str2double(fread(fid,8,'*char')');
reserved          = sprintf(fread(fid,44,'*char'));
header.reserved   = reserved; % also save to header
header.records    = str2double(fread(fid,8,'*char')');
header.duration   = str2double(fread(fid,8,'*char')');
header.ns         = str2double(fread(fid,4,'*char')');

% Check file validity
if header.records == -1
    error(['Invalid file - File not closed when recording was finished. ' ...
        'The record number is set to -1, which should have been changed ' ...
        'immediately when the measurement was completed.']);
elseif isnan(header.ns)
    error('Invalid EDF header. Number of signals is not a number.')
end

% Check EDF version
if printopt
    disp(['Filename: ', fname])
    if contains(reserved, 'EDF+C')
        fprintf('EDF type: EDF compatible EDF+ \n')
    elseif contains(reserved, 'EDF+D')
        fprintf('EDF type: Discontinuous EDF file \n')
    else
        fprintf('EDF type: Regular EDF file \n')
    end
end

% Get signal labels
for ii = 1:header.ns
    header.label{ii} = regexprep(fread(fid,16,'*char')','\W','');
end

% Select either requested signals or all signals included if none is given
if isempty(targetSignals)
    targetSignals = 1:numel(header.label);
elseif iscell(targetSignals)||ischar(targetSignals)
    targetSignals = find(ismember(header.label,regexprep(targetSignals,'\W','')));

    if isempty(targetSignals) && size(header.label, 1) > 1
        error('The requested signal(s) were not detected. \n')
    end
end

% Transducer type
for ii = 1:header.ns
    header.transducer{ii} = fread(fid,80,'*char')';
end
% Physical dimension
for ii = 1:header.ns
    header.units{ii} = fread(fid,8,'*char')';
end
% Physical minimum
for ii = 1:header.ns
    header.physicalMin(ii) = str2double(fread(fid,8,'*char')');
end
% Physical maximum
for ii = 1:header.ns
    header.physicalMax(ii) = str2double(fread(fid,8,'*char')');
end
% Digital minimum
for ii = 1:header.ns
    header.digitalMin(ii) = str2double(fread(fid,8,'*char')');
end
% Digital maximum
for ii = 1:header.ns
    header.digitalMax(ii) = str2double(fread(fid,8,'*char')');
end
% Filtering
for ii = 1:header.ns
    header.prefilter{ii} = fread(fid,80,'*char')';
end
% Number of samples
for ii = 1:header.ns
    header.samples(ii) = str2double(fread(fid,8,'*char')');
end
% Second reserve
for ii = 1:header.ns
    header.reserved2 = fread(fid,32,'*char')';
end

% Save also most relevant info for requested signals
ch_info.label = header.label(targetSignals);
ch_info.label = regexprep(ch_info.label,'\W','');
ch_info.fs = header.samples(targetSignals)./header.duration;
ch_info.units = header.units(targetSignals);
ch_info.prefilter = header.prefilter(targetSignals);
end

%%
% Reads the signals of the edf file.
% Input:
% fid: The file id of the edf file
% header: all header information got from get_header
% targetSignals: indices of the selected signals
% printopt = toggles printing info to command window
%
% Output:
% signals: the stored signals in a cell
% annotationsStart: the startpoint of annotations

function [signals, annotationsStart] = get_signals(fid, header, targetSignals, printopt)
% In EDFs, the signals are stored in the following format. Signals are
% divided into small parts, which duration is specified by
% header.duration. Each record (R1-R6) consists of header.duration
% sized part of all signals (S1-S5) which number is specified by
% header.ns. Each record of signal(ii) consists of header.samples(ii)*2
% sample points. Each sample point is presented by 2 bytes (therefore
% the *2). To skip a whole record, one must skip sum(header.samples)*2
% number of bytes.
%
%    S1 S2 S3 S4 S5
% R1 S  S  S  S  S
% R2 i  i  i  i  i
% R3 g  g  g  g  g
% R4 n  n  n  n  n
% R5 a  a  a  a  a
% R6 l  l  l  l  l

% These will be used to locate signal records later
startOfRecordings = ftell(fid);
sumOfSampleBytes = sum(header.samples)*2;
targetSamples = header.samples(targetSignals);

% For signal scaling
scalefac = (header.physicalMax - header.physicalMin)./(header.digitalMax - header.digitalMin);
dc = header.physicalMax - scalefac .* header.digitalMax;
scalefac = scalefac(targetSignals);
dc = dc(targetSignals);

% Preallocate the signals to avoid reallocation during the reading of the file
signals = cell(size(targetSignals, 2), 1);

for iSignal = 1:size(signals, 1)
    signals{iSignal} = zeros(1, header.records*header.samples(targetSignals(iSignal)));
end

% Preallocate annotations startpoint in file as empty
annotationsStart=[]; 

lineLength = 0;
for iSignal = 1:size(targetSignals, 2)
    if printopt
        fprintf(repmat('\b', 1, lineLength));
        lineLength = fprintf(['Reading signal number: ', num2str(iSignal), '/', num2str(size(targetSignals, 2))]);
    end

    if targetSignals(iSignal) == 1
        startOfSignal = startOfRecordings;
    else
        startOfSignal = startOfRecordings + sum(header.samples(1:targetSignals(iSignal) - 1)*2);
    end

    % Find start spots of annotations
    if strcmp(header.label(iSignal),'EDFAnnotations')
        % Append for each found annotation signal
        % Only the first one should be needed but all are saved anyway
        % Perforamnce impact on matrix size changes should be negligible as the mat is so small
        annotationsStart=[annotationsStart;startOfSignal];
    end

    fseek(fid, startOfSignal, 'bof');

    % A signal is read record by record into 'signals'. After a signal
    % segment is read, the loop will skip a whole record forward to
    % capture the next record of the corresponding signal.

    for iRecord = 1:header.records
        signals{iSignal}((iRecord-1)*targetSamples(iSignal) + 1: ...
            iRecord*targetSamples(iSignal)) ...
            = fread(fid, targetSamples(iSignal), 'int16') * scalefac(iSignal) + dc(iSignal);

        fseek(fid, sumOfSampleBytes-targetSamples(iSignal)*2, 'cof');
    end

end
if printopt
    fprintf('\n')
end
end

%%
% Gets the events of an EDF+ file. This functionality should work but is not tested very much.
% Input:
% fid: file id of the event file
% annotationsStart: the startpoint of annotations
%
% Output:
% events: all events or annotations

function events = get_events(fid,annotationsStart)
% fprintf('Reading annotations \n');

% Find startpoint of annotations
if isempty(annotationsStart)
    % If there are no annotations found or file is annotations only
    startOfEvents = ftell(fid);
else
    % All annotions should be in a row so should be ok to just read from
    % first found annotation. Most files will only have one annotation channel anyway
    startOfEvents = annotationsStart(1);
end

fseek(fid, startOfEvents, 'bof');

% The events are formatted as following:
%
% start_of_recording{20}{20}description{20}{0}
% onset{21}event_duration{20}description{20}{0}
%
% Where event onsets are seconds after start_of_recording and can be
% negative or positive numbers. {x} represents the corresponding ASCII
% character. The event_duration may be omitted, which will also omit
% the leading {21} character. Textscan function is not very versatile
% with missing information, when also the delimiter of the
% corresponding column is missing. Therefore, the events must be
% scanned in format float-string-string-string. The first column will
% always be onset, and may therefore be scanned as float. The second
% column is usually duration, except when the duration is ommited. In
% this situation the textscan moves already to the next field and will
% put the description in the duration column (reading this as float
% would cause an error). This must be corrected afterwards.

events = textscan(fid, '%.3f%s%s%s', 'EndOfLine', char(0), 'Delimiter', [char(21), char(20)], 'EmptyValue', NaN);
empty_descriptions = cellfun('isempty', events{3});

% Move the desciption to correct column
events{3}(empty_descriptions) = events{2}(empty_descriptions);
events{2}(empty_descriptions) = {NaN};

% Duration to double
events{2} = str2double(events{2});

% Add descriptions to event cell for clarity
% Technically unnecessary as this is same for all EDFs
% as annotations should always be stored in this order in the EDFs
titles = cell(1,size(events,2));  
titles{1}='onset';
titles{2}='duration';
titles{3}='label';

events=[titles;events];

end






%%-------------------------------------------------------------------------
% ---------------------- SLF TO MAT CONVERSION ----------------------------
%%-------------------------------------------------------------------------


function loadedRawSignals = slf_mat_converter(Path, filename, channels_to_load, mainDialog)

    loadedRawSignals = struct();   % FINAL output directly

    items = dir(Path);

    matchingFolder = items([items.isdir] & strcmp({items.name}, filename));

    % if isempty(matchingFolder)
    %     error('No folder with the name "%s" found in "%s".', filename, Path);
    % end

    specificFolderPath = fullfile(Path, matchingFolder(1).name);

    subItems = dir(specificFolderPath);
    subFolders = subItems([subItems.isdir] & ~ismember({subItems.name}, {'.', '..'}));
    subFolderNames = {subFolders.name};

    % Find matching channel folders
    matchingSubFolders = subFolderNames(contains(subFolderNames, channels_to_load));

    for i = 1:length(matchingSubFolders)

        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end

        subFolderName = matchingSubFolders{i};
        currentSubFolderPath = fullfile(specificFolderPath, subFolderName);

        npyFiles  = dir(fullfile(currentSubFolderPath,'*.npy'));
        jsonFiles = dir(fullfile(currentSubFolderPath,'*.json'));

        if length(npyFiles) ~= 1
            error('The %s folder must contain exactly one .npy file.', subFolderName);
        end

        if length(jsonFiles) ~= 1
            error('The %s folder must contain exactly one .json file.', subFolderName);
        end

        dataFilePath      = fullfile(currentSubFolderPath, npyFiles.name);
        attributeFilePath = fullfile(currentSubFolderPath, jsonFiles.name);

        data = readNPY(dataFilePath);

        attributes = jsondecode(fileread(attributeFilePath));
        samplingRate = attributes.sampling_rate;

        % ---- STRUCT STORAGE (instead of cell + eval) ----
        chanName  = subFolderName;   % channel name comes from folder
        fieldName = matlab.lang.makeValidName(strtrim(chanName));

        loadedRawSignals.(fieldName) = struct( ...
            'data', data, ...
            'fs',   samplingRate ...
        );
    end
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




%%-------------------------------------------------------------------------
% ---------------------- FORMATTING THE SIGNALS ---------------------------
%%-------------------------------------------------------------------------




function [mat_signal_filtered, rawEEGSignals] = EEG_formatting(fs_new, rawSignals, channels_of_interest, mainDialog)


    % Construct selected signals (with reference subtraction if needed)
    loadedSignals = struct();
    
    for s = 1:numel(channels_of_interest)
        sigName = channels_of_interest{s};
        
        if contains(sigName, '-')
            parts = split(sigName, '-');
            chan1 = parts{1};
            chan2 = parts{2};
    
            % Subtract signals for referenced channel
            sigVec = rawSignals.(chan1).data - rawSignals.(chan2).data;
    
            % Choose a fs (if you want, just take fs of first channel)
            fsVal = rawSignals.(chan1).fs;
        else
            sigVec = rawSignals.(sigName).data;
            fsVal  = rawSignals.(sigName).fs;
        end
    
        % Convert field name to safe struct field
        fieldName = strrep(sigName, '-', '_');
    
        % Store both data and fs in a sub-struct
        loadedSignals.(fieldName) = struct('data', sigVec, 'fs', fsVal);
    end


    sigNames = fieldnames(loadedSignals);
    origFs = loadedSignals.(sigNames{1}).fs;  % take Fs from first channel

    hpFilt = designfilt( ...
        'highpassiir', ...
        'FilterOrder', 3, ...
        'PassbandFrequency', 0.1, ...
        'PassbandRipple', 0.001, ...
        'SampleRate', origFs, ...
        'DesignMethod', 'cheby1');

    lpFilt = designfilt( ...
        'lowpassiir', ...
        'FilterOrder', 20, ...
        'PassbandFrequency', 30, ...
        'PassbandRipple', 0.001, ...
        'SampleRate', origFs, ...
        'DesignMethod', 'cheby1');

    % ------------------------------------------------------------
    % Apply filtering + resampling
    % ------------------------------------------------------------
    mat_signal_filtered = struct();

    for i = 1:numel(sigNames)
        sigName = sigNames{i};
        channel = loadedSignals.(sigName).data;

        if isrow(channel)
            channel = channel.';
        end

        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end

        channel_out = filter_and_resample_EEG( ...
            channel, ...
            lpFilt, ...
            hpFilt, ...
            fs_new, ...
            origFs);

        mat_signal_filtered.(sigName) = channel_out;
    end

    % ------------------------------------------------------------
    % Resample the raw EEG
    % ------------------------------------------------------------

    rawEEGSignals = struct();

    for i = 1:numel(sigNames)
        sigName = sigNames{i};
        channel = loadedSignals.(sigName).data;

        if isrow(channel)
            channel = channel.';
        end

        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end

        channel_out = resample_raw_EEG( ...
            channel, ...
            fs_new, ...
            origFs);
        
        rawEEGSignals.(sigName) = channel_out;
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
% Ratio for the resample function
[p,q] = rat(fs_new / sampling_frequency);
%resample
resampled_channel = resample(channel,p,q);
end


function [resampled_channel] = resample_raw_EEG(channel,fs_new,sampling_frequency)
if isa(channel, 'single')
    channel = double(channel);
end
% Ratio for the resample function
[p,q] = rat(fs_new / sampling_frequency);
%resample
resampled_channel = resample(channel,p,q);
end









%%-------------------------------------------------------------------------
% ------------------------- K-COMPLEX DETECTION ---------------------------
%%-------------------------------------------------------------------------




function kComplexes = Detect_Kcomplex(app, Fs, mat_signal_filtered, filename, save_path, rawEEGSignals, Thresholds, mainDialog)




channels_to_analyse = fieldnames(mat_signal_filtered);


% Default thresholds
thrs_power = Thresholds.Power_Kcomplex;           % Delta band
thrs_rms = Thresholds.RMS_Kcomplex;               % Filtered signal (base frequency)
thrs_envelope = Thresholds.Envelope_Kcomplex;     % Delta band
thrs_cwt = Thresholds.CWT_Kcomplex;               % Filtered signal (base frequency)
thrs_corr = Thresholds.Correlation_Kcomplex;      

% Alternative annotation name
aasmEvent_name = 'UNSURE';

% Define column names
% columnNames = channels_to_analyse;

% Create an empty cell array with specified column names
dataTableKcomplex = cell(1, numel(channels_to_analyse)); % Start with one row

% Assign column names to the first row of the cell array
dataTableKcomplex(1, :) = channels_to_analyse;

event = 'Kcomplex';


% ------------ Create canonical Kcomplexes ---------- %

durations = [1.0, 1.2, 1.5, 2.0];         % include durations s
p200_flags = [true, false];               % P200 sometimes absent
n550_widths = [0.10, 0.13, 0.18, 0.22];
p900_widths = [0.10, 0.15, 0.22];
p200_widths = [0.02, 0.04, 0.06];

idx = 0; templates = struct();
for d = durations
  for pf = p200_flags
    for nw = n550_widths
      for pw = p900_widths
        for p2w = p200_widths
          idx = idx + 1;
          [t, w] = makeKcomplexTemplate(Fs, d, ...
              0.15, p2w, 0.20, ...   % P200: small amp, latency ~200ms
              1.0, nw, 0.55, ...     % N550: main neg, latency ~550ms
              0.5, pw, 0.90, ...     % P900: late positive ~900ms
              pf);
          templates(idx).t = t;
          templates(idx).wave = w;
          templates(idx).meta.duration = d;
          templates(idx).meta.includeP200 = pf;
          templates(idx).meta.n550_width = nw;
          templates(idx).meta.p900_width = pw;
        end
      end
    end
  end
end



% ---------- ANALYSE THE CHANNELS OF INTEREST ---------- %

for channel = 1 : length(channels_to_analyse)

    chan = cell2mat(channels_to_analyse(channel));
    EEGbf = mat_signal_filtered.(chan);
    rawSignal = rawEEGSignals.(chan);


    % Window size and window step for caluculating powers in signal
    WinL = round(0.4 * Fs);  % Window length (0.2 seconds)
    WinS = floor(0.2 * Fs);  % Step size in samples (0.1 seconds)

    % Filter the data for delta band
    f_low = 0.1; % Lower cutoff frequency
    f_high = 4.5; % Upper cutoff frequency
    order = 2; % Filter order

    % Design the bandpass filter using the Butterworth filter 
    [b, a] = butter(order, [f_low, f_high]/(Fs/2), 'bandpass');

    % Apply the bandpass filter to the signal
    EEG_delta = filter(b, a, EEGbf);

    % Smooth the EEG_delta with moving average filter
    EEG_delta_smooth = movmean(EEG_delta, 0.1*Fs);


% ------------ CALCULATE DELTA POWER ------------ %

    delta_power = zeros(size(EEG_delta));

    padL = floor((0.1/2)*Fs); % Zero padding length for the window

    for i = 1:WinS:length(EEG_delta) - WinL + 1
         % Define window indices
        Wstart = max(1, i - padL);
        Wend   = min(length(EEG_delta), i + WinL - 1);
    
        % Extract and pad window symmetrically if needed
        window = EEG_delta(Wstart:Wend);
        if length(window) < WinL
            window = padarray(window,padL,'pre'); % pad at the end if short
        end
    
        % Compute power in the window
        pw = sum(window.^2) / WinL;
    
        % Assign the same value to all samples within the window
        delta_power(i:i+WinL-1) = max(log10(pw));
    end



% ------------ CALCULATE RMS ------------ %

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

    % Check main dialog cancel
    if isvalid(mainDialog) && mainDialog.CancelRequested
        if isvalid(mainDialog), close(mainDialog); end
        error('User requested cancel. Exiting analysis.');  % propagate cancellation
    end







% ------------ KCOMPLEX DETECTION AND EVALUATION ------------ %


    % Detect kcomplex candidates with point system using thresholds. If all
    % thresholds are exceeded Kcomplex is most likely detected. kcomplex should
    % be at least 0.5 s in length.



    % Update the WinL and make the window overlap to be 50%
    WinL = round(15*Fs);
    overlap = round(0.5*WinL);

    % Calculate the total number of windows
    total_windows = floor((length(EEGbf) - WinL) / (WinL - overlap)) + 1;

    kComplexes = cell(total_windows, 1);

    dd = uiprogressdlg(app.UIFigure, ...
    'Title','K-complex Analysis', ...
    'Message','Detecting Kcomplexes over night...', ...
    'Cancelable','on');

    ripple_band_high = [20 30];
    order = 2;

    for i = 1:total_windows
        % Calculate start and end indices for the current window
        startIndex = 1 + (i - 1) * (WinL - overlap);
        endIndex = startIndex + WinL - 1;

        % Ensure endIndex does not exceed the length of the signal
        endIndex = min(endIndex, length(EEGbf));

        % Extract a windows from different signal
        power_window = delta_power(startIndex:endIndex);


        % Find indices over the threshold
        indices_over_threshold = find(power_window >= thrs_power);

        if isvalid(dd) && dd.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end

        % Check if there are any consecutive runs within the desired length range
        if isempty(indices_over_threshold)
            % If no consecutive indices, break out of the loop and move to the next iteration
            dd.Value = i/total_windows;
            drawnow
            continue;
        end

        % Find consecutive indices
        candidate_indices = find_consecutive_kcomplexes(indices_over_threshold,Fs);
        % Check if there are any consecutive runs within the desired length range
        if isempty(candidate_indices)
            % If no consecutive indices, break out of the loop and move to the next iteration
            dd.Value = i/total_windows;
            drawnow
            continue;
        end




        % CALCULATING THE CWT FROM broad band 
        [cwtResult, frequencies] = cwt(EEGbf(startIndex:endIndex), 'amor', Fs);

        % Find regions where CWT values exceed the threshold
        lowFreqIndices = frequencies >= 0.1 & frequencies <= 4.5;
        cwtLowFreq = abs(cwtResult(lowFreqIndices, :));
        % Create masks for threshold exceeding regions
        mask = cwtLowFreq > thrs_cwt;


        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end



        % Map indices from power_window to kComplex_power for each consecutive run
        filtered_indices = cell(size(candidate_indices));

        for j = 1:length(candidate_indices)
            % Map indices from power_window
            real_indices = candidate_indices{j} + startIndex - 1;

            first_index = real_indices(1,1);

            % Calculate starting point in seconds
            starting_time_seconds = first_index / Fs;

            % Store filtered indices and starting time as a struct or tuple
            filtered_indices{j} = struct('indices', real_indices, 'start_sec', starting_time_seconds);
        end





        event_indices = cell(size(filtered_indices));

        % Pick up the indices to speed up the process
        for j = 1:length(filtered_indices)
            event_indices{j} = filtered_indices{j,1}.indices;
        end




        % ------------ CANDIDATE EVALUATION ------------ %

        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end


        % Initialize points for the current window
        points_per_candidate = zeros(size(event_indices));

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
                points_per_candidate(j) = points_per_candidate(j) + 1;
            end

            % Process the result based on the rms values
            if any(subset_rms_bf >= 180)
                points_per_candidate(j) = points_per_candidate(j) - 1;
            end
        end


        
        % Calculate the mean FFT to remove artefacts
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
                points_per_candidate(j) = points_per_candidate(j) - 3;
            end 

        end



        % Check if there are areas that exceeds the CWT threshold.
        for j = 1:size(event_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = event_indices{j}(1) - startIndex + 1;
            end_index = event_indices{j}(end) - startIndex + 1;
            timeIndices = start_index:end_index;
            % Pick up the values from cwt results for candidate area
            candidate_areas = bwconncomp(mask(:,timeIndices));
            % Initialize a variable for the number of areas in the specified time range
            numAreasInRange = candidate_areas.NumObjects;
            % Check if values in candidate area exceeds the cwt
            % threshold
            if numAreasInRange > 0
                % Process the result based on the average value
                points_per_candidate(j) = points_per_candidate(j) + 1;
            end
        end






         % Check the amplitude differrence
         for j = 1:size(event_indices, 1)
             % Get start and end indices of the current consecutive run
             start_index = event_indices{j}(1);
             end_index = event_indices{j}(end);

             signal_part = EEGbf(start_index:end_index);
             % Calculate the peak-to-peak value
             peak_to_peak = max(signal_part) - min(signal_part);
             threshold = 70; % According to the AASM

             % Check if any value in envelope is over the threshold
             if peak_to_peak >= threshold
                 points_per_candidate(j) = points_per_candidate(j) + 1;
             end
         end



         
         % Canonical Kcomplex comparison
         for j = 1:numel(event_indices)
            idx = event_indices{j};
            seg = EEGbf(idx(1):idx(end));
            seg = seg - mean(seg);
            seg = seg / (max(abs(seg)) + eps);  % match template amplitude
            % Scale the segment to start from 0
            seg = seg - seg(1);
                                
            best_corr = 0;
        
            for k = 1:numel(templates)
                tmpl = templates(k).wave;
                % resample segment to template length
                seg_resampled = interp1(1:length(seg), seg, ...
                    linspace(1,length(seg), length(tmpl)), 'pchip', 'extrap');

                corrval = corr(seg_resampled(:), tmpl(:));
        
                if corrval > best_corr
                    best_corr = corrval;
                end
            end
            
            % Check the correlation between best template and the candidate
            if best_corr > thrs_corr  
                points_per_candidate(j) = points_per_candidate(j) + 2;
            else
                points_per_candidate(j) = points_per_candidate(j) - 3;
            end
        end




         % Peak detection based on Hilbert transform
         for j = 1 : size(event_indices,1)  
             numAreas = zeros(1,3);            
             start_index = event_indices{j}(1);
             end_index = event_indices{j}(end);

             extension = 5*Fs;

             safe_start = max(1, start_index - extension);
             safe_end = min(length(EEGbf), end_index + extension);

             [b2,a2] = butter(order, ripple_band_high/(Fs/2),'bandpass');
             filt_bf_high = filtfilt(b2, a2, EEGbf(safe_start:safe_end));

             env_bf = abs(hilbert(EEGbf(safe_start:safe_end)));
             env_delta = abs(hilbert(EEG_delta_smooth(safe_start:safe_end)));
             env_bf_high = abs(hilbert(filt_bf_high));
            
             baseline_mean1 = mean(env_bf);
             baseline_std1  = std(env_bf);
             thr1_global    = baseline_mean1 + 2.2*baseline_std1;

             baseline_mean1_high = mean(env_bf_high);
             baseline_std1_high  = std(env_bf_high);
             thr1_global_high    = baseline_mean1_high + 1.6*baseline_std1_high;
                
             baseline_mean2 = mean(env_delta);
             baseline_std2  = std(env_delta);
             thr2_global    = baseline_mean2 + thrs_envelope*baseline_std2;
         

             signal_part_bf = EEGbf(start_index:end_index);
             signal_part_delta = EEG_delta_smooth(start_index:end_index);             

             filt_event_bf_high = filtfilt(b2,a2, signal_part_bf);
    
             env1 = abs(hilbert(signal_part_bf));
             env1_high = abs(hilbert(filt_event_bf_high));
             env2 = abs(hilbert(signal_part_delta));

             rippleMask1 = env1 > thr1_global;  
             rippleMask1_high = env1_high > thr1_global_high;
             rippleMask2 = env2 > thr2_global;

             % --- Cut 5% from both ends ---
             cut_fraction = 0.05;  % 5%
             N = length(rippleMask2);
            
             start_cut = round(cut_fraction * N) + 1;
             end_cut   = round((1 - cut_fraction) * N);
            
             % Keep only the middle 90%
             trimmedMask = rippleMask2(start_cut:end_cut);
             trimmedMaskBf = rippleMask1(start_cut:end_cut);
            
             CCripple1 = bwconncomp(trimmedMaskBf);
             CCripple2 = bwconncomp(trimmedMask);
             CCripple3 = bwconncomp(rippleMask1_high);
                
             % Get the number of areas
             numAreas(1) = CCripple1.NumObjects;
             numAreas(2) = CCripple2.NumObjects;
             numAreas(3) = CCripple3.NumObjects;


             if numAreas(1) == 1
                 points_per_candidate(j) = points_per_candidate(j) + 1;
             end

             if numAreas(2) == 1 || numAreas(2) == 2
                 points_per_candidate(j) = points_per_candidate(j) + 2;
             end

             if numAreas(1) == 0 && numAreas(2) == 0
                 points_per_candidate(j) = points_per_candidate(j) - 1;
             end


             if numAreas(2) > 2
                 points_per_candidate(j) = points_per_candidate(j) - 1;
             end

             if numAreas(3) >= 2
                 points_per_candidate(j) = points_per_candidate(j) - 1;
             end


         end         



   

         % ------------ POINT CALCULATION ------------ %

         % Check if any consecutive part has accumulated six points
         if any(points_per_candidate >= 5)
             % Find indices of consecutive parts with six points
             indices_six_points = find(points_per_candidate >= 5);

             % Create a cell array to store Kcomplex indices for the current window
             kComplex_candidates = cell(length(indices_six_points), 1);

             % Iterate over indices and store corresponding segments in current_kcomplex_indices
             for k = 1:length(indices_six_points)
                 current_consecutive_part = filtered_indices{indices_six_points(k)};
                 kComplex_indices = current_consecutive_part.indices(1,1):current_consecutive_part.indices(end,1);
                 current_kcomplex_starting_time = current_consecutive_part.start_sec;
                 current_kcomplex_duration = length(kComplex_indices) / Fs;

                 signal_part = EEGbf(kComplex_indices);
                 [pks_min_neg, locs_min] = findpeaks(-signal_part, 'MinPeakProminence', 20);
                pks_min = -pks_min_neg;
                
                if isempty(locs_min)
                    [~, locs_min] = min(signal_part);   % global min
                    pks_min = signal_part(locs_min);
                end
                
                [~, idx_deep] = min(pks_min);
                min_loc = locs_min(idx_deep);
                [~, locs_max] = findpeaks(signal_part, 'MinPeakProminence', 10);
                W = round(0.5 * length(signal_part));

                % --- Left side ---
                idxL = locs_max < min_loc & locs_max >= min_loc - W;
                if any(idxL)
                    % nearest maximum
                    [~, rel] = min(abs(locs_max(idxL) - min_loc));
                    left_loc = locs_max(idxL);
                    left_loc = left_loc(rel);
                else
                    % fallback: take left window boundary
                    left_loc = max(1, min_loc - W);
                end
                
                % --- Right side ---
                idxR = locs_max > min_loc & locs_max <= min_loc + W;
                if any(idxR)
                    % nearest maximum
                    [~, rel] = min(abs(locs_max(idxR) - min_loc));
                    right_loc = locs_max(idxR);
                    right_loc = right_loc(rel);
                else
                    % fallback: take right window boundary
                    right_loc = min(length(signal_part), min_loc + W);
                end

                slope_width = (right_loc - left_loc)/Fs;

                 % Calculate the peak-to-peak value
                 peak_to_peak = max(signal_part) - min(signal_part);
                 kComplex_candidates{k} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', chan, 'sampling_frequency', Fs, ...
                                                    'start_sec', current_kcomplex_starting_time, 'duration', current_kcomplex_duration, ...
                                                    'indices', kComplex_indices, 'amplitude', peak_to_peak, 'slope_width_seconds', slope_width);
             end

             % Store Kcomplex indices for the current window in
             % kcomplex_candidates
             kComplexes{i} = kComplex_candidates;

             if isvalid(dd) && dd.CancelRequested
                if isvalid(mainDialog), close(mainDialog); end
                error('User requested cancel. Exiting analysis.');  % propagate cancellation
            end

         end

        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end

         % Update waitbar
         dd.Value = i/total_windows;
         drawnow

     end

     close(dd);

     % Remove empty cells from kcomplex (indices)
     non_empty_indices = ~cellfun('isempty', kComplexes);
     kComplexes = kComplexes(non_empty_indices);



     % ------------ OVERLAP CORRECTION ------------ %


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


     dataTableKcomplex{2,channel} = kComplexes;
end

    % Save
    save(fullfile(save_path,[filename '_' 'kcomplexes']), 'dataTableKcomplex');
end




% ------------ FUNCTIONS USED IN KCOMPLEX DETECTION ------------ %





function filtered_consecutive_indices = find_consecutive_kcomplexes(indices, Fs)
    % This function identifies and filters consecutive runs of indices

    % Initialize variables
    consecutive_indices = {};
    current_run = [];

    % Initialize the minimum and maximum length of consecutive runs in data points
    min_length = 0.5 * Fs;
    max_length = 3 * Fs;

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


function [t, tmpl] = makeKcomplexTemplate(Fs, duration, ...
    p200_amp, p200_width, p200_latency, ...
    n550_amp, n550_width, n550_latency, ...
    p900_amp, p900_width, p900_latency, ...
    includeP200)
% makeKcomplexTemplate 3-component K-complex template
%
% Inputs:
%   Fs            - sampling rate (Hz)
%   duration      - total template length (s) (should be >= max latency)
%   p200_amp      - amplitude of early positive (set small, e.g. 0.2). sign positive
%   p200_width    - std (s) of early positive (e.g. 0.03)
%   p200_latency  - latency (s) of early positive (e.g. 0.20)
%   n550_amp      - amplitude of main negative (use positive number, function adds negative sign) (e.g. 1.0)
%   n550_width    - std (s) of main negative (e.g. 0.13)
%   n550_latency  - latency (s) of main negative (e.g. 0.55)
%   p900_amp      - amplitude of late positive (e.g. 0.5)
%   p900_width    - std (s) of late positive (e.g. 0.15)
%   p900_latency  - latency (s) of late positive (e.g. 0.90)
%   includeP200   - boolean, whether to include the early P200 (true/false)
%
% Outputs:
%   t    - time vector (s)
%   tmpl - normalized template (max abs = 1)
%
% Example:
%  [t, tmpl] = makeKcomplexTemplate(256, 1.2, 0.2, 0.03, 0.2, 1.0, 0.13, 0.55, 0.5, 0.15, 0.90, true);

if nargin < 11, includeP200 = true; end

% time vector (ensure duration covers the latencies)
t = linspace(0, duration, round(duration*Fs));

% build lobes (Gaussians). n550 is negative lobe (we create as negative)
p200 = zeros(size(t));
if includeP200
    p200 = p200_amp .* exp(-((t - p200_latency).^2) / (2 * p200_width^2));
end

n550 = - n550_amp .* exp(-((t - n550_latency).^2) / (2 * n550_width^2));

p900 = p900_amp .* exp(-((t - p900_latency).^2) / (2 * p900_width^2));

% combine
tmpl = p200 + n550 + p900;

% baseline remove and normalize
tmpl = tmpl - mean(tmpl);
tmpl = tmpl / (max(abs(tmpl)) + eps);
end








%%-------------------------------------------------------------------------
% ------------------------- SPINDLE DETECTION ----------------------------- 
%%-------------------------------------------------------------------------





function spindles = Detect_Spindles(app, Fs, mat_signal_filtered, filename, save_path, rawEEGSignals, Thresholds, mainDialog)

% Default thresholds
thrs_sigma = Thresholds.Power_spindle;             % Sigma band power
thrs_relative = Thresholds.RelativePower_spindle;  % Relative power (sigma band - delta band)
thrs_envelope = Thresholds.Envelope_spindle;       % Sigma band
thrs_cwt = Thresholds.CWT_spindle;                 % Sigma band
thrs_rms_sigma = Thresholds.RMS_spindle;           % Sigma band
thrs_P2P = Thresholds.P2P_spindle;                 % Sigma band

% Alternative annotation name
aasmEvent_name = 'UNSURE';

channels_to_analyse = fieldnames(mat_signal_filtered);

% Create a table where to save all detected spindles from different
% channels

% Define column names
% columnNames = channels_to_analyse;

% Create an empty cell array with specified column names
dataTableSpindle = cell(1, numel(channels_to_analyse)); % Start with one row

% Assign column names to the first row of the cell array
dataTableSpindle(1, :) = channels_to_analyse;

event = 'spindle';

% Clip the signals to be in range -100,100
for channel = 1:length(channels_to_analyse)
    chan = cell2mat(channels_to_analyse(channel));
    clipped_signal = clip(mat_signal_filtered.(chan), -100, 100);
    mat_signal_filtered.(chan) = clipped_signal;
end



for channel = 1 : length(channels_to_analyse)

    chan = cell2mat(channels_to_analyse(channel));
    EEGbf = mat_signal_filtered.(chan);
    rawSignal = rawEEGSignals.(chan);

    % ------------ FILTERING ------------ %


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

    % Check main dialog cancel
    if isvalid(mainDialog) && mainDialog.CancelRequested
        if isvalid(mainDialog), close(mainDialog); end
        error('User requested cancel. Exiting analysis.');  % propagate cancellation
    end


    % ------------ POWER CALCULATIONS ------------ %

    % Calculate the absolut powers of sigma and delta filtered signals and
    % calculate the ratio for relative sigma power
    sigma_power = zeros(size(EEG_sigma));

    padL = floor((0.1/2)*Fs); % Zero padding length for the window


    for i = 1:WinS:length(EEG_sigma) - WinL + 1
        % Define window indices
        Wstart = max(1, i - padL);
        Wend   = min(length(EEG_sigma), i + WinL - 1);
        
        % Extract and pad window symmetrically if needed
        window = EEG_sigma(Wstart:Wend);
        if length(window) < WinL
            window = padarray(window,padL,'pre'); % pad at the end if short
        end
        
        % Compute power in the window
        pw = sum(window.^2) / WinL;
        
        % Assign the same value to all samples within the window
        sigma_power(i:i+WinL-1) = max(log10(pw));            
    end


    delta_power = zeros(size(EEG_delta));

    for i = 1:WinS:length(EEG_delta) - WinL + 1
        % Define window indices
        Wstart = max(1, i - padL);
        Wend   = min(length(EEG_delta), i + WinL - 1);
        
        % Extract and pad window symmetrically if needed
        window = EEG_delta(Wstart:Wend);
        if length(window) < WinL
            window = padarray(window,padL,'pre'); % pad at the end if short
        end
        
        % Compute power in the window
        pw = sum(window.^2) / WinL;
        
        % Assign the same value to all samples within the window
        delta_power(i:i+WinL-1) = max(log10(pw));            
    end

    ratio_power = sigma_power./delta_power;



    % Calculate the first IMF 
    [imf_s,~,~] = emd(EEG_sigma);
    first_imf = imf_s(:,1);


    % ------------ SPINDLE DETECTION AND EVALUATION ------------ %


    % If all thresholds are exceeded spindle is most likely detected. 
    % Detected spindles are at least 0.5 second in length.

    % Check main dialog cancel
    if isvalid(mainDialog) && mainDialog.CancelRequested
        if isvalid(mainDialog), close(mainDialog); end
        error('User requested cancel. Exiting analysis.');  % propagate cancellation
    end


    % Update the WinL and make the window overlap to be 50%
    WinL = 30*Fs;
    WinL_rms = round(Fs * 0.3);
    overlap = round(0.5*WinL);

    % Calculate the total number of windows
    total_windows = floor((length(EEGbf) - WinL) / (WinL - overlap)) + 1;

    spindles = cell(total_windows, 1);

    dd = uiprogressdlg(app.UIFigure, ...
    'Title','Spindle Analysis', ...
    'Message','Detecting spindles over night...', ...
    'Cancelable','on');

   
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
        power_window = sigma_power(startIndex:endIndex);
        % ratio_window = ratio_power(startIndex:endIndex);

        % Find indices over the threshold
        indices_over_threshold = find(power_window >= thrs_sigma);


        if isvalid(dd) && dd.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end        
      

        % Check if there are any consecutive runs within the desired length range
        if isempty(indices_over_threshold)
            % If no consecutive indices, break out of the loop and move to the next iteration
            dd.Value = i/total_windows;
            drawnow
            continue;
        end


        % Find consecutive indices
        candidate_indices = find_consecutive_spindles(indices_over_threshold,Fs);
        % Check if there are any consecutive runs within the desired length range
        if isempty(candidate_indices)
            % If no consecutive indices, break out of the loop and move to the next iteration
            dd.Value = i/total_windows;
            drawnow
            continue;
        end


        % Map indices to correspond the original signals indices
        % for each consecutive run
        filtered_candidate_indices = cell(size(candidate_indices));

        for j = 1:length(candidate_indices)
            % Map indices 
            absolute_sigma_power_indices = candidate_indices{j} + startIndex - 1;
            start_time_seconds = absolute_sigma_power_indices(1,1) / Fs;
            % Store filtered indices and starting time as a struct
            filtered_candidate_indices{j} = struct('indices', absolute_sigma_power_indices, 'start_sec', start_time_seconds);
        end

        for j = 1:length(filtered_candidate_indices)
            candidate_indices{j} = filtered_candidate_indices{j}.indices;
        end


        % CALCULATING THE CWT FROM FIRST MODE
        [cwtResult, frequencies] = cwt(first_imf(startIndex:endIndex), 'amor', Fs);

        % Find regions where CWT values exceed the threshold
        highFreqIndices = frequencies >= 11 & frequencies <= 16;
        cwtHighFreq = abs(cwtResult(highFreqIndices, :));
        % Create masks for threshold exceeding regions
        mask = cwtHighFreq > thrs_cwt;



        ripple_band = [11 16];



        % ------------ CANDIDATE EVALUATION ------------ %



        % Initialize points for the current window
        points_per_candidate = zeros(size(candidate_indices));

        % Check the amplitude
        for j = 1:size(candidate_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = candidate_indices{j}(1);
            end_index = candidate_indices{j}(end);
            % Extract the corresponding part amplitude cahnge values
            Peak2Peak = max(EEG_sigma(start_index:end_index)) - min(EEG_sigma(start_index:end_index));
            Peak2Peak_bf = max(EEGbf(start_index:end_index)) - min(EEGbf(start_index:end_index));
            % Process the result based on peak to peak values
            if any(Peak2Peak >= thrs_P2P)
                % Increment the counter for the current consecutive part
                points_per_candidate(j) = points_per_candidate(j) + 1;
            end
            % Process the result based on peak to peak values.
            if any(Peak2Peak >= 100)
                % Increment the counter for the current consecutive part
                points_per_candidate(j) = points_per_candidate(j) - 1;
            end
            if Peak2Peak_bf == 200
                % Increment the counter for the current consecutive part
                points_per_candidate(j) = points_per_candidate(j) - 5;
            end
        end



        % Check values for rms.
        for j = 1:size(candidate_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = candidate_indices{j}(1) - startIndex + 1;
            end_index = candidate_indices{j}(end) - startIndex + 1;
            % Subset rms_values array between start_index and end_index
            subset_rms_values = rms_values_bf(start_index:end_index)';
            % Find indices where the RMS values exceed the threshold
            indices_over_threshold_rms = find(subset_rms_values >= 18,1);
            % Process the result based on the rms values
            if ~isempty(indices_over_threshold_rms)
                % Increment the counter for the current consecutive part
                points_per_candidate(j) = points_per_candidate(j) - 2;
            end
        end


         % Artifact detection through FFT
         for j = 1:size(candidate_indices, 1)
            start_index = candidate_indices{j}(1);
            end_index = candidate_indices{j}(end);
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
                points_per_candidate(j) = points_per_candidate(j) - 3;
            end 

        end




        % Iterate over each consecutive run and check values in ratio_window
        for j = 1:size(candidate_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = candidate_indices{j}(1);
            end_index = candidate_indices{j}(end);
            % Extract the corresponding part from ratio_window
            ratio_consecutive_part = ratio_power(start_index:end_index);
            % Process the result based on the values in ratio_window
            if any(ratio_consecutive_part >= thrs_relative)
                % Increment the counter for the current consecutive part
                points_per_candidate(j) = points_per_candidate(j) + 1;
            end
        end



        % Iterate over each consecutive run and check rms values in
        % sigma band
        for j = 1:size(candidate_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = candidate_indices{j}(1) - startIndex + 1;
            end_index = candidate_indices{j}(end) - startIndex + 1;
            % Extract the corresponding part from ratio_window
            max_rms_sigma = max(rms_values_sigma(start_index:end_index));
            % Process the result based on the values in ratio_window
            if max_rms_sigma >= thrs_rms_sigma
                % Increment the counter for the current consecutive part
                points_per_candidate(j) = points_per_candidate(j) + 1;
            end
        end





        % Iterate over each consecutive run and check if there are
        % areas that exceeds the CWT threshold.
        for j = 1:size(candidate_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = candidate_indices{j}(1) - startIndex + 1;
            end_index = candidate_indices{j}(end) - startIndex + 1;
            timeIndices = start_index:end_index;
            % Pick up the values from cwt results for candidate area
            candidate_areas = bwconncomp(mask(:,timeIndices));
            % Initialize a variable for the number of areas in the specified time range
            numAreasInRange = candidate_areas.NumObjects;
            % Check if values in candidate area exceeds the cwt
            % threshold
            if numAreasInRange > 0
                % Process the result based on the average value
                points_per_candidate(j) = points_per_candidate(j) + 1;
            end
        end



        % Peak detection based on Hilbert transform
        for j = 1 : size(candidate_indices,1)  
            numAreas = zeros(1,2);
            start_index = candidate_indices{j}(1);
            end_index = candidate_indices{j}(end);

            extension = 2*Fs;

            safe_start = max(1, start_index - extension);
            safe_end = min(length(EEGbf), end_index + extension);

            [b,a] = butter(order, ripple_band/(Fs/2),'bandpass');

            filt_bf = filtfilt(b, a, EEGbf(safe_start:safe_end));

            env_bf = abs(hilbert(EEGbf(safe_start:safe_end)));
            env_sigma = abs(hilbert(filt_bf));
            
            baseline_mean1 = mean(env_bf);
            baseline_std1  = std(env_bf);
            thr1_global    = baseline_mean1 + 2*baseline_std1;
                
            baseline_mean2 = mean(env_sigma);
            baseline_std2  = std(env_sigma);
            thr2_global    = baseline_mean2 + thrs_envelope*baseline_std2;

            signal_part_bf = EEGbf(start_index:end_index);
            signal_part_sigma = EEG_sigma(start_index:end_index);
                      
    
            env1 = abs(hilbert(signal_part_bf));
            env2 = abs(hilbert(signal_part_sigma));

            rippleMask1 = env1 > thr1_global;  
            rippleMask2 = env2 > thr2_global;
            CCripple1 = bwconncomp(rippleMask1);
            CCripple2 = bwconncomp(rippleMask2);
                
            % Get the number of areas
            numAreas(1) = CCripple1.NumObjects;
            numAreas(2) = CCripple2.NumObjects;

            if numAreas(1) >= 1
                points_per_candidate(j) = points_per_candidate(j) + 1;
            end

            if numAreas(2) >= 1
                points_per_candidate(j) = points_per_candidate(j) + 2;
            end

         end




        % ------------ EVALUATION ------------ %

        % Check if any consecutive part has accumulated five points
        if any(points_per_candidate >= 5)
            % Find indices of consecutive parts with five points
            indices_three_points = find(points_per_candidate >= 5);
            % Create a cell array to store spindle indices for the current window
            current_spindle_indices = cell(length(indices_three_points), 1);
            % Iterate over indices and store corresponding segments in current_spindle_indices
            for k = 1:length(indices_three_points)
                current_consecutive_part_indices = filtered_candidate_indices{indices_three_points(k)};
                current_spindle_indices{k} = current_consecutive_part_indices.indices(1,1):current_consecutive_part_indices.indices(end,1);
                current_spindle_starting_time = filtered_candidate_indices{indices_three_points(k)};
                current_spindle_starting_time = current_spindle_starting_time.start_sec;
                current_spindle_duration = length(current_spindle_indices{k,1}) / Fs;
                current_spindle_amplitude = max(EEGbf(current_spindle_indices{k,1})) - min(EEGbf(current_spindle_indices{k,1}));
                
                
                current_spindle_indices{k} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', chan, 'sampling_frequency', Fs, ...
                                                        'indices', current_spindle_indices{k}, 'start_sec', current_spindle_starting_time, ...
                                                        'duration', current_spindle_duration, 'amplitude', current_spindle_amplitude);
            end
            % Store spindle indices for the current window in all_spindle_indices
            spindles{i} = current_spindle_indices;

            if isvalid(dd) && dd.CancelRequested
                if isvalid(mainDialog), close(mainDialog); end
                error('User requested cancel. Exiting analysis.');  % propagate cancellation
            end
        end

        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end

        % Update waitbar
        dd.Value = i/total_windows;
        drawnow

    end
    close(dd);

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

    startTimes = cellfun(@(s) s.start_sec, spindles);
    [~, idx] = sort(startTimes);
    spindles = spindles(idx);



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
            max_power = max(sigma_power(combined_indices));

            segment = EEGbf(combined_indices);
            [pxx, f] = pwelch(segment, [], [], [], Fs);
            
            % Limit to spindle range (9–16 Hz)
            valid = f >= 9 & f <= 16;
            f_valid = f(valid);
            pxx_valid = pxx(valid);
            
            % Identify peak frequency
            [~, idx_max] = max(pxx_valid);
            f_peak = f_valid(idx_max);
            
            connectedEvents{i} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', chan, 'sampling_frequency', Fs, ...
                                        'indices', combined_indices, 'start_sec', current_event_start, ...
                                        'duration', combined_duration, 'amplitude', new_amplitude, 'spindle_power', max_power, ...
                                        'max_freq_fft', f_peak);
            i = i + 2;
        else
            max_power = max(sigma_power(current_event_indices));
    
            segment = EEGbf(current_event_indices);
            [pxx, f] = pwelch(segment, [], [], [], Fs);
            
            % Limit to spindle range (9–16 Hz)
            valid = f >= 9 & f <= 16;
            f_valid = f(valid);
            pxx_valid = pxx(valid);
            
            % Identify peak frequency
            [~, idx_max] = max(pxx_valid);
            f_peak = f_valid(idx_max);

            connectedEvents{i} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', chan, 'sampling_frequency', Fs, ...
                                        'indices', current_event_indices, 'start_sec', current_event_start, ...
                                        'duration', current_event_dur, 'amplitude', current_event_amp, ...
                                        'spindle_power', max_power, 'max_freq_fft', f_peak);
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

end





% ------------ FUNCTIONS USED IN SPINDLE DETECTION ------------ %




function filtered_consecutive_indices = find_consecutive_spindles(indices, Fs)
    % This function identifies and filters consecutive runs of indices

    % Initialize variables
    consecutive_indices = {};
    current_run = [];

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
% -------------------- DELTA WAVE / SWS DETECTION -------------------------
%%-------------------------------------------------------------------------





function delta_waves = Detect_delta_waves(app, Fs, mat_signal_filtered, filename, save_path, rawEEGSignals, Thresholds, mainDialog)



channels_to_analyse = fieldnames(mat_signal_filtered);

% Default thresholds
thrs_power = Thresholds.Power_SWS;            % Delta band 
thrs_envelope = Thresholds.Envelope_SWS;      % Delta band 
thrs_rms = Thresholds.RMS_SWS;                % Delta band
thrs_wave_length = Thresholds.WaveLength_SWS; % Delta band 
thrs_cwt = Thresholds.CWT_SWS;                % Delta band 

% Alternative annotation name
aasmEvent_name = 'UNSURE';


% Create a table where to save all detected slow waves from different
% channels

% Define column names
% columnNames = channels_to_analyse;
% Create an empty cell array with specified column names
dataTableSWS = cell(1, numel(channels_to_analyse)); % Start with one row
% Assign column names to the first row of the cell array
dataTableSWS(1, :) = channels_to_analyse;

event = 'Delta_waves';


for channel = 1 : length(channels_to_analyse)
    chan = cell2mat(channels_to_analyse(channel));
    EEGbf = mat_signal_filtered.(chan);
    rawSignal = rawEEGSignals.(chan);


    % ------------ FILTERING AND POWER CALCULATIONS ------------ %


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
    delta_power = zeros(size(EEG_delta));
    padL = floor((1/2)*Fs); % Zero padding length for the window


    for i = 1:WinS:length(EEG_delta) - WinL + 1
        % Define window indices
        Wstart = max(1, i - padL);
        Wend   = min(length(EEG_delta), i + WinL - 1);
        
        % Extract and pad window symmetrically if needed
        window = EEG_delta(Wstart:Wend);
        if length(window) < WinL
            window = padarray(window,padL,'pre'); % pad at the end if short
        end
        
        % Compute power in the window
        pw = sum(window.^2) / WinL;
        
        % Assign the same value to all samples within the window
        delta_power(i:i+WinL-1) = max(log10(pw));            
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

    % Check main dialog cancel
    if isvalid(mainDialog) && mainDialog.CancelRequested
        if isvalid(mainDialog), close(mainDialog); end
        error('User requested cancel. Exiting analysis.');  % propagate cancellation
    end




    % ------------ DELTA WAVE DETECTION AND EVALUATION ------------ %



    % Update the WinL and make the window overlap to be 50%
    WinL = 30*Fs;
    overlap = round(0.5*WinL);

    % Calculate the total number of windows
    total_windows = floor((length(EEGbf) - WinL) / (WinL - overlap)) + 1;

    % Initialize the score table for possible candidates for later
    % evaluation.
    delta_waves = cell(total_windows, 1);

    dd = uiprogressdlg(app.UIFigure, ...
    'Title','Delta Wave Analysis', ...
    'Message','Detecting delta waves over night...', ...
    'Cancelable','on');

    ripple_band = [0.1 4.5];
    ripple_band_high = [20 30];


    for i = 1:total_windows
        % Calculate start and end indices for the current window
        startIndex = 1 + (i - 1) * (WinL - overlap);
        endIndex = startIndex + WinL - 1;
        % Ensure endIndex does not exceed the length of the signal
        endIndex = min(endIndex, length(EEGbf));


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

        if isvalid(dd) && dd.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end


        if isempty(filtered_consecutive_indices)
            dd.Value = i/total_windows;
            drawnow
            continue;
        end

        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end



        event_indices = cell(size(filtered_consecutive_indices));

        % Pick up the indices to speed up the process
        for j = 1:length(filtered_consecutive_indices)
            event_indices{j} = filtered_consecutive_indices{j,1}.indices;
        end



        % ------------ CANDIDATE EVALUATION ------------ %



        % Initialize points for the current window
        points_per_candidate = zeros(size(event_indices));


        % Check if the power exceeds the threshold in candidate area
        for j = 1:size(event_indices, 1)
            % Get start and end indices of the current consecutive run
            start_index = event_indices{j}(1);
            end_index = event_indices{j}(end);
            % Extract the corresponding part from EEG_sigma
            power_part = delta_power(start_index:end_index);
            % Process the result based on the wave length
            if any(power_part >= thrs_power)
               points_per_candidate(j) = points_per_candidate(j) + 1;
            end
        end



        % Wave length check. In this part the code checks the wave
        % length of the signal in detected area.
        for j = 1:size(event_indices, 1)
             % Get start and end indices of the current consecutive run
             start_index = event_indices{j}(1);
             end_index = event_indices{j}(end);

             % Extract the corresponding part from EEG_delta
             delta_consecutive_part = EEG_delta(start_index:end_index);

             [~, locs] = findpeaks(delta_consecutive_part);
             if length(locs) > 1
                 period_samples = diff(locs); % Differences between consecutive peaks
                 periods = period_samples / Fs; % Convert to time periods
                 % Count the number of periods greater than the threshold
                 num_long_periods = sum(periods >= thrs_wave_length);
                 wavelength = 0;
             else
                 % Calculate the wavelength extrapolating in scenarios
                 % where only one peak is detected.
                 res = wavelength_phase_with_fallback(delta_consecutive_part,Fs,'MaxAllowedPeriod', 3.0);
                 wavelength = res.wavelength;
                 num_long_periods = length(locs);
             end

             % Process the result based on the count of long periods
             if num_long_periods >= 3 || wavelength >= 0.55
                 points_per_candidate(j) = points_per_candidate(j) + 1;
             end
        end





        for j = 1:size(event_indices, 1)
            start_index = event_indices{j}(1);
            end_index = event_indices{j}(end);
            % Extract the corresponding part from EEG_sigma
            rms_part_delta = rms_values_delta(start_index:end_index);
            % rms_part_raw = rms_raw(start_index:end_index);
            if max(rms_part_delta) >= thrs_rms
                % Process the result based on the values in envelope
                points_per_candidate(j) = points_per_candidate(j) + 1;
            end
        end



        % Artefact detection
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
                points_per_candidate(j) = points_per_candidate(j) - 1;
            end 

        end




        % Peak detection based on Hilbert transform
        for j = 1 : size(event_indices,1)  
            numAreas = zeros(1,3);
            start_index = event_indices{j}(1);
            end_index = event_indices{j}(end);

            extension = 10*Fs;

            safe_start = max(1, start_index - extension);
            safe_end = min(length(EEGbf), end_index + extension);

            [b,a] = butter(order, ripple_band/(Fs/2),'bandpass');
            [b2,a2] = butter(order, ripple_band_high/(Fs/2),'bandpass');

            filt_bf = filtfilt(b, a, EEGbf(safe_start:safe_end));
            filt_bf_high = filtfilt(b2, a2, EEGbf(safe_start:safe_end));

            env_bf = abs(hilbert(EEGbf(safe_start:safe_end)));
            env_delta = abs(hilbert(filt_bf));

            env_bf_high = abs(hilbert(filt_bf_high));
               
            baseline_mean1 = mean(env_bf);
            baseline_std1  = std(env_bf);
            thr1_global    = baseline_mean1 + 2*baseline_std1;

            baseline_mean1_high = mean(env_bf_high);
            baseline_std1_high  = std(env_bf_high);
            thr1_global_high    = baseline_mean1_high + 2.6*baseline_std1_high;
                
            baseline_mean2 = mean(env_delta);
            baseline_std2  = std(env_delta);
            thr2_global    = baseline_mean2 + thrs_envelope*baseline_std2;

            signal_part_bf = EEGbf(start_index:end_index);
            signal_part_delta = EEG_delta(start_index:end_index);
                      
            filt_event_bf_high = filtfilt(b2,a2, signal_part_bf);
    
            env1 = abs(hilbert(signal_part_bf));
            env1_high = abs(hilbert(filt_event_bf_high));
            env2 = abs(hilbert(signal_part_delta));

            rippleMask1 = env1 > thr1_global;  
            rippleMask1_high = env1_high > thr1_global_high;
            rippleMask2 = env2 > thr2_global;
            CCripple1 = bwconncomp(rippleMask1);
            CCripple2 = bwconncomp(rippleMask2);
            CCripple3 = bwconncomp(rippleMask1_high);
                
            % Get the number of areas
            numAreas(1) = CCripple1.NumObjects;
            numAreas(2) = CCripple2.NumObjects;
            numAreas(3) = CCripple3.NumObjects;

            if numAreas(1) >= 2
                points_per_candidate(j) = points_per_candidate(j) + 1;
            end

            if numAreas(2) >= 2
                points_per_candidate(j) = points_per_candidate(j) + 1;
            end

            if numAreas(3) >= 2
                points_per_candidate(j) = points_per_candidate(j) - 1;
            end

         end



         % ------------ SAVING ACCEPTED CANDIDATES ------------ %

        

        % Check if any consecutive part has accumulated four points
        if any(points_per_candidate >= 4)
            % Find indices of consecutive parts with four points
            indices_four_points = find(points_per_candidate >= 4);
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

            if isvalid(dd) && dd.CancelRequested
                if isvalid(mainDialog), close(mainDialog); end
                error('User requested cancel. Exiting analysis.');  % propagate cancellation
            end
        end

        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end

        % Update waitbar
        dd.Value = i/total_windows;
        drawnow
    end

    close(dd);

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

    startTimes = cellfun(@(s) s.start_sec, delta_waves);
    [~, idx] = sort(startTimes);
    delta_waves = delta_waves(idx);



    % ------------ OVERLAP CORRECTION ------------ %


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
            maxDeltaPower = max(delta_power(combined_indices));
            sigPart = EEG_delta(combined_indices);
            N = length(sigPart);
            Y = fft(sigPart .* hann(N));
            P = abs(Y/N).^2;
            f = (0:N-1)*(Fs/N);
            
            valid = f >= 0.1 & f <=5;

            f_valid = f(valid);
            pxx_valid = P(valid);
            
            % Identify peak frequency
            [~, idx_max] = max(pxx_valid);
            f_peak = f_valid(idx_max);

            connectedEvents{i} = struct('indices', combined_indices, 'start_sec', current_event_start, ...
                                        'duration', combined_duration, 'amplitude', new_amplitude, ...
                                        'filtered_signal', EEGbf(combined_indices), 'delta_power', maxDeltaPower, ...
                                        'max_fft', f_peak);
            i = i + 2;
        else
            maxDeltaPower = max(delta_power(current_event_indices));
            sigPart = EEG_delta(current_event_indices);
            N = length(sigPart);
            Y = fft(sigPart .* hann(N));
            P = abs(Y/N).^2;
            f = (0:N-1)*(Fs/N);
            
            valid = f >= 0.1 & f <=5;

            f_valid = f(valid);
            pxx_valid = P(valid);
            
            % Identify peak frequency
            [~, idx_max] = max(pxx_valid);
            f_peak = f_valid(idx_max);
            connectedEvents{i} = struct('indices', current_event_indices, 'start_sec', current_event_start, ...
                                        'duration', current_event_dur, 'amplitude', current_event_amp, ...
                                        'filtered_signal', current_event_filt, 'delta_power', maxDeltaPower, ...
                                        'max_fft', f_peak);
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
            maxDeltaPower = max(delta_power(combined_indices));
            sigPart = EEG_delta(combined_indices);
            N = length(sigPart);
            Y = fft(sigPart .* hann(N));
            P = abs(Y/N).^2;
            f = (0:N-1)*(Fs/N);
            
            valid = f >= 0.1 & f <=5;

            f_valid = f(valid);
            pxx_valid = P(valid);
            
            % Identify peak frequency
            [~, idx_max] = max(pxx_valid);
            f_peak = f_valid(idx_max);
            connectedEvents{i} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', chan, 'sampling_frequency', Fs, ...
                                        'indices', combined_indices, 'start_sec', current_event_start, ...
                                        'duration', combined_duration, 'amplitude', new_amplitude, 'delta_power', maxDeltaPower, ...
                                        'max_fft', f_peak);
            i = i + 2;
        else
            maxDeltaPower = max(delta_power(current_event_indices));
            sigPart = EEG_delta(current_event_indices);
            N = length(sigPart);
            Y = fft(sigPart .* hann(N));
            P = abs(Y/N).^2;
            f = (0:N-1)*(Fs/N);
            
            valid = f >= 0.1 & f <=5;

            f_valid = f(valid);
            pxx_valid = P(valid);
            
            % Identify peak frequency
            [~, idx_max] = max(pxx_valid);
            f_peak = f_valid(idx_max);
            connectedEvents{i} = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'input_channel', chan, 'sampling_frequency', Fs, ...
                                        'indices', current_event_indices, 'start_sec', current_event_start, ...
                                        'duration', current_event_dur, 'amplitude', current_event_amp, 'delta_power', maxDeltaPower, ...
                                         'max_fft', f_peak);
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

end




% ------------ FUNCTIONS USED IN DELTA WAVE DETECTION------------ %


function out = wavelength_phase_with_fallback(x, Fs, varargin)
    % wavelength_phase_with_fallback  Estimate wavelength (seconds) from signal x.

    % parse inputs
    p = inputParser();
    addParameter(p, 'MaxAllowedPeriod', Inf, @(v) isnumeric(v) && isscalar(v) && v>0);
    addParameter(p, 'RefIndex', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
    addParameter(p, 'Bandpass', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
    parse(p, varargin{:});
    maxAllowed = p.Results.MaxAllowedPeriod;
    refIdx = p.Results.RefIndex;
    bp = p.Results.Bandpass;

    % prepare
    x = x(:);
    N = length(x);
    t = (0:N-1)'/Fs;

    % optional bandpass
    if ~isempty(bp)
        [b,a] = butter(2, bp./(Fs/2), 'bandpass');
        x = filtfilt(b,a,x);
    end

    % analytic signal and phase
    z = hilbert(x);
    phi = unwrap(angle(z));

    % choose reference index
    if isempty(refIdx)
        % prefer a local peak near start of window; otherwise global max
        [~, locs_p] = findpeaks(x);
        if ~isempty(locs_p)
            refIdx = locs_p(1);
        else
            [~, refIdx] = max(x);
        end
    end
    t_ref = t(refIdx);
    phi0 = phi(refIdx);

    % target phase (next same phase -> +2*pi)
    target = phi0 + 2*pi;

    % search forward from reference
    phi_forward = phi(refIdx:end);
    t_forward = t(refIdx:end);

    out = struct('wavelength', NaN, ...
                 'confidence', 0, ...
                 'method', 'none', ...
                 'lower_bound', NaN, ...
                 't_ref', t_ref, ...
                 't_cross', NaN);

    % If we already have a full cycle inside the window:
    idx_full = find(phi_forward >= target, 1, 'first');
    if ~isempty(idx_full) && idx_full > 1
        % linear interpolate between idx_full-1 and idx_full
        i1 = idx_full-1; i2 = idx_full;
        phi1 = phi_forward(i1); phi2 = phi_forward(i2);
        t1 = t_forward(i1); t2 = t_forward(i2);
        if phi2 == phi1
            t_cross = t2;
        else
            frac = (target - phi1) / (phi2 - phi1);
            t_cross = t1 + frac*(t2 - t1);
        end
        out.t_cross = t_cross;
        out.wavelength = t_cross - t_ref;
        out.confidence = 1;
        out.method = 'observed';
        out.lower_bound = out.wavelength; % full cycle observed -> lower bound == estimate
        return;
    end

    % No full cycle observed -> compute guaranteed lower bound:
    t_end = t(end);
    lower_bound = (t_end - t_ref);         % time remaining in window after reference
    out.lower_bound = lower_bound;

    % Try linear-extrapolation of phase slope to estimate crossing time
    % Fit straight line to phi_forward vs t_forward
    dphi = diff(phi_forward);
    dt = 1/Fs;
    % use median derivative for robustness
    slope = median(dphi)/dt;   % rad/s

    % If slope is near zero or negative -> cannot extrapolate reasonably
    if slope <= 0 || ~isfinite(slope)
        out.method = 'lower_bound_only';
        out.wavelength = NaN;    % or could set to lower_bound if you want a numeric fallback
        out.confidence = 0;
        return;
    end

    % time required for phase to increase by 2*pi
    time_needed = (2*pi) / slope;
    t_cross_est = t_ref + time_needed;

    % compute estimated wavelength
    wavelength_est = time_needed;

    % If extrapolated crossing lies extremely far beyond allowed max, clamp and set low confidence
    if wavelength_est > maxAllowed
        % clamp to maxAllowed (user decision) and mark low confidence
        out.t_cross = t_ref + maxAllowed;
        out.wavelength = maxAllowed;
        out.confidence = 0;
        out.method = 'extrapolated_clamped';
        return;
    end

    % Accept extrapolated estimate
    out.t_cross = t_cross_est;
    out.wavelength = wavelength_est;
    out.confidence = 0;        % extrapolated -> lower confidence than observed
    out.method = 'extrapolated';
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
% --------------------------- CAP DETECTION -------------------------------
%%-------------------------------------------------------------------------







function cap = Detect_CAP(app, Fs, mat_signal_filtered, filename, save_path, Thresholds, mainDialog)

channels_to_analyse = fieldnames(mat_signal_filtered);

dd = uiprogressdlg(app.UIFigure, ...
    'Title','CAP Analysis', ...
    'Message','Detecting CAPs...', ...
    'Cancelable','on');

% Default thresholds
thrs_power = Thresholds.Power_CAP;         % Filtered signal (base frequency)
thrs_envelope = Thresholds.Envelope_CAP;
thrs_descriptor = Thresholds.Desc_CAP;     % For all band frequencies
thrs_cwt = Thresholds.CWT_CAP;             % Filtered signal (base frequency)
amp_diff = Thresholds.AmpDiff_CAP;         % Filtered signal (base frequency)

% Alternative annotation name
aasmEvent_name = 'UNSURE';



% ------------ FILTERING ------------ %


bands = {'Delta', 'Theta', 'Alpha', 'Sigma', 'Beta'};
band_freq = [0.1, 4; 4, 8; 8, 12; 12, 16; 16, 30]; % Frequency ranges for each band

% Initialize cell array to store filtered signals
filtered_signals = cell(numel(channels_to_analyse), numel(bands));

for i = 1 : numel(channels_to_analyse)
    channel = cell2mat(channels_to_analyse(i));
    signal = mat_signal_filtered.(channel);
    % Filter the signals for each frequency band
    for j = 1:numel(bands)
        % Design the bandpass filter using the Butterworth filter
        [b, a] = butter(2, band_freq(j, :)/(Fs/2), 'bandpass');
        % Apply the bandpass filter to the signal
        filtered_signals{i, j} = filter(b, a, signal);
    end

    if i == numel(channels_to_analyse)
        EEGbf = mat_signal_filtered.(channel);
    end
end





% ------------ RMS CALCULATION ------------ %

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
    % Assign the RMS value to the corresponding positions in rms_values
    rms_bf(start_idx:end_idx) = rms_value_bf;
end

dd.Value = 0.2;
drawnow;


% ------------ POWER CALCULATIONS ------------ %


% Create a table with named columns
Signal_info = table('Size', [1, numel(channels_to_analyse) + 1], ...
                    'VariableTypes', [{'cell'}, repmat({'cell'}, 1, numel(channels_to_analyse))], ...
                    'VariableNames', [{'Type'}, channels_to_analyse(:).']);

% Assign the signal band names to the first row and channel names to the second row
Signal_info{1, 'Type'} = {'Signal power'};


powerWinL = floor(2*Fs);
overlap = round(0.5*powerWinL);

for i = 1:numel(channels_to_analyse)
        channel = cell2mat(channels_to_analyse(i));
        signal = mat_signal_filtered.(channel);

        total_windows = floor((length(signal) - powerWinL) / (powerWinL - overlap)) + 1;

        % Initialize arrays to hold power and contribution count
        power = zeros(size(signal));
        contribution_count = zeros(size(signal));


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

        end

        % Normalize the power by the contribution count
        power = power ./ max(contribution_count, 1);

        Signal_info.(channel){1} = power;

end

% Check main dialog cancel
if isvalid(mainDialog) && mainDialog.CancelRequested
    if isvalid(mainDialog), close(mainDialog); end
    error('User requested cancel. Exiting analysis.');  % propagate cancellation
end

dd.Value = 0.3;
drawnow;

        if isvalid(dd) && dd.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end


% ------------ SHORT-LONG DESCRIPTOR CALCULATIONS ------------ %

    descriptors = cell(numel(channels_to_analyse),numel(bands));


    for n = 1 : numel(channels_to_analyse)
        [imf_delta,~,~] = emd(filtered_signals{n,1});
        [imf_theta,~,~] = emd(filtered_signals{n,2});
        [imf_alpha,~,~] = emd(filtered_signals{n,3});
        [imf_sigma,~,~] = emd(filtered_signals{n,4});
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

            for k = 1:total_windows
                startIndex = 1 + (k - 1) * (WinL_long - overlap_long);
                endIndex = startIndex + WinL_long - 1;

                % Ensure endIndex does not exceed the length of the signal
                endIndex = min(endIndex, length(imf_signal));

                segment = imf_signal(startIndex:endIndex);
                mean_segment = mean(abs(segment));
                long_average(startIndex:endIndex) = mean_segment;
            end


            mean_amplitudes_long{j} = long_average;

            % Total windows for shorter window average
            total_windows = floor((length(imf_signal) - WinL_short) / (WinL_short - overlap_short)) + 1;


            for k = 1:total_windows
                startIndex = 1 + (k - 1) * (WinL_short - overlap_short);
                endIndex = startIndex + WinL_short - 1;

                % Ensure endIndex does not exceed the length of the signal
                endIndex = min(endIndex, length(imf_signal));

                segment = imf_signal(startIndex:endIndex);
                mean_segment = mean(abs(segment));
                short_average(startIndex:endIndex) = mean_segment;
            end

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
            end
        end
    end

        if isvalid(dd) && dd.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end

    % ------------ BINARY IMAGE RECOGNITION ------------ %


dd.Value = 0.5;
drawnow;


    % Parameters
    epochL = 60 * Fs; % Length of epoch in samples
    epoch_overlap = round(0.5 * epochL); % Overlap in samples
    epochs = floor((length(EEGbf) - epochL) / (epochL - epoch_overlap)) + 1;

    merge_distance = 2 * Fs;
    min_region_length = 2 * Fs;

    A_phase_candidates = cell(0,numel(channels_to_analyse));

    for i = 1:epochs
        startIndex = 1 + (i - 1) * (epochL - epoch_overlap);
        endIndex = startIndex + epochL - 1;

        % Ensure endIndex does not exceed the length of the signal
        endIndex = min(endIndex, length(EEGbf));

        all_binary_images = zeros(length(descriptors), epochL);


        if size(descriptors, 1) == 1
            descriptors(2,:) = descriptors(1,:);
        end


        % Iterate over each descriptor result to plot masks
        for j = 1:length(descriptors)
            descriptor_1 = descriptors{1,j};
            descriptor_2 = descriptors{2,j};

            % Find regions where descriptors exceed the threshold
            exceed_1 = descriptor_1(startIndex:endIndex) >= thrs_descriptor;
            exceed_2 = descriptor_2(startIndex:endIndex) >= thrs_descriptor;

            % Find the intersection of exceedance regions
            combined_exceed = exceed_1 & exceed_2;

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
                all_binary_images(j,:) = merged_binary_image;
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

        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end
        if isvalid(dd) && dd.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end



        % Combine all binary images into a single binary image
        combined_binary_image = any(all_binary_images, 1);


% ------------ SUB SET EVALUATION (A1, A2, A3) ------------ %



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
    end



dd.Value = 0.7;
drawnow;



% ------------ A-PHASE DETECTION AND EVALUATION ------------ %


    % Descriptors show the local amplitude change that points out from the
    % background EEG. Next thing is to find the indices for each cahnnel where the
    % indices exceeds the thresholds. Both channels must exceed the
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



    % ------------ A-PHASE EVALUATION ------------ %

    % Check main dialog cancel
    if isvalid(mainDialog) && mainDialog.CancelRequested
        if isvalid(mainDialog), close(mainDialog); end
        error('User requested cancel. Exiting analysis.');  % propagate cancellation
    end


    points_per_candidate = zeros(size(A_phase_candidates(:,1)));
    
    if ~size(channels_to_analyse,1) == 1
        first_signal_power = Signal_info.(channels_to_analyse{1}){1,1};
        second_signal_power = Signal_info.(channels_to_analyse{2}){1,1};
    else
        first_signal_power = Signal_info.(channels_to_analyse{1}){1,1};
        second_signal_power = Signal_info.(channels_to_analyse{1}){1,1};
    end


    ripple_band = [20 30];
    order = 2;

    for i = 1:size(A_phase_candidates, 1)

        % Pick up the indices for current A phase candidate
        indices = A_phase_candidates{i,1};
        phase_start = indices(1);
        phase_end = indices(end);


        % Preallocate space for storing number of areas and other results
        numAreasCWT = zeros(1, 2);  % For storing number of areas for both signals
        isEqualLength = zeros(1, 2);  % For checking the equal length condition

        % Check the candidates with CWT
        for j = 1:numel(channels_to_analyse)
            channel = cell2mat(channels_to_analyse(j));
            % Extract the current signal
            current_signal = mat_signal_filtered.(channel);

            % Perform CWT
            [cwtResult, frequencies] = cwt(current_signal(phase_start:phase_end), 'amor', Fs);

            highFreqIndices = frequencies >= 0.1 & frequencies <= 30;

            cwtHighFreq = abs(cwtResult(highFreqIndices, :));

            % Create masks for threshold exceeding regions
            mask = cwtHighFreq >= thrs_cwt;

            % Find connected components in the masks
            CC = bwconncomp(mask);

            % Get the number of areas
            numAreasCWT(j) = CC.NumObjects;

            % Get the plot length (width)
            area_length = size(mask, 2);  % Length of the plot

            % Get the properties of the connected components
            props = regionprops(CC, 'BoundingBox');

             % Extract the widths (lengths) of the bounding boxes
            widths = arrayfun(@(x) x.BoundingBox(3), props);

            % Check if any bounding box width equals the plot length
            isEqualLength(j) = any(widths == area_length);
        end

        if isscalar(channels_to_analyse)
            numAreasCWT(2) = numAreasCWT(1);
        end

        % Now compare results for both signals
        if numAreasCWT(1) >= 2 && numAreasCWT(2) >= 2 || isEqualLength(1) >= 1 && numAreasCWT(2) >= 2 || isEqualLength(2) >= 1 && numAreasCWT(1) >= 2
            points_per_candidate(i) = points_per_candidate(i) + 1;
        end



        numAreas = zeros(1, 2);
        
        % Peak detection based on Hilbert transform
        for j = 1:numel(channels_to_analyse)               
            channel = cell2mat(channels_to_analyse(j));
            % Extract the current signal
            current_signal = mat_signal_filtered.(channel);

            Aphase = current_signal(phase_start:phase_end);
    
            [b,a] = butter(order, ripple_band/(Fs/2),'bandpass');
            filt_sig = filtfilt(b,a, Aphase);
            
            if j == 1
                env1 = abs(hilbert(filt_sig));
                thr1 = mean(env1) + thrs_envelope*std(env1);   % constant threshold
                rippleMask1 = env1 > thr1;  
                CCripple1 = bwconncomp(rippleMask1);
                % Get the number of areas
                numAreas(j) = CCripple1.NumObjects;
            else
                env2 = abs(hilbert(filt_sig));
                thr2 = mean(env2) + thrs_envelope*std(env2);   % constant threshold
                rippleMask2 = env2 > thr2;  
                CCripple2 = bwconncomp(rippleMask2);
                % Get the number of areas
                numAreas(j) = CCripple2.NumObjects;
            end
            
        end


        
        

        if isscalar(channels_to_analyse)
            numAreas(2) = numAreas(1);
        end

        if numAreas(1) >= 1
            points_per_candidate (i) = points_per_candidate(i) + 1;
        end

        if numAreas(2) >= 1
            points_per_candidate (i) = points_per_candidate(i) + 1;
        end



        % Check if the emd threshold is exceeded in both channels
        rms_values = rms_bf(phase_start:phase_end);

        if any(rms_values > 100)
            points_per_candidate(i) = points_per_candidate(i) - 1;
        end



        % Check if the amplitude difference in the candidate area exceeds
        % the threshold in both channels

        % Preallocate space for storing the amplitude differences
        amplitude_differences = zeros(1,2);

        for j = 1 : numel(channels_to_analyse)
            channel = cell2mat(channels_to_analyse(j));
            % Extract the current signal
            current_signal = mat_signal_filtered.(channel);

            % Calculate the amplitude difference for the current candidate
            amplitude_difference = max(current_signal(phase_start:phase_end)) - min(current_signal(phase_start:phase_end));

            % Add the amplitude difference to the preallocated variable
            amplitude_differences(j) = amplitude_difference;
        end

        if isscalar(channels_to_analyse)
            amplitude_differences(2) = amplitude_differences(1);
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

        % Check main dialog cancel
        if isvalid(mainDialog) && mainDialog.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end
        if isvalid(dd) && dd.CancelRequested
            if isvalid(mainDialog), close(mainDialog); end
            error('User requested cancel. Exiting analysis.');  % propagate cancellation
        end
    end


dd.Value = 0.9;
drawnow;


    % Find all the A phase candidates that has 3 points or more and add them  to
    % the A_phases
    A_phases = {};
    for i = 1:numel(points_per_candidate)
        if points_per_candidate(i) >= 3
            A_phases{end+1,1} = A_phase_candidates{i,1};
            A_phases{end, 2} = A_phase_candidates{i,2};
        end
    end


    % In case there are no identified A phases during analysis
    if isempty(A_phases)
        sequences = {};
        cycles = {};
        cap = struct('A_phases', {A_phases},'Sequences', {sequences}, 'Cycles', {cycles});

        % Save
        save(fullfile(save_path,[filename '_' 'CAP']), 'cap');
        return;
    end



% ------------ CYCLE AND SEQUENCE FORMATION ------------ %


    % Extract the first element from each cell in the first column
    first_elements = cellfun(@(x) x(1), A_phases(:,1));

    % Sort the rows of A_phases based on the extracted first elements
    [~, sorted_idx] = sort(first_elements);

    % Apply the sorting index to A_phases
    A_phases = A_phases(sorted_idx, :);



    difference = 2*Fs;

    [rowCount, ~] = size(A_phases);


    j = 1;

    while j <= rowCount - 1
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
    % as B phase. Definition of one CAP-cycle is: A B

    cycle_limit = 60 * Fs; % Define the threshold

    cycles = {}; % Initialize the cell array to hold all cycles

    currentCycle = {};

    [rowCount, ~] = size(A_phases);


    for i = 1:rowCount - 1
        % Get the current and next cells
        currentCell = A_phases{i,1};
        currentSubtype = A_phases{i,2};
        nextCell = A_phases{i+1,1};

        % Get the last element of the current cell and the first element of the next cell
        lastElementCurrent = currentCell(end);
        firstElementNext = nextCell(1);


        if (firstElementNext - lastElementCurrent) > 0 && ...
            (firstElementNext - lastElementCurrent) <= cycle_limit

            B_phase =(lastElementCurrent+1:firstElementNext-1)';

            currentCycle(1,end+1) = {currentCell};
            currentCycle{2,end} = 'A';
            currentCycle{3,end} = currentSubtype;

            currentCycle{1,end+1} = B_phase;
            currentCycle{2,end} = 'B';

        else
            currentCycle(1,end+1) = {currentCell};
            currentCycle{2,end} = 'A';
            currentCycle{3,end} = currentSubtype;

            cycles{end+1} = currentCycle';
            currentCycle = {};

        end

    end



    % Remove all those cells from cycles that are too short to be a cycle.
    % In other words, remove the isolated A-phases that are non-CAP.
    for i = length(cycles):-1:1
        current_cycle = cycles{i};

        if size(current_cycle,1) == 1
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

        if size(current_cycle,1) < 5
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



dd.Value = 0.95;
drawnow;
    

    % Create a data table for CAPs
    event = 'A_phase';

    for i = 1:rowCount
        current_event = A_phases{i,1};

        subtype = A_phases{i,2};
        starting_time = current_event(1)/Fs;
        duration = length(current_event)/Fs;

        if size(channels_to_analyse, 1) == 1
            channels_to_analyse{2} = channels_to_analyse{1};
        end

        current_struct = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'sampling_frequency', Fs, ...
                                'subtype', subtype, 'start_sec', starting_time, 'duration', duration, ...
                                'indices', current_event, 'channels', {channels_to_analyse{1};channels_to_analyse{2}});

        A_phases{i} = current_struct;

    end

    event = 'sequence';

    for i = 1:length(sequences)
        current_event = sequences(i);

        starting_time = current_event{1,1}{1,1}(1)/Fs;
        ending_time = current_event{1,1}{end,1}(end)/Fs;
        duration = ending_time - starting_time;

        if size(channels_to_analyse, 1) == 1
            channels_to_analyse{2} = channels_to_analyse{1};
        end

        current_struct = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'sampling_frequency', Fs, ...
                                'start_sec', starting_time, 'duration', duration, ...
                                'Sequence', sequences(i), 'channels', {channels_to_analyse{1};channels_to_analyse{2}});

        sequences{i} = current_struct;

    end


    event = 'cycle';

    for i = 1:length(cycles)
        current_event = cycles(i);

        starting_time = current_event{1,1}{1,1}(1)/Fs;
        ending_time = current_event{1,1}{2,1}(end)/Fs;
        duration = ending_time - starting_time;

        if size(channels_to_analyse, 1) == 1
            channels_to_analyse{2} = channels_to_analyse{1};
        end

        current_struct = struct('File_ID', filename, 'name', aasmEvent_name, 'event_name', event, 'sampling_frequency', Fs, ...
                                'start_sec', starting_time, 'duration', duration, ...
                                'cycle', cycles(i), 'channels', {channels_to_analyse{1};channels_to_analyse{2}});

        cycles{i} = current_struct;

    end





    cap = struct('A_phases', {A_phases},'Sequences', {sequences}, 'Cycles', {cycles});

    % Save
    save(fullfile(save_path,[filename '_' 'CAP']), 'cap');
    dd.Value = 1.0;
    drawnow;
    close(dd);

end

