function RAVEN_UI
    % Create the figure
    fig = uifigure('Name', 'RAVEN UI', 'Position', [200 200 600 600]);


    % Default threshold values
    Thresholds.WaveLength_Kcomplex = 0.35;
    Thresholds.Hilbert_Kcomplex = 60;
    Thresholds.RMS_Kcomplpex = 25;
    Thresholds.Power_Kcomplex = 2.3;
    Thresholds.CWT_Kcomplex = 30;
    Thresholds.P2P_spindle = 10;
    Thresholds.RelativePower_spindle = 0.8;
    Thresholds.Hilbert_spindle = 5;
    Thresholds.RMS_spindle = 6;
    Thresholds.Power_spindle = 1.2;
    Thresholds.CWT_spindle = 12;
    Thresholds.WaveLength_SWS = 0.35;
    Thresholds.Hilbert_SWS = 80;
    Thresholds.Power_SWS = 1.5;
    Thresholds.RMS_SWS = 50;
    Thresholds.CWT_SWS = 10;
    Thresholds.Desc_CAP = 0.18;
    Thresholds.Power_CAP = 2.5;
    Thresholds.AmpDiff_CAP = 120;
    Thresholds.CWT_CAP = 18;
    setappdata(fig,'Thresholds', Thresholds);

    % Sampling frequency for Spindles, K-complexes, and Delta waves
    uilabel(fig, 'Position', [10 570 300 22], 'Text', 'DEFINE SAMPLING FREQUENCIES FOR ANALYSIS', 'FontWeight', 'bold');
    uilabel(fig, 'Position', [10 550 305 22], 'Text', 'Sampling Frequency (Spindle, K-complex, Delta waves):');
    Fsmicro = uieditfield(fig, 'numeric', 'Position', [315 550 100 22], 'Value', 512);

    % Sampling frequency for CAP
    uilabel(fig, 'Position', [10 520 300 22], 'Text', 'Sampling Frequency (CAP cycle):');
    FsCAP = uieditfield(fig, 'numeric', 'Position', [315 520 100 22], 'Value', 128);

    % Input information for EDF path
    uilabel(fig, 'Position', [10 470 300 22], 'Text', 'DEFINE FILE FORMAT, PATH, AND MODE', 'FontWeight', 'bold');
    uilabel(fig, 'Position', [10 440 200 22], 'Text', 'Choose The File Format:');
    fileFormat = uidropdown(fig, 'Position', [150 440 100 22], ...
        'Items', {'EDF', 'SLF'}, 'Value', 'EDF');

    uilabel(fig, 'Position', [10 410 200 22], 'Text', 'Directory Path For File(s):');
    inputPathField = uieditfield(fig, 'text', 'Position', [150 410 300 22]);

    browserButtonInput = uibutton(fig, 'push', ...
        'Text', 'Browse', ...
        'Position', [455 410 55 22],...
        'ButtonPushedFcn', @(btn, event) openDirectoryExplorer(inputPathField),...
        'FontWeight','bold');

    

    % File Name
    uilabel(fig, 'Position', [10 380 100 22], 'Text', 'Filename:');
    uilabel(fig, 'Position', [305 380 250 22], 'Text', '(Fill only when using "Individual Focus" mode.)');
    uilabel(fig, 'Position', [10 350 150 22], 'Text', 'Select Analysis Mode:');
    edfName = uieditfield(fig, 'text', 'Position', [150 380 150 22], 'Placeholder','.edf / ID number (SLF)');
    mode = uidropdown(fig, 'Position', [150 350 150 22], ...
        'Items', {'Individual Focus', 'Group Analysis'}, 'Value', 'Group Analysis');

    % Channels of interest and reference subtraction
    uilabel(fig, 'Position', [10 300 300 22], 'Text', 'DEFINE CHANNELS OF INTEREST', 'FontWeight', 'bold');
    uilabel(fig, 'Position', [10 270 300 22], 'Text', 'Channels (Spindle, K-complex, Delta waves):');
    uilabel(fig, 'Position', [10 240 300 22], 'Text', 'Reference Channel Subtraction (optional):');
    uilabel(fig, 'Position', [10 210 300 22], 'Text', 'Channels (CAP cycle, two channels):');
    uilabel(fig, 'Position', [10 180 300 22], 'Text', 'Reference Channel Subtraction (optional):');
    
    channelField = uitextarea(fig, 'Position', [260 270 300 22],'Placeholder','Example: C3,M1,...');
    RefChannelField = uitextarea(fig, 'Position', [260 240 300 22],'Placeholder','Example: C3-M1,...');
    channelFieldCAP = uitextarea(fig, 'Position', [260 210 300 22],'Placeholder','Example: C3,M1,C4,M2');
    RefChannelFieldCAP = uitextarea(fig, 'Position', [260 180 300 22],'Placeholder','Example: C3-M1,C4-M2');


    % Define microstructures to be analysed
    uilabel(fig, 'Position', [10 120 300 22], 'Text', 'SELECT MICROSTRUCTURES FOR ANALYSIS', 'FontWeight', 'bold');
    uilabel(fig, 'Position', [10 90 100 22], 'Text', 'Microtructures:');
    cb_spindle = uicheckbox(fig, 'Position', [100 90 150 22], 'Text', 'Spindle');
    cb_kcomplex = uicheckbox(fig, 'Position', [165 90 150 22], 'Text', 'K-Complex');
    cb_delta_wave = uicheckbox(fig, 'Position', [250 90 150 22], 'Text', 'Delta Wave');
    cb_cap = uicheckbox(fig, 'Position', [335 90 150 22], 'Text', 'CAP');

    % Set appdata to store selected microstructures
    setappdata(fig, 'microstructures', {}); % Initialize with an empty cell array
    

    % Add a button to update the selected microstructures
    btn_micro = uibutton(fig, 'Position', [390 90 60 22],'FontWeight', 'bold', 'Text', 'Submit', ...
            'ButtonPushedFcn', @(btn_micro, event) updateMicrostructures(btn_micro,fig, cb_spindle, cb_kcomplex, cb_delta_wave, cb_cap));


    % Output Path
    uilabel(fig, 'Position', [10 30 100 22], 'Text', 'Save Path:');
    outputPathField = uieditfield(fig, 'text', 'Position', [80 30 300 22]);

     browserButtonSave = uibutton(fig, 'push', ...
        'Text', 'Browse', ...
        'Position', [390 30 55 22],...
        'ButtonPushedFcn', @(btn, event) openDirectoryExplorer(outputPathField),...
        'FontWeight','bold');

    % Run Button
    uibutton(fig, 'Position', [500 35 60 60], 'Text', 'Run', ...
        'FontWeight','bold',...
        'BackgroundColor','g',...
        'ButtonPushedFcn', @(btn, event) runRAVEN(btn_micro,Fsmicro, FsCAP, ...
        inputPathField, edfName, channelField, RefChannelField, ...
        channelFieldCAP, RefChannelFieldCAP, ...
        getappdata(fig, 'microstructures'), outputPathField, mode,fileFormat, getappdata(fig,'Thresholds')));

        % Create button for the thresholds
    uibutton(fig,'Position', [460 525 100 40],'Text','Thresholds',...
        'FontWeight','bold',...
        'ButtonPushedFcn',@(btn, event) openNewFigure(fig));
end

function openDirectoryExplorer(inputPathField)
    % Open directory explorer
    directory = uigetdir(pwd, 'Select a Directory');
    if directory ~= 0
        % Update the text field with the selected directory path
        inputPathField.Value = directory;
    end
end




function openNewFigure(fig)
    % Create the new figure
    newFig = uifigure('Name', 'Adjustable Thresholds', 'Position', [250, 250, 500, 600]);

    uilabel(newFig, 'Text', 'K-COMPLEX DETECTION', 'Position', [10, 570, 200, 20],...
        'FontWeight','bold');
    uilabel(newFig,"Text",'Wavelength: ','Position',[10, 550, 200, 20]);
    uilabel(newFig,"Text",'(seconds)','Position',[170, 550, 200, 20]);
    uilabel(newFig,"Text",'Hilbert transform: ','Position',[10, 530, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[170, 530, 200, 20]);
    uilabel(newFig,"Text",'Root mean square: ','Position',[10, 510, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[170, 510, 200, 20]);
    uilabel(newFig,"Text",'Delta band power: ','Position',[10, 490, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V" + char(178) + "/Hz)",'Position',[170, 490, 200, 20]);
    uilabel(newFig,"Text",'Continuous wavelet transform: ','Position',[10, 470, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[230, 470, 200, 20]);
    WaveLength_Kcomplex = uieditfield(newFig,"numeric",'Position',[115, 550, 50, 18],'Value',0.35);
    Hilbert_Kcomplex = uieditfield(newFig,"numeric",'Position',[115, 530, 50, 18],'Value',60);
    RMS_Kcomplpex = uieditfield(newFig,"numeric",'Position',[115, 510, 50, 18],'Value',25);
    Power_Kcomplex = uieditfield(newFig,"numeric",'Position',[115, 490, 50, 18],'Value',2.3);
    CWT_Kcomplex = uieditfield(newFig,"numeric",'Position',[175, 470, 50, 18],'Value',30);


    uilabel(newFig, 'Text', 'SPINDLE DETECTION', 'Position', [10, 440, 200, 20],...
        'FontWeight','bold');
    uilabel(newFig,"Text",'Peak to peak: ','Position',[10, 420, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[170, 420, 200, 20]);
    uilabel(newFig,"Text",'Relative power: ','Position',[10, 400, 200, 20]);
    uilabel(newFig,"Text",'Hilbert transform: ','Position',[10, 380, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[170, 380, 200, 20]);
    uilabel(newFig,"Text",'Root mean square: ','Position',[10, 360, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[170, 360, 200, 20]);
    uilabel(newFig,"Text",'Sigma band power: ','Position',[10, 340, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V" + char(178) + "/Hz)",'Position',[170, 340, 200, 20]);
    uilabel(newFig,"Text",'Continuous wavelet transform: ','Position',[10, 320, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[230, 320, 200, 20]);
    P2P_spindle = uieditfield(newFig,"numeric",'Position',[115, 420, 50, 18],'Value',10);
    RelativePower_spindle = uieditfield(newFig,"numeric",'Position',[115, 400, 50, 18],'Value',0.8);
    Hilbert_spindle = uieditfield(newFig,"numeric",'Position',[115, 380, 50, 18],'Value',5);
    RMS_spindle = uieditfield(newFig,"numeric",'Position',[115, 360, 50, 18],'Value',6);
    Power_spindle = uieditfield(newFig,"numeric",'Position',[115, 340, 50, 18],'Value',1.2);
    CWT_spindle = uieditfield(newFig,"numeric",'Position',[175, 320, 50, 18],'Value',12);



    uilabel(newFig, 'Text', 'DELTA WAVE DETECTION', 'Position', [10, 290, 200, 20],...
        'FontWeight','bold');
    uilabel(newFig,"Text",'Wavelength: ','Position',[10, 270, 200, 20]);
    uilabel(newFig,"Text",'(seconds)','Position',[170, 270, 200, 20]);
    uilabel(newFig,"Text",'Hilbert transform: ','Position',[10, 250, 200, 20]);
    uilabel(newFig,"Text", "(" + char(181) + "V)",'Position',[170, 250, 200, 20]);
    uilabel(newFig,"Text",'Delta band power: ','Position',[10, 230, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V" + char(178) + "/Hz)",'Position',[170, 230, 200, 20]);
    uilabel(newFig,"Text",'Root mean square: ','Position',[10, 210, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[170, 210, 200, 20]);
    uilabel(newFig,"Text",'Continuous wavelet transform: ','Position',[10, 190, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[230, 190, 200, 20]);
    WaveLength_SWS = uieditfield(newFig,"numeric",'Position',[115, 270, 50, 18],'Value',0.35);
    Hilbert_SWS = uieditfield(newFig,"numeric",'Position',[115, 250, 50, 18],'Value',80);
    Power_SWS = uieditfield(newFig,"numeric",'Position',[115, 230, 50, 18],'Value',1.5);
    RMS_SWS = uieditfield(newFig,"numeric",'Position',[115, 210, 50, 18],'Value',50);
    CWT_SWS = uieditfield(newFig,"numeric",'Position',[175, 190, 50, 18],'Value',10);


    uilabel(newFig, 'Text', 'CAP DETECTION', 'Position', [10, 160, 200, 20],...
        'FontWeight','bold');
    uilabel(newFig,"Text",'Short-Long term ratio: ','Position',[10, 100, 200, 20]);
    uilabel(newFig,"Text",'Signal power: ','Position',[10, 140, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V" + char(178) + "/Hz)",'Position',[145, 140, 200, 20]);
    uilabel(newFig,"Text",'Amplitude difference: ','Position',[10, 120, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[180, 120, 200, 20]);
    uilabel(newFig,"Text",'Continuous wavelet transform: ','Position',[10, 80, 200, 20]);
    uilabel(newFig,"Text","(" + char(181) + "V)",'Position',[230, 80, 200, 20]);
    Desc_CAP = uieditfield(newFig,"numeric",'Position',[130, 100, 50, 18],'Value',0.18);
    Power_CAP = uieditfield(newFig,"numeric",'Position',[90, 140, 50, 18],'Value',2.5);
    AmpDiff_CAP = uieditfield(newFig,"numeric",'Position',[125, 120, 50, 18],'Value',120);
    CWT_CAP = uieditfield(newFig,"numeric",'Position',[175, 80, 50, 18],'Value',18);

    uibutton(newFig,'Text','Apply','FontWeight','bold',...
        'Position',[350, 70, 100, 40],...
        'ButtonPushedFcn',@(btn, event) submitValues(fig,WaveLength_Kcomplex,Hilbert_Kcomplex,RMS_Kcomplpex,Power_Kcomplex,CWT_Kcomplex,...
        P2P_spindle,RelativePower_spindle,Hilbert_spindle,RMS_spindle,Power_spindle,CWT_spindle,...
        WaveLength_SWS,Hilbert_SWS,Power_SWS,RMS_SWS,CWT_SWS,...
        Desc_CAP,Power_CAP,AmpDiff_CAP,CWT_CAP, newFig));



    
end



function submitValues(fig,WaveLength_Kcomplex,Hilbert_Kcomplex,RMS_Kcomplpex,Power_Kcomplex,CWT_Kcomplex,...
        P2P_spindle,RelativePower_spindle,Hilbert_spindle,RMS_spindle,Power_spindle,CWT_spindle,...
        WaveLength_SWS,Hilbert_SWS,Power_SWS,RMS_SWS,CWT_SWS,...
        Desc_CAP,Power_CAP,AmpDiff_CAP,CWT_CAP, newFig)
    Thresholds.WaveLength_Kcomplex = WaveLength_Kcomplex.Value;
    Thresholds.Hilbert_Kcomplex = Hilbert_Kcomplex.Value;
    Thresholds.RMS_Kcomplpex = RMS_Kcomplpex.Value;
    Thresholds.Power_Kcomplex = Power_Kcomplex.Value;
    Thresholds.CWT_Kcomplex = CWT_Kcomplex.Value;
    Thresholds.P2P_spindle = P2P_spindle.Value;
    Thresholds.RelativePower_spindle = RelativePower_spindle.Value;
    Thresholds.Hilbert_spindle = Hilbert_spindle.Value;
    Thresholds.RMS_spindle = RMS_spindle.Value;
    Thresholds.Power_spindle = Power_spindle.Value;
    Thresholds.CWT_spindle = CWT_spindle.Value;
    Thresholds.WaveLength_SWS = WaveLength_SWS.Value;
    Thresholds.Hilbert_SWS = Hilbert_SWS.Value;
    Thresholds.Power_SWS = Power_SWS.Value;
    Thresholds.RMS_SWS = RMS_SWS.Value;
    Thresholds.CWT_SWS = CWT_SWS.Value;
    Thresholds.Desc_CAP = Desc_CAP.Value;
    Thresholds.Power_CAP = Power_CAP.Value;
    Thresholds.AmpDiff_CAP = AmpDiff_CAP.Value;
    Thresholds.CWT_CAP = CWT_CAP.Value;

    setappdata(fig, 'Thresholds', Thresholds);
    
    close(newFig);
end


function updateMicrostructures(btn,fig, cb1, cb2, cb3, cb4)
    % Clear the previous selections
    microstructures = {};

    % Check each checkbox and update the microstructures list
    if cb1.Value
        microstructures{end+1} = 'spindle';
    end
    if cb2.Value
        microstructures{end+1} = 'kcomplex';
    end
    if cb3.Value
        microstructures{end+1} = 'delta_wave';
    end
    if cb4.Value
        microstructures{end+1} = 'cap';
    end

    % Store the selected microstructures in appdata
    setappdata(fig, 'microstructures', microstructures);

    btn.BackgroundColor = [0 1 0];
end



function runRAVEN(btn,srField, FsCAP, inputPathField, edfName, channelField, RefChannelField,...
    channelFieldCAP, RefChannelFieldCAP, microstructures, outputPathField, Patient, fileFormat,Thresholds)
    % Read input values
    if isequal(btn.BackgroundColor, [1 1 1])
        error('Please, make sure you submitted selected Microstructures before running the analysis.');
    end
    btn.BackgroundColor = [1 1 1];
    FsMicrostructures = srField.Value;
    FsCAPcycle = FsCAP.Value;
    inputPath = inputPathField.Value;
    edfFileName = edfName.Value;
    channelsmicro = split(channelField.Value, ',')'; 
    RefForMicro = split(RefChannelField.Value,',')';
    channelsCAP = split(channelFieldCAP.Value,',')';
    RefForCAP = split(RefChannelFieldCAP.Value,',')';
    if isempty(microstructures)
        error('Please, select the microstructures for analysis.');
    else
        MicroForAnalysis = microstructures;
    end
    outputPath = outputPathField.Value;
    ID = Patient.Value;
    format = fileFormat.Value;

    % Call the RAVEN function
        RAVEN(FsMicrostructures, FsCAPcycle, ...
        inputPath, edfFileName, channelsmicro,RefForMicro, ...
        channelsCAP, RefForCAP, ...
        MicroForAnalysis, outputPath, ID, format, Thresholds);

end

