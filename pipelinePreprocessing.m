%% =========================================================================
%% PREPROCESSING EEG PIPELINE - BIOLOGICAL SIGNAL AND IMAGE PROCESSING PROJECT
%% =========================================================================
clear; clc; close all;

%% 1. PATH CONFIGURATION
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
cd(script_dir);

data_folder   = fullfile(script_dir, 'Data');
helper_folder = fullfile(script_dir, 'helper_code');
output_folder = fullfile(script_dir, 'Output');

if exist(helper_folder, 'dir')
    addpath(helper_folder);
else
    error('Folder "helper_code" not found!');
end

if ~exist(output_folder, 'dir'), mkdir(output_folder); end

%% 2. AUTOMATIC EEGLAB INITIALIZATION
if ispc, user_home = getenv('USERPROFILE'); else, user_home = getenv('HOME'); end
desktop_dir = fullfile(user_home, 'Desktop');

eeglab_in_project = dir(fullfile(script_dir, 'eeglab*'));
eeglab_in_desktop = dir(fullfile(desktop_dir, 'eeglab*'));

if ~isempty(eeglab_in_project)
    eeglab_path = fullfile(script_dir, eeglab_in_project(1).name);
elseif ~isempty(eeglab_in_desktop)
    eeglab_path = fullfile(desktop_dir, eeglab_in_desktop(1).name);
else
    eeglab_path = '';
end

if exist('eeglab', 'file')
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');
elseif ~isempty(eeglab_path) && exist(eeglab_path, 'dir')
    addpath(eeglab_path);
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui'); 
else
    error('EEGLAB not found!');
end

%% 3. LOADING CHANLOCS (Electrode coordinates)
chanlocs_path = fullfile(helper_folder, 'chanlocs.mat');
if exist(chanlocs_path, 'file')
    chanlocs_data = load(chanlocs_path);
    fields = fieldnames(chanlocs_data);
    chanlocs = chanlocs_data.(fields{1});
    fprintf('chanlocs.mat loaded successfully.\n');
else
    error('File chanlocs.mat not found in the helper_code folder!');
end

%% 4. PREPROCESSING PARAMETERS
srate     = 500;  % Hz
epoch_dur = 2;    % seconds per sub-epoch
v_thresh  = 80;   % uV

subject_files = dir(fullfile(data_folder, 'Subject*.mat'));
if isempty(subject_files), error('No .mat files found in the Data folder!'); end

fprintf('Found %d data files to process.\n\n', length(subject_files));
rejection_summary = table();
eeglab; close;

%% 5. PREPROCESSING LOOP (with multi-channel reconstruction)
for i = 1:length(subject_files)
    file_name = subject_files(i).name;
    file_path = fullfile(data_folder, file_name);
    [~, base_name, ~] = fileparts(file_name);
    
    if contains(base_name, '_1')
        condition_type = 'RESTING'; 
        target_duration = 180; % 180s Resting State
    else
        condition_type = 'TASK'; 
        target_duration = 60;  % 60s Mental Task
    end
    
    fprintf('[%d/%d] Preprocessing [%s]: %s...\n', i, length(subject_files), condition_type, file_name);
    
    % --- 1. RECONSTRUCTING SEPARATE CHANNELS INTO A 2D MATRIX  ---
    subj_data = load(file_path);
    var_names = fieldnames(subj_data);
    
    data_matrix = [];
    for ch = 1:length(var_names)
        ch_data = double(subj_data.(var_names{ch}));
        if size(ch_data, 1) > size(ch_data, 2)
            ch_data = ch_data';
        end
        data_matrix = [data_matrix; ch_data]; %#ok<AGROW>
    end
    
    % --- 2. IMPORTING INTO EEGLAB WITH CHANLOCS CORRESPONDENCE ---
    EEG = pop_importdata('dataformat', 'array', ...
                         'data', data_matrix, ...
                         'srate', srate, ...
                         'nbchan', size(data_matrix, 1));

    % CHANNEL DATA ASSIGNMENT AND CORRESPONDENCE CHECK
    EEG.chanlocs = chanlocs(1:size(data_matrix, 1));
    EEG = eeg_checkset(EEG);

    % --- 3. PREPROCESSING STEPS ---
    % Re-referencing (Average Reference)
    EEG = pop_reref(EEG, []);

    % Band-pass (0.5 - 45 Hz) and Notch (50 Hz)
    EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5, 'hicutoff', 45);
    EEG = pop_eegfiltnew(EEG, 'locutoff', 49, 'hicutoff', 51, 'revfilt', 1);


    % --- 4. DATA TRIMMING AND EPOCHING (2s Epochs) ---
    max_pnts = min(EEG.pnts, target_duration * srate);
    data_trimmed = EEG.data(:, 1:max_pnts);
    
    pnts_epoch = epoch_dur * srate;
    n_epochs = floor(max_pnts / pnts_epoch);
    data_trimmed = data_trimmed(:, 1:(n_epochs * pnts_epoch));
    
    data_3d = reshape(data_trimmed, [EEG.nbchan, pnts_epoch, n_epochs]);
    
    % Rejection of sub-epochs > +/- 80 uV
    valid_idx = squeeze(all(all(abs(data_3d) <= v_thresh, 1), 2));
    clean_data = data_3d(:, :, valid_idx);
    
    num_epochs_kept = size(clean_data, 3);
    rejected_epochs = n_epochs - num_epochs_kept;
    rejection_rate  = (rejected_epochs / n_epochs) * 100;

    % RECORDING STATISTICS
    new_entry = table({base_name}, {condition_type}, n_epochs, num_epochs_kept, ...
                      rejected_epochs, rejection_rate, ...
                      'VariableNames', {'Subject', 'Condition', 'InitialEpochs', 'KeptEpochs', 'RejectedEpochs', 'RejectionRatePct'});
    rejection_summary = [rejection_summary; new_entry];

    % --- 5. SAVING PREPROCESSED DATA  ---
    mat_output_path = fullfile(output_folder, sprintf('%s_preprocessed.mat', base_name));
    save(mat_output_path, 'clean_data', 'base_name', 'condition_type', 'srate', 'num_epochs_kept');

    % Optional Visual Quality Control (QC): Generate and save mean EEG power topoplot
    if num_epochs_kept > 0
        mean_power = mean(mean(clean_data.^2, 2), 3);
        f = figure('Visible', 'off');
        topoplot(mean_power, EEG.chanlocs, 'maplimits', 'maxmin', 'electrodes', 'on');
        title(sprintf('EEG Power [%s] - %s', condition_type, base_name), 'Interpreter', 'none');
        colorbar;
        saveas(f, fullfile(output_folder, sprintf('%s_topoplot.png', base_name)));
        close(f);
    end

    fprintf('  -> Completed! %d clean epochs out of %d (Rejected: %.1f%%).\n\n', ...
            num_epochs_kept, n_epochs, rejection_rate);
end

% Saving report
writetable(rejection_summary, fullfile(output_folder, 'Rejection_Epochs_Report.csv'));

fprintf('=================================================================\n');
fprintf('Preprocessing completed successfully for all subjects!\n');
fprintf('=================================================================\n');
