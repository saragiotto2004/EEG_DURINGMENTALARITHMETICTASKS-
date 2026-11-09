% MATLAB SCRIPT: Spectral Feature Extraction (PSD) for EEG
% Computes the average power in the Delta, Theta, Alpha, Beta1,
% Beta2, Gamma bands for each channel, epoch, subject, and condition.
clear; clc;

chanlocs_file         = 'chanlocs.mat';    % path to the file containing channel names
output_dir            = 'Output';          % folder containing the *_preprocessed.mat files
n_keep                = 19;                % number of channels present in clean_data
use_labels_in_columns = true;              % true -> column type Fp1_Alpha | false -> Ch1_Alpha

% 1. Frequency bands (Unified Theta instead of Theta1/Theta2)
bands.Delta  = [1 4];
bands.Theta  = [4 8];    % Standard Theta band (4 - 8 Hz)
bands.Alpha  = [8 13];
bands.Beta1  = [13 19.9];
bands.Beta2  = [20 25];
bands.Gamma  = [24 45];

band_names = fieldnames(bands);

% Channel names from chanlocs.mat 
S = load(chanlocs_file);
chanlocs = [];
if isfield(S, 'chanlocs')
    chanlocs = S.chanlocs;
elseif isfield(S, 'EEG') && isfield(S.EEG, 'chanlocs')
    chanlocs = S.EEG.chanlocs;
else
    fn = fieldnames(S);
    for k = 1:numel(fn)
        v = S.(fn{k});
        if (isstruct(v) && isfield(v, 'labels')) || ...
           (istable(v) && ismember('labels', v.Properties.VariableNames))
            chanlocs = v; break
        end
    end
end

if isempty(chanlocs)
    error('No variable with the "labels" found in %s.', chanlocs_file);
end

% Labels of all channels
if istable(chanlocs)
    all_labels = cellstr(string(chanlocs.labels))';
else
    all_labels = cellstr(string({chanlocs.labels}));
end
all_labels = all_labels(:)';

if numel(all_labels) < n_keep
    error('chanlocs has only %d channels, at least %d are needed.', numel(all_labels), n_keep);
end

% Use only the first n_keep channels
if numel(all_labels) > n_keep
    fprintf('Channels discarded from chanlocs: %s\n', strjoin(all_labels(n_keep+1:end), ', '));
end

if istable(chanlocs)
    chanlocs = chanlocs(1:n_keep, :);
else
    chanlocs = chanlocs(1:n_keep);
end
channel_labels = all_labels(1:n_keep);

% Valid and unique names for table fields
valid_labels = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(channel_labels));

% Correspondence table index -> electrode
T_channels = table((1:n_keep)', string(channel_labels(:)), ...
                  'VariableNames', {'ChannelIndex','Electrode'});
disp('Channels used:');
disp(T_channels);

% 2. Preprocessed files
mat_files = dir(fullfile(output_dir, '*_preprocessed.mat'));
if isempty(mat_files)
    error('No _preprocessed.mat files found in the folder "%s".', output_dir);
end

all_features = [];
fprintf('Starting processing of %d files...\n', length(mat_files));

% 3. Loop over files (one subject and one condition per file)
for i = 1:length(mat_files)
    filename = fullfile(output_dir, mat_files(i).name);
    data = load(filename);
    
    clean_data = data.clean_data;              % [channels, samples, epochs]
    srate      = double(data.srate);
    subject    = char(data.base_name(1));
    condition  = char(data.condition_type(1));
    
    [n_channels, ~, n_epochs] = size(clean_data);
    fprintf('Processing: %s - Condition: %s (%d epochs)\n', subject, condition, n_epochs);
    
    % Check: number of channels consistent with names
    if n_channels ~= n_keep
        error('"%s" has %d channels, expected %d.', mat_files(i).name, n_channels, n_keep);
    end
    
    % 4. Loop over epochs
    for ep = 1:n_epochs
        epoch_data = clean_data(:, :, ep);
        
        feat_row = struct();
        feat_row.Subject   = subject;
        feat_row.Condition = condition;
        feat_row.Epoch     = ep;
        
        % 5. Loop over channels
        for ch = 1:n_channels
            % PSD using Welch's method, 1-second window (srate samples)
            [pxx, freqs] = pwelch(epoch_data(ch, :), srate, [], [], srate);
            
            % 6. Average power per band
            for b = 1:length(band_names)
                b_name = band_names{b};
                f_low  = bands.(b_name)(1);
                f_high = bands.(b_name)(2);
                
                idx_band   = (freqs >= f_low) & (freqs <= f_high);
                band_power = mean(pxx(idx_band));
                
                if use_labels_in_columns
                    field_str = sprintf('%s_%s', valid_labels{ch}, b_name);
                else
                    field_str = sprintf('Ch%d_%s', ch, b_name);
                end
                
                feat_row.(field_str) = band_power;
            end
        end
        
        if isempty(all_features)
            all_features = feat_row;
        else
            all_features(end+1) = feat_row; %#ok<SAGROW>
        end
    end
end

% 7. Conversion to table
T_features = struct2table(all_features);
disp('First rows of extracted features:');
disp(head(T_features, 5));

% 8. Saving
save(fullfile(output_dir, 'Extracted_Features_Theta.mat'), ...
     'T_features', 'channel_labels', 'valid_labels', 'chanlocs', 'T_channels');
writetable(T_features, fullfile(output_dir, 'Extracted_Features_Theta.csv'));
writetable(T_channels, fullfile(output_dir, 'Channel_Labels_Theta.csv'));

fprintf('\nExtraction completed. Files saved in "%s".\n', output_dir);
