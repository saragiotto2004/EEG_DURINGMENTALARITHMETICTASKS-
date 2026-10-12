% =========================================================================
% MATLAB SCRIPT: Spectral Feature Extraction (PSD) for EEG
% Objective: Calculate the power in frequency bands (Delta, Theta, Alpha, Beta, Gamma)
% for each channel, epoch, subject, and condition in the Output folder.
% =========================================================================

clear; clc;

% 1. Definition of frequency bands of interest (in Hz)
bands.Delta = [1 4];    % Deep sleep / artifacts monitoring
bands.Theta1 = [4.1 5.8];   % General attention / relaxed alertness
bands.Theta2 = [5.9 7.4];   % Working memory and active cognitive control
bands.Alpha = [8 13];   % Attention and task-related modulation (e.g., mental arithmetic)
bands.Beta1  = [13 19.9];  % Sustained alert attention and focused concentration
bands.Beta2 = [20 25];  % High cognitive load / intense problem solving
bands.Gamma = [24 45];  % Higher cognitive processes

band_names = fieldnames(bands);

% 2. Definition of the path to the preprocessing output folder
output_dir = 'Output'; 

% Search for all preprocessed files ending with '_preprocessed.mat'
mat_files = dir(fullfile(output_dir, '*_preprocessed.mat'));

if isempty(mat_files)
    error('No _preprocessed.mat files found in folder "%s". Check the path.', output_dir);
end

% Initialization of the structure that will hold all extracted features
all_features = [];

fprintf('Starting processing of %d files...\n', length(mat_files));

% 3. Iteration loop over each file (corresponding to a subject and a condition)
for i = 1:length(mat_files)
    filename = fullfile(output_dir, mat_files(i).name);
    data = load(filename);
    
    % Extraction of variables saved in the .mat file (converted to char for safety)
    clean_data = data.clean_data; % 3D Matrix: [n_channels, n_samples, n_epochs]
    srate = double(data.srate);   % Sampling rate
    subject = char(data.base_name(1));  % Subject identifier (e.g., 'Subject01_1')
    condition = char(data.condition_type(1)); % Condition (e.g., 'RESTING' or 'TASK')
    
    [n_channels, n_samples, n_epochs] = size(clean_data);
    fprintf('Processing: %s - Condition: %s (%d epochs)\n', subject, condition, n_epochs);
    
    % 4. Loop over each available epoch for the subject/condition
    for ep = 1:n_epochs
        epoch_data = clean_data(:, :, ep); % Extraction of current epoch data [n_channels x n_samples]
        
        % Creation of data structure for the current row
        feat_row.Subject = subject;
        feat_row.Condition = condition;
        feat_row.Epoch = ep;
        
        % 5. Loop over each EEG channel
        for ch = 1:n_channels
            % Calculation of Power Spectral Density (PSD) using Welch's method
            % Analysis window equal to 1 second (srate samples)
            [pxx, freqs] = pwelch(epoch_data(ch, :), srate, [], [], srate);
            
            % 6. Extraction of mean power for each frequency band
            for b = 1:length(band_names)
                b_name = band_names{b};
                f_low = bands.(b_name)(1);
                f_high = bands.(b_name)(2);
                
                % Selection of frequencies falling within the current band
                idx_band = (freqs >= f_low) & (freqs <= f_high);
                band_power = mean(pxx(idx_band)); % Mean power in the band
                
                % Dynamic saving of the feature (e.g., Ch1_Delta, Ch1_Alpha, etc.)
                field_str = sprintf('Ch%d_%s', ch, b_name);
                feat_row.(field_str) = band_power;
            end
        end
        
        % Appending the row to the main structure array
        if isempty(all_features)
            all_features = feat_row;
        else
            all_features(end+1) = feat_row;
        end
    end
end

% 7. Conversion of the structure into a MATLAB Table for easier statistical analysis
T_features = struct2table(all_features);

% Display the first few rows of the resulting table in the Command Window
disp('First rows of extracted features:');
disp(head(T_features, 5));

% 8. Saving the feature dataset in both .mat and .csv formats inside the Output folder
save(fullfile(output_dir, 'Extracted_Features.mat'), 'T_features');
writetable(T_features, fullfile(output_dir, 'Extracted_Features.csv'));

fprintf('\nExtraction completed successfully! Files saved in "%s/Extracted_Features.csv".\n', output_dir);