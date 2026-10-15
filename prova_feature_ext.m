
% MATLAB SCRIPT: Spectral Feature Extraction (PSD) for EEG
% Calcola la potenza media nelle bande Delta, Theta1, Theta2, Alpha, Beta1,
% Beta2, Gamma per ogni canale, epoca, soggetto e condizione.



clear; clc;


chanlocs_file         = 'chanlocs.mat';  % percorso del file con i nomi dei canali
output_dir            = 'Output';        % cartella con i file *_preprocessed.mat
n_keep                = 19;              % numero di canali presenti in clean_data
use_labels_in_columns = true;            % true -> colonne tipo Fp1_Alpha | false -> Ch1_Alpha

% 1. Bande di frequenza
bands.Delta  = [1 4];
bands.Theta1 = [4.1 5.8];
bands.Theta2 = [5.9 7.4];
bands.Alpha  = [8 13];
bands.Beta1  = [13 19.9];
bands.Beta2  = [20 25];
bands.Gamma  = [24 45];

band_names = fieldnames(bands);

%Nomi dei canali da chanlocs.mat 
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
    error('Nessuna variabile con il campo "labels" trovata in %s.', chanlocs_file);
end

% Etichette di tutti i canali
if istable(chanlocs)
    all_labels = cellstr(string(chanlocs.labels))';
else
    all_labels = cellstr(string({chanlocs.labels}));
end
all_labels = all_labels(:)';

if numel(all_labels) < n_keep
    error('chanlocs ha solo %d canali, ne servono almeno %d.', numel(all_labels), n_keep);
end

% Uso solo i primi n_keep canali
if numel(all_labels) > n_keep
    fprintf('Canali scartati da chanlocs: %s\n', strjoin(all_labels(n_keep+1:end), ', '));
end
if istable(chanlocs)
    chanlocs = chanlocs(1:n_keep, :);
else
    chanlocs = chanlocs(1:n_keep);
end
channel_labels = all_labels(1:n_keep);

% Nomi validi e unici per i campi della tabella
valid_labels = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(channel_labels));

% Tabella di corrispondenza indice -> elettrodo
%T_channels = table((1:n_keep)', string(channel_labels(:)), ...
%                  'VariableNames', {'ChannelIndex','Electrode'});
%disp('Canali usati:');
%disp(T_channels);

% 2. File preprocessati
mat_files = dir(fullfile(output_dir, '*_preprocessed.mat'));
if isempty(mat_files)
    error('Nessun file _preprocessed.mat trovato nella cartella "%s".', output_dir);
end

all_features = [];
fprintf('Inizio elaborazione di %d file...\n', length(mat_files));

% 3. Ciclo sui file (un soggetto e una condizione per file)
for i = 1:length(mat_files)
    filename = fullfile(output_dir, mat_files(i).name);
    data = load(filename);

    clean_data = data.clean_data;               % [canali, campioni, epoche]
    srate      = double(data.srate);
    subject    = char(data.base_name(1));
    condition  = char(data.condition_type(1));

    [n_channels, ~, n_epochs] = size(clean_data);
    fprintf('Elaboro: %s - Condizione: %s (%d epoche)\n', subject, condition, n_epochs);

    % Controllo: numero di canali coerente con i nomi
    if n_channels ~= n_keep
        error('"%s" ha %d canali, attesi %d.', mat_files(i).name, n_channels, n_keep);
    end

    % 4. Ciclo sulle epoche
    for ep = 1:n_epochs
        epoch_data = clean_data(:, :, ep);

        feat_row = struct();
        feat_row.Subject   = subject;
        feat_row.Condition = condition;
        feat_row.Epoch     = ep;

        % 5. Ciclo sui canali
        for ch = 1:n_channels
            % PSD con metodo di Welch, finestra di 1 secondo (srate campioni)
            [pxx, freqs] = pwelch(epoch_data(ch, :), srate, [], [], srate);

            % 6. Potenza media per banda
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

% 7. Conversione in tabella
T_features = struct2table(all_features);

disp('Prime righe delle feature estratte:');
disp(head(T_features, 5));

% 8. Salvataggio
save(fullfile(output_dir, 'Extracted_Features.mat'), ...
     'T_features', 'channel_labels', 'valid_labels', 'chanlocs', 'T_channels');
writetable(T_features, fullfile(output_dir, 'Extracted_Features.csv'));
writetable(T_channels, fullfile(output_dir, 'Channel_Labels.csv'));

fprintf('\nEstrazione completata. File salvati in "%s".\n', output_dir);