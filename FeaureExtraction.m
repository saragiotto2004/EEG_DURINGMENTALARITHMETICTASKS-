% =========================================================================
% SCRIPT MATLAB: Estrazione delle Feature Spettrali (PSD) per EEG
% Obiettivo: Calcolare la potenza nelle bande di frequenza (Theta, Alpha, Beta, Gamma)
% per ciascun canale, epoca, soggetto e condizione nella cartella Output.
% =========================================================================

clear; clc;

% 1. Definizione delle bande di frequenza di interesse (in Hz)
bands.Theta = [4 8];    % Memoria di lavoro e controllo cognitivo
bands.Alpha = [8 13];   % Attenzione e modulazione legata al compito (es. calcolo mentale)
bands.Beta  = [13 30];  % Attivazione corticale e concentrazione
bands.Gamma = [30 45];  % Processi cognitivi superiori

band_names = fieldnames(bands);

% 2. Definizione del percorso della cartella di output del preprocessing
output_dir = 'Output'; 

% Ricerca di tutti i file prepreprocessati che terminano con '_preprocessed.mat'
mat_files = dir(fullfile(output_dir, '*_preprocessed.mat'));

if isempty(mat_files)
    error('Nessun file _preprocessed.mat trovato nella cartella "%s". Verifica il percorso.', output_dir);
end

% Inizializzazione della struttura che conterrà tutte le feature estratte
all_features = [];

fprintf('Inizio elaborazione di %d file...\n', length(mat_files));

% 3. Ciclo di iterazione su ogni file (corrispondente a un soggetto e a una condizione)
for i = 1:length(mat_files)
    filename = fullfile(output_dir, mat_files(i).name);
    data = load(filename);
    
    % Estrazione delle variabili salvate nel file .mat (convertite in stringa/char per sicurezza)
    clean_data = data.clean_data; % Matrice 3D: [n_channels, n_samples, n_epochs]
    srate = double(data.srate);   % Frequenza di campionamento
    subject = char(data.base_name(1));  % Identificativo del soggetto (es. 'Subject01_1')
    condition = char(data.condition_type(1)); % Condizione (es. 'RESTING' o 'TASK')
    
    [n_channels, n_samples, n_epochs] = size(clean_data);
    fprintf('Elaborazione: %s - Condizione: %s (%d epoche)\n', subject, condition, n_epochs);
    
    % 4. Ciclo su ogni epoca disponibile per il soggetto/condizione
    for ep = 1:n_epochs
        epoch_data = clean_data(:, :, ep); % Estrazione dei dati dell'epoca corrente [n_channels x n_samples]
        
        % Creazione della struttura dati per la riga corrente
        feat_row.Subject = subject;
        feat_row.Condition = condition;
        feat_row.Epoch = ep;
        
        % 5. Ciclo su ogni canale EEG
        for ch = 1:n_channels
            % Calcolo della Densità Spettrale di Potenza (PSD) tramite il metodo di Welch
            % Finestra di analisi pari a 1 secondo (srate campioni)
            [pxx, freqs] = pwelch(epoch_data(ch, :), srate, [], [], srate);
            
            % 6. Estrazione della potenza media per ciascuna banda di frequenza
            for b = 1:length(band_names)
                b_name = band_names{b};
                f_low = bands.(b_name)(1);
                f_high = bands.(b_name)(2);
                
                % Selezione delle frequenze che rientrano nella banda corrente
                idx_band = (freqs >= f_low) & (freqs <= f_high);
                band_power = mean(pxx(idx_band)); % Potenza media nella banda
                
                % Salvataggio dinamico della feature (es. Ch1_Alpha, Ch2_Beta, ecc.)
                field_str = sprintf('Ch%d_%s', ch, b_name);
                feat_row.(field_str) = band_power;
            end
        end
        
        % Accodamento della riga all'array di strutture principale
        if isempty(all_features)
            all_features = feat_row;
        else
            all_features(end+1) = feat_row;
        end
    end
end

% 7. Conversione della struttura in una Tabella MATLAB per facilitare l'analisi statistica
T_features = struct2table(all_features);

% Mostra le prime righe della tabella risultante nella Command Window
disp('Prime righe delle feature estratte:');
disp(head(T_features, 5));

% 8. Salvataggio del dataset di feature sia in formato .mat che .csv all'interno della cartella Output
save(fullfile(output_dir, 'Extracted_Features.mat'), 'T_features');
writetable(T_features, fullfile(output_dir, 'Extracted_Features.csv'));

fprintf('\nEstrazione completata con successo! File salvati in "%s/Extracted_Features.csv".\n', output_dir);