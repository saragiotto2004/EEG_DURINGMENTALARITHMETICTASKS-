clear;
clc;

% L'obiettivo finale è processare i segnali EEG, estrarre delle feature significative
% (caratteristiche fisiologiche) e confrontare le due condizioni per identificare come 
% cambia l'attività cerebrale durante il carico cognitivo.

% ANALISI STATISTICA E CONFRONTO TRA LE CONDIZIONI
% Variazione di Alpha e Theta tra Rest (R) e Task (T): rapporto R/T,
% per soggetto e per canale, con indicazione del canale che presenta la massima variazione.
% N.B. per Alpha indico il valore più grande di R/T (dove la potenza di T è diminuita di più),
% mentre per Theta indico il valore più piccolo di R/T (dove T è aumentata di più).

output_dir = 'Output';
load(fullfile(output_dir, 'Extracted_Features_Theta.mat'), 'T_features', 'valid_labels');
T = T_features;

% Numero di canali 
vn = T.Properties.VariableNames;
n_ch = sum(startsWith(vn, 'Ch') & endsWith(vn, '_Delta'));

% Condizione: R o T 
cond = extractBefore(upper(string(T.Condition)), 2);
cond = cond(:);

% Divido i dati nei soggetti
isR = cond == "R";
newSubj = [true; isR(2:end) & ~isR(1:end-1)];
subj_num = cumsum(newSubj);
n_subj = max(subj_num);

getBand = @(ch, b) T.(sprintf('%s_%s', valid_labels{ch}, b));
% Calcolo e visualizzazione
results = cell(n_subj, 1);
summary = table();

for s = 1:n_subj
    iR = (subj_num == s) & (cond == "R");
    iT = (subj_num == s) & (cond == "T");
    
    Channel   = strings(n_ch, 1);
    Alpha_RT  = nan(n_ch, 1);
    Theta_RT  = nan(n_ch, 1);
    
    for ch = 1:n_ch
        alpha  = getBand(ch, 'Alpha');
        theta  = getBand(ch, 'Theta');
        
        Channel(ch)  = "Ch" + ch;
        Alpha_RT(ch) = mean(alpha(iR)) / mean(alpha(iT));    
        Theta_RT(ch) = mean(theta(iR)) / mean(theta(iT));  
    end
    
    % Alpha: R/T massimo, Theta: R/T minimo
    [~, iA]  = max(Alpha_RT);
    [~, iTh] = min(Theta_RT);
    
    Alpha_Max = strings(n_ch, 1);  Alpha_Max(iA)   = "MAX";
    Theta_Min = strings(n_ch, 1);  Theta_Min(iTh)  = "MIN";
    
    results{s} = table(Channel, Alpha_RT, Alpha_Max, Theta_RT, Theta_Min);
    
    fprintf('\n---------- SOGGETTO %d ----------\n', s);
    disp(results{s});
    
    % Controllo di sicurezza e creazione riga di sommario
    if isempty(iA) || isempty(iTh)
        ch_a = ""; val_a = NaN;
        ch_t = ""; val_t = NaN;
    else
        ch_a = Channel(iA);  val_a = Alpha_RT(iA);
        ch_t = Channel(iTh); val_t = Theta_RT(iTh);
    end
    
    current_row = table(s, {ch_a}, val_a, {ch_t}, val_t, ...
        'VariableNames', {'Subject', 'Alpha_CanaleMax', 'Alpha_RT_Max', ...
                          'Theta_CanaleMin', 'Theta_RT_Min'});
                      
    % Accumulo sicuro nella tabella di riepilogo (corretto qui)
    if isempty(summary)
        summary = current_row;
    else
        summary = [summary; current_row]; %#ok<AGROW>
    end
end

fprintf('\n---------- RIEPILOGO PER SOGGETTO ----------\n');
disp(summary);

% Salvataggio dei file di output
xlsx_file = fullfile(output_dir, 'Alpha_Theta_RT.xlsx');
if exist(xlsx_file, 'file'), delete(xlsx_file); end

all_tbl = table();
for s = 1:n_subj
    writetable(results{s}, xlsx_file, 'Sheet', sprintf('Soggetto_%d', s));
    tmp = results{s};
    tmp.Subject = repmat(s, height(tmp), 1);
    all_tbl = [all_tbl; movevars(tmp, 'Subject', 'Before', 'Channel')]; %#ok<AGROW>
end

writetable(summary, xlsx_file, 'Sheet', 'Riepilogo');
writetable(all_tbl, fullfile(output_dir, 'Alpha_Theta_RT.csv'));
fprintf('\nSalvato in "%s": Alpha_Theta_RT.xlsx e .csv\n', output_dir);

% VISUALIZZAZIONE RISULTATI (Theta unico e Alpha)
ch_alpha = [11 10 13 19];
ch_theta = [10 9 19]; % Centrali e Parietali
smooth_w = 1;               

% Potenza per epoca, media sui canali selezionati (usando getBand)
theta_p = zeros(height(T), 1);
for ch = ch_theta
    theta_p = theta_p + getBand(ch, 'Theta') / numel(ch_theta);
end

alpha_p = zeros(height(T), 1);
for ch = ch_alpha
    alpha_p = alpha_p + getBand(ch, 'Alpha') / numel(ch_alpha);
end

pot    = {theta_p, alpha_p};
nomi   = {'Theta', 'Alpha'};
canali = {ch_theta, ch_alpha};

colR = [0.2 0.5 0.8];   
colT = [0.85 0.35 0.2];

% Numero minimo di epoche per ogni soggetto (tra R e T) 
n_min = zeros(n_subj, 1);
for s = 1:n_subj
    nR = sum(subj_num == s & cond == "R");
    nT = sum(subj_num == s & cond == "T");
    n_min(s) = min(nR, nT);
    fprintf('Soggetto %d: R = %d epoche, T = %d epoche -> uso le prime %d\n', s, nR, nT, n_min(s));
end

% Figure per Theta e Alpha
for k = 1:numel(pot)
    figure('Color', 'w', 'Position', [100 100 1100 600], 'Name', nomi{k});
    for s = 1:n_subj
        subplot(2, 3, s); hold on;
        pR = pot{k}(subj_num == s & cond == "R");
        pT = pot{k}(subj_num == s & cond == "T");
        pR = pR(1:n_min(s));
        pT = pT(1:n_min(s));
        
        plot(movmean(pR, smooth_w), '-o', 'Color', colR, 'MarkerSize', 3, 'DisplayName', 'R');
        plot(movmean(pT, smooth_w), '-o', 'Color', colT, 'MarkerSize', 3, 'DisplayName', 'T');
        yline(mean(pR), '--', 'Color', colR, 'HandleVisibility', 'off');
        yline(mean(pT), '--', 'Color', colT, 'HandleVisibility', 'off');
        
        xlim([1 n_min(s)]);
        title(sprintf('Soggetto %d (%d epoche)', s, n_min(s)));
        xlabel('Epoca'); ylabel('Potenza (unit^2/Hz)');
        if s == 1, legend('Location', 'best'); end
        grid on; box off;
    end
    
    if numel(canali{k}) == n_ch
        ch_str = 'tutti i canali';
    else
        selected_names = valid_labels(canali{k});
        ch_str = ['canali ' char(strjoin(string(selected_names), ', '))];
    end
    
    sgtitle(sprintf('%s - potenza per epoca, R vs T (%s)', nomi{k}, ch_str));
end