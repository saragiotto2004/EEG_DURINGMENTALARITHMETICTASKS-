clear;
clc;

% The ultimate goal is to process EEG signals, extract significant features
% (physiological characteristics) and compare the two conditions to identify how 
% brain activity changes during cognitive load, specifically examining frontal(F) and 
% frontopolar(Fp) Theta channels alongside Alpha.

% STATISTICAL ANALYSIS AND COMPARISON BETWEEN CONDITIONS
% Variation of Alpha and Theta between Rest (R) and Task (T): R/T ratio,
% per subject and per channel, indicating the channel with the maximum variation.
% Note: for Alpha, I indicate the largest value of R/T (where T power decreased the most),
% while for Theta, I indicate the smallest value of R/T (where T increased the most).

output_dir='Output';
load(fullfile(output_dir, 'Extracted_Features_Theta.mat'), 'T_features', 'valid_labels');
T=T_features;

% Number of channels 
vn = T.Properties.VariableNames;
n_ch = sum(startsWith(vn, 'Ch') & endsWith(vn, '_Delta'));

% Condition: R or T 
cond = extractBefore(upper(string(T.Condition)), 2);
cond = cond(:);

% Split data into subjects
isR = cond == "R";
newSubj = [true; isR(2:end) & ~isR(1:end-1)];
subj_num = cumsum(newSubj);
n_subj = max(subj_num);

getBand = @(ch, b) T.(sprintf('%s_%s', valid_labels{ch}, b));
% Calculation and visualization
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
    
    % Alpha: maximum R/T, Theta: minimum R/T
    [~, iA]  = max(Alpha_RT);
    [~, iTh] = min(Theta_RT);
    
    Alpha_Max = strings(n_ch, 1);  Alpha_Max(iA)   = "MAX";
    Theta_Min = strings(n_ch, 1);  Theta_Min(iTh)  = "MIN";
    
    results{s} = table(Channel, Alpha_RT, Alpha_Max, Theta_RT, Theta_Min);
    
    fprintf('\n---------- SUBJECT %d ----------\n', s);
    disp(results{s});
    
    % Safety check and summary row creation
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
                      
    % Safe accumulation in the summary table
    if isempty(summary)
        summary = current_row;
    else
        summary = [summary; current_row]; 
    end
end

fprintf('\n---------- SUMMARY BY SUBJECT ----------\n');
disp(summary);

% Saving output files
xlsx_file = fullfile(output_dir, 'Alpha_Theta_RT.xlsx');
if exist(xlsx_file, 'file'), delete(xlsx_file); end

all_tbl = table();
for s = 1:n_subj
    writetable(results{s}, xlsx_file, 'Sheet', sprintf('Subject_%d', s));
    tmp = results{s};
    tmp.Subject = repmat(s, height(tmp), 1);
    all_tbl = [all_tbl; movevars(tmp, 'Subject', 'Before', 'Channel')]; 
end

writetable(summary, xlsx_file, 'Sheet', 'Summary');
writetable(all_tbl, fullfile(output_dir, 'Alpha_Theta_RT.csv'));
fprintf('\nSaved in "%s": Alpha_Theta_RT.xlsx and .csv\n', output_dir);

% RESULTS VISUALIZATION (Unified Theta [Frontal/Frontopolar] and Alpha)
ch_alpha = [11 10 13 19];
ch_theta = [2 1 6]; % Frontal and Frontopolar Theta channels
smooth_w = 1;               

% Power per epoch, averaged over selected channels (using getBand)
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

% Minimum number of epochs for each subject (between R and T) 
n_min = zeros(n_subj, 1);
for s = 1:n_subj
    nR = sum(subj_num == s & cond == "R");
    nT = sum(subj_num == s & cond == "T");
    n_min(s) = min(nR, nT);
    fprintf('Subject %d: R = %d epochs, T = %d epochs -> using first %d\n', s, nR, nT, n_min(s));
end

% Median and interquartile range
med_v = nan(numel(pot), n_subj, 2);    
iqr_v = nan(numel(pot), n_subj, 2);
cond_list = ["R", "T"];
stats_tbl = table();
for k = 1:numel(pot)
    for s = 1:n_subj
        for c = 1:2
            x  = pot{k}(subj_num == s & cond == cond_list(c));
            x  = x(1:n_min(s));
            q1 = prctile(x, 25);
            q3 = prctile(x, 75);
            med_v(k, s, c) = median(x);
            iqr_v(k, s, c) = q3 - q1;
            row = table(string(nomi{k}), s, cond_list(c), numel(x), med_v(k,s,c), q1, q3, iqr_v(k,s,c), ...
                'VariableNames', {'Banda','Soggetto','Condizione','N_epoche','Mediana','Q1','Q3','IQR'});
            stats_tbl = [stats_tbl; row];
        end
    end
end

fprintf('\n---------- MEDIAN AND IQR OF POWER (unit^2/Hz) ----------\n');
for k = 1:numel(pot)
    fprintf('\n=== %s ===\n', nomi{k});
    for s = 1:n_subj
        fprintf('Subject %d | R: median = %.4g, IQR = %.4g | T: median = %.4g, IQR = %.4g\n', ...
            s, med_v(k,s,1), iqr_v(k,s,1), med_v(k,s,2), iqr_v(k,s,2));
    end
end
disp(stats_tbl);
writetable(stats_tbl, fullfile(output_dir, 'Mediana_IQR_FFp.csv'));

% Figures for Theta and Alpha
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
        title({sprintf('Subject %d (%d epochs)', s, n_min(s)), ...
            sprintf('R: med = %.3g, IQR = %.3g', med_v(k,s,1), iqr_v(k,s,1)), ...
            sprintf('T: med = %.3g, IQR = %.3g', med_v(k,s,2), iqr_v(k,s,2))}, ...
            'FontSize', 8, 'FontWeight', 'normal');
        xlabel('Epoch'); ylabel('Power (unit^2/Hz)');
        if s == 1, legend('Location', 'best'); end
        grid on; box off;
    end

    if numel(canali{k}) == n_ch
        ch_str = 'all channels';
    else
        selected_names = valid_labels(canali{k});
        ch_str = ['channels ' char(strjoin(string(selected_names), ', '))];
    end

    sgtitle(sprintf('%s - power per epoch, R vs T (%s)', nomi{k}, ch_str));
end
