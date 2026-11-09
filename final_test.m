clear;
clc;
% The ultimate goal is to process EEG signals, extract significant features
% (physiological characteristics) and compare the two conditions to identify how 
% brain activity changes during cognitive load.

% STATISTICAL ANALYSIS AND COMPARISON BETWEEN CONDITIONS
% Variation of Alpha, Theta1, and Theta2 between Rest (R) and Task (T): R/T ratio,
% per subject and per channel, indicating the channel with the maximum variation.
% Note: for alpha, I indicate the largest value of R/T because it tells me where the power
% of T decreased the most, the opposite for theta1 and theta2.

output_dir='Output';
load(fullfile(output_dir, 'Extracted_Features.mat'), 'T_features', 'channel_labels', 'valid_labels');
T=T_features;

% Channels 
n_ch=numel(channel_labels);
getBand= @(ch, b) T.(sprintf('%s_%s', valid_labels{ch}, b));

% Condition: R or T
cond=extractBefore(upper(string(T.Condition)), 2);
cond=cond(:);

% Split data into subjects 
isR      =cond == "R";
newSubj  =[true; isR(2:end) & ~isR(1:end-1)];
subj_num =cumsum(newSubj);
n_subj   =max(subj_num);
fprintf('Subjects identified: %d\n', n_subj);

% Calculation and visualization
results=cell(n_subj, 1);
summary=table();

for s = 1:n_subj
    iR=(subj_num == s) & (cond == "R");
    iT=(subj_num == s) & (cond == "T");

    Channel   =strings(n_ch, 1);
    Alpha_RT  =nan(n_ch, 1);
    Theta1_RT =nan(n_ch, 1);
    Theta2_RT =nan(n_ch, 1);

    for ch = 1:n_ch
        alpha  =getBand(ch, 'Alpha');
        theta1 =getBand(ch, 'Theta1');
        theta2 =getBand(ch, 'Theta2');

        Channel(ch)   =string(channel_labels{ch});
        Alpha_RT(ch)  =mean(alpha(iR))  / mean(alpha(iT));
        Theta1_RT(ch) =mean(theta1(iR)) / mean(theta1(iT));
        Theta2_RT(ch) =mean(theta2(iR)) / mean(theta2(iT));
    end

    % Alpha: maximum R/T | Theta1 and Theta2: minimum R/T
    [~, iA]   = max(Alpha_RT);
    [~, iTh1] = min(Theta1_RT);
    [~, iTh2] = min(Theta2_RT);

    Alpha_Max  = strings(n_ch, 1);  Alpha_Max(iA)    = "MAX";
    Theta1_Min = strings(n_ch, 1);  Theta1_Min(iTh1) = "MIN";
    Theta2_Min = strings(n_ch, 1);  Theta2_Min(iTh2) = "MIN";

    results{s} = table(Channel, Alpha_RT, Alpha_Max, ...
                       Theta1_RT, Theta1_Min, Theta2_RT, Theta2_Min);

    fprintf('\n---------- SUBJECT %d ----------\n', s);
    disp(results{s});

    summary = [summary; table(s, Channel(iA), Alpha_RT(iA), ...
        Channel(iTh1), Theta1_RT(iTh1), Channel(iTh2), Theta2_RT(iTh2), ...
        'VariableNames', {'Subject','Alpha_CanaleMax','Alpha_RT_Max', ...
                          'Theta1_CanaleMin','Theta1_RT_Min', ...
                          'Theta2_CanaleMin','Theta2_RT_Min'})];
end

fprintf('\n---------- SUMMARY BY SUBJECT ----------\n');
disp(summary);

% Saving
xlsx_file =fullfile(output_dir, 'Alpha_Theta_RT.xlsx');
if exist(xlsx_file, 'file'), delete(xlsx_file); end

all_tbl=table();
for s = 1:n_subj
    writetable(results{s}, xlsx_file, 'Sheet', sprintf('Subject_%d', s));
    tmp=results{s};
    tmp.Subject=repmat(s, height(tmp), 1);
    all_tbl=[all_tbl; movevars(tmp, 'Subject', 'Before', 'Channel')]; 
end
writetable(summary, xlsx_file, 'Sheet', 'Summary');
writetable(all_tbl, fullfile(output_dir, 'Alpha_Theta_RT.csv'));

fprintf('\nSaved in "%s": Alpha_Theta_RT.xlsx and .csv\n', output_dir);



% Trend of band power across epochs, R vs T, per subject.
% For each subject, the R and T series are truncated to the minimum number of epochs
% between its two conditions.
% Selected channels (from Seleznov literature)
%   Alpha  : Pz, P3, P4, O1, O2
%   Theta1 : Fz, F3, F4, Cz
%   Theta2 : same as alpha

names_alpha={'Pz','P3','P4','O1','O2'};
names_theta1={'Fz','F3','F4','Cz'};
names_theta2={'Pz','P3','P4','O1','O2'};
smooth_w=1; % moving average window = 1 -> no smoothing applied

% Name -> index conversion, with a check that all channels exist
nameToIdx= @(nm) cellfun(@(x) find(strcmpi(channel_labels, x), 1), nm, 'UniformOutput', false);
chk= @(nm) assert(all(~cellfun(@isempty, nameToIdx(nm))), ...
    'Channel not found among: %s', strjoin(channel_labels, ', '));
chk(names_alpha); chk(names_theta1); chk(names_theta2);

ch_alpha=cell2mat(nameToIdx(names_alpha));
ch_theta1=cell2mat(nameToIdx(names_theta1));
ch_theta2=cell2mat(nameToIdx(names_theta2));

% Power per epoch, averaged over selected channels: for each band, for each epoch, computes the arithmetic mean
% of the selected channels
theta1_p=zeros(height(T), 1);
for ch=ch_theta1
    theta1_p=theta1_p + getBand(ch, 'Theta1') / numel(ch_theta1);
end
theta2_p=zeros(height(T), 1);
for ch=ch_theta2
    theta2_p=theta2_p + getBand(ch, 'Theta2') / numel(ch_theta2);
end
alpha_p=zeros(height(T), 1);
for ch=ch_alpha
    alpha_p=alpha_p + getBand(ch, 'Alpha') / numel(ch_alpha);
end

pot    ={theta1_p, theta2_p, alpha_p};
nomi   ={'Theta1', 'Theta2', 'Alpha'}; % Outputs three vectors, theta1_p, theta2_p, alpha_p, as long as the rows of the table.
canali ={names_theta1, names_theta2, names_alpha};
colR=[0.2 0.5 0.8];
colT=[0.85 0.35 0.2];

% Minimum number of epochs for each subject (between R and T) 
n_min=zeros(n_subj, 1);
for s= 1:n_subj
    nR=sum(subj_num == s & cond == "R");
    nT=sum(subj_num == s & cond == "T");
    n_min(s)=min(nR, nT);
    fprintf('Subject %d: R = %d epochs, T = %d epochs -> using first %d\n', s, nR, nT, n_min(s));
end

% Figures 
for k= 1:3
    figure('Color','w','Position',[100 100 1100 600], 'Name', nomi{k});
    for s= 1:n_subj
        subplot(2,3,s); hold on;
        pR=pot{k}(subj_num==s & cond=="R");
        pT=pot{k}(subj_num==s & cond=="T");
        pR=pR(1:n_min(s));
        pT=pT(1:n_min(s));

        plot(movmean(pR, smooth_w), '-o', 'Color', colR, 'MarkerSize', 3, 'DisplayName', 'R');
        plot(movmean(pT, smooth_w), '-o', 'Color', colT, 'MarkerSize', 3, 'DisplayName', 'T');
        yline(mean(pR), '--', 'Color', colR, 'HandleVisibility', 'off');
        yline(mean(pT), '--', 'Color', colT, 'HandleVisibility', 'off');

        xlim([1 max(n_min(s), 2)]);
        title(sprintf('Subject %d (%d epochs)', s, n_min(s)));
        xlabel('Epoch'); ylabel('Power (unit^2/Hz)');
        if s == 1, legend('Location', 'best'); end
        grid on; box off;
    end
    sgtitle(sprintf('%s - power per epoch, R vs T (channels %s)', ...
        nomi{k}, strjoin(canali{k}, ', ')));
end
