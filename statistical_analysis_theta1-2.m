
%% TEST THETA 1 AND THETA 2 

%% STATISTICAL ANALYSIS SCRIPT: THETA1, THETA2, AND ALPHA
% Wilcoxon Rank-Sum (intra-subject across epochs) and Signed-Rank (group N = 6)
clear; clc; close all;

%% 1. Data Loading and Environment Setup
output_dir = 'Output';
mat_file   = fullfile(output_dir, 'Extracted_Features.mat');

if ~exist(mat_file, 'file')
    error('File not found in "%s". Please run feature extraction first.', mat_file);
end

% Load extracted features table and channel labels
load(mat_file, 'T_features', 'channel_labels', 'valid_labels');
T = T_features;

n_ch    = numel(channel_labels);
getBand = @(ch, b) T.(sprintf('%s_%s', valid_labels{ch}, b));

% Condition extraction: 'R' (Rest) or 'T' (Task)
cond = extractBefore(upper(string(T.Condition)), 2);
cond = cond(:);

% Identify and segment individual subjects based on condition transitions
isR      = cond == "R";
newSubj  = [true; isR(2:end) & ~isR(1:end-1)];
subj_num = cumsum(newSubj);
n_subj   = max(subj_num);

%% 2. Region of Interest (ROI) Definitions from Seleznov Literature
names_theta1 = {'Fz', 'F3', 'F4', 'Cz'};
names_theta2 = {'Pz', 'P3', 'P4', 'O1', 'O2'};
names_alpha  = {'Pz', 'P3', 'P4', 'O1', 'O2'};

% Map channel labels to column indices
nameToIdx = @(nm) cellfun(@(x) find(strcmpi(channel_labels, x), 1), nm, 'UniformOutput', false);
ch_theta1 = cell2mat(nameToIdx(names_theta1));
ch_theta2 = cell2mat(nameToIdx(names_theta2));
ch_alpha  = cell2mat(nameToIdx(names_alpha));

%% 3. Average Epoch-Wise Power Computation per ROI
% Compute spatial average across selected channels for each epoch

% Theta1 (Fronto-Central)
theta1_p = zeros(height(T), 1);
for ch = ch_theta1
    theta1_p = theta1_p + getBand(ch, 'Theta1') / numel(ch_theta1);
end

% Theta2 (Parieto-Occipital)
theta2_p = zeros(height(T), 1);
for ch = ch_theta2
    theta2_p = theta2_p + getBand(ch, 'Theta2') / numel(ch_theta2);
end

% Alpha (Parieto-Occipital)
alpha_p = zeros(height(T), 1);
for ch = ch_alpha
    alpha_p = alpha_p + getBand(ch, 'Alpha') / numel(ch_alpha);
end

%% 4. Minimum Epoch Matching per Subject (Rest vs Task)
% Ensure balanced pair comparisons by truncating to the minimum available epochs
n_min = zeros(n_subj, 1);
for s = 1:n_subj
    nR = sum(subj_num == s & cond == "R");
    nT = sum(subj_num == s & cond == "T");
    n_min(s) = min(nR, nT);
end

%% 5. Intra-Subject Analysis (Wilcoxon Rank-Sum / Mann-Whitney Test per Epoch)
% Preallocate arrays to store subject-level medians for subsequent group analysis
med_R_th1 = zeros(n_subj, 1);  med_T_th1 = zeros(n_subj, 1);
med_R_th2 = zeros(n_subj, 1);  med_T_th2 = zeros(n_subj, 1);
med_R_alp = zeros(n_subj, 1);  med_T_alp = zeros(n_subj, 1);

p_intra_th1 = zeros(n_subj, 1);
p_intra_th2 = zeros(n_subj, 1);
p_intra_alp = zeros(n_subj, 1);

fprintf('========================================================================================\n');
fprintf('     1. INTRA-SUBJECT ANALYSIS (Mann-Whitney / Wilcoxon Rank-Sum per Epoch)\n');
fprintf('========================================================================================\n');

for s = 1:n_subj
    % Extract balanced epochs (1:n_min(s)) for Rest and Task
    r_th1 = theta1_p(subj_num == s & cond == "R");  r_th1 = r_th1(1:n_min(s));
    t_th1 = theta1_p(subj_num == s & cond == "T");  t_th1 = t_th1(1:n_min(s));
    
    r_th2 = theta2_p(subj_num == s & cond == "R");  r_th2 = r_th2(1:n_min(s));
    t_th2 = theta2_p(subj_num == s & cond == "T");  t_th2 = t_th2(1:n_min(s));
    
    r_alp = alpha_p(subj_num == s & cond == "R");   r_alp = r_alp(1:n_min(s));
    t_alp = alpha_p(subj_num == s & cond == "T");   t_alp = t_alp(1:n_min(s));
    
    % Store median band power values
    med_R_th1(s) = median(r_th1);  med_T_th1(s) = median(t_th1);
    med_R_th2(s) = median(r_th2);  med_T_th2(s) = median(t_th2);
    med_R_alp(s) = median(r_alp);  med_T_alp(s) = median(t_alp);
    
    % Execute intra-subject non-parametric rank-sum test across epochs
    p_intra_th1(s) = ranksum(t_th1, r_th1);
    p_intra_th2(s) = ranksum(t_th2, r_th2);
    p_intra_alp(s) = ranksum(t_alp, r_alp);
    
    fprintf('SUBJECT %d (%d epochs compared per condition):\n', s, n_min(s));
    fprintf('  - Theta1 (Fronto-Central)   : Med R = %.4f | Med T = %.4f | p = %.4e\n', ...
        med_R_th1(s), med_T_th1(s), p_intra_th1(s));
    fprintf('  - Theta2 (Parieto-Occipital): Med R = %.4f | Med T = %.4f | p = %.4e\n', ...
        med_R_th2(s), med_T_th2(s), p_intra_th2(s));
    fprintf('  - Alpha  (Parieto-Occipital): Med R = %.4f | Med T = %.4f | p = %.4e\n\n', ...
        med_R_alp(s), med_T_alp(s), p_intra_alp(s));
end

%% 6. Group-Level Analysis (Wilcoxon Signed-Rank Test, Paired N = 6)
fprintf('========================================================================================\n');
fprintf('     2. GROUP-LEVEL STATISTICAL ANALYSIS (Wilcoxon Signed-Rank Test, N = 6)\n');
fprintf('========================================================================================\n');

% Execute paired Wilcoxon signed-rank test on subject median values (Task vs Rest)
[p_grp_th1, ~, stats_th1] = signrank(med_T_th1, med_R_th1);
[p_grp_th2, ~, stats_th2] = signrank(med_T_th2, med_R_th2);
[p_grp_alp, ~, stats_alp] = signrank(med_T_alp, med_R_alp);

diff_th1 = med_T_th1 - med_R_th1;
diff_th2 = med_T_th2 - med_R_th2;
diff_alp = med_T_alp - med_R_alp;

% Display Theta1 results
fprintf('A. THETA1 (4.1 - 5.8 Hz) — FRONTO-CENTRAL ROI (%s):\n', strjoin(names_theta1, ', '));
fprintf('   - Subjects with increase (T > R): %d out of %d\n', sum(diff_th1 > 0), n_subj);
fprintf('   - Signed rank statistic (W)     : %d\n', stats_th1.signedrank);
fprintf('   - p-value                        : %.5f  %s\n\n', p_grp_th1, getSignifTag(p_grp_th1));

% Display Theta2 results
fprintf('B. THETA2 (5.9 - 7.4 Hz) — PARIETO-OCCIPITAL ROI (%s):\n', strjoin(names_theta2, ', '));
fprintf('   - Subjects with decrease/ERD (T < R): %d out of %d\n', sum(diff_th2 < 0), n_subj);
fprintf('   - Signed rank statistic (W)          : %d\n', stats_th2.signedrank);
fprintf('   - p-value                            : %.5f  %s\n\n', p_grp_th2, getSignifTag(p_grp_th2));

% Display Alpha results
fprintf('C. ALPHA (8 - 13 Hz) — PARIETO-OCCIPITAL ROI (%s):\n', strjoin(names_alpha, ', '));
fprintf('   - Subjects with desynchronization (T < R): %d out of %d\n', sum(diff_alp < 0), n_subj);
fprintf('   - Signed rank statistic (W)              : %d\n', stats_alp.signedrank);
fprintf('   - p-value                                : %.5f  %s\n', p_grp_alp, getSignifTag(p_grp_alp));
fprintf('========================================================================================\n');

%% 7. Summary Table Generation and Export
T_group_subbands = table(...
    {'Theta1_FrontoCentral'; 'Theta2_ParietoOccipital'; 'Alpha_ParietoOccipital'}, ...
    [sum(diff_th1 > 0); sum(diff_th2 < 0); sum(diff_alp < 0)], ...
    [stats_th1.signedrank; stats_th2.signedrank; stats_alp.signedrank], ...
    [p_grp_th1; p_grp_th2; p_grp_alp], ...
    string({getSignifTag(p_grp_th1); getSignifTag(p_grp_th2); getSignifTag(p_grp_alp)}), ...
    'VariableNames', {'Band_ROI', 'N_Subjects_ExpectedDirection', 'Statistic_W', 'p_value', 'Outcome'});

disp(T_group_subbands);
writetable(T_group_subbands, fullfile(output_dir, 'Wilcoxon_Group_Subbands.csv'));
fprintf('\nSummary table successfully exported to "%s/Wilcoxon_Group_Subbands.csv"\n', output_dir);

%% Auxiliary Function for Statistical Significance Tags
function tag = getSignifTag(p)
    if p < 0.001
        tag = '--> HIGHLY SIGNIFICANT (*** p < 0.001)';
    elseif p < 0.01
        tag = '--> VERY SIGNIFICANT (** p < 0.01)';
    elseif p < 0.05
        tag = '--> STATISTICALLY SIGNIFICANT (* p < 0.05)';
    else
        tag = '--> NOT SIGNIFICANT (p >= 0.05)';
    end
end
