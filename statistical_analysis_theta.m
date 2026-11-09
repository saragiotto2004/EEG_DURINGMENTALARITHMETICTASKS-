
%% TEST THETA 

%% STATISTICAL ANALYSIS SCRIPT THETA (Wilcoxon Signed-Rank and Rank-Sum Tests)
% Based on the exact channel configurations:
% - Theta CP  : ch_theta_cp  = [10 9 19]     (Central-Parietal ROI: C4, C3, Pz)
% - Theta FFp : ch_theta_ffp = [2 1 6]       (Frontal-Frontopolar ROI: Fp2, Fp1, F8)
% - Alpha     : ch_alpha     = [11 10 13 19]  (Centro-Parietal/Posterior ROI: T5, C4, P3, Pz)

clear; clc; close all;

%% 1. Data Loading and Environment Setup
output_dir = 'Output';
mat_file   = fullfile(output_dir, 'Extracted_Features_Theta.mat');

if ~exist(mat_file, 'file')
    error('File not found in "%s". Please run feature extraction first.', mat_file);
end

% Load extracted features table and channel identifiers
load(mat_file, 'T_features', 'valid_labels');
T = T_features;

% Condition extraction: 'R' (Rest) or 'T' (Task)
cond = extractBefore(upper(string(T.Condition)), 2);
cond = cond(:);

% Identify and segment individual subjects based on condition transitions
isR = cond == "R";
newSubj = [true; isR(2:end) & ~isR(1:end-1)];
subj_num = cumsum(newSubj);
n_subj = max(subj_num);

% Anonymous helper function to dynamically fetch band power for a given channel
getBand = @(ch, b) T.(sprintf('%s_%s', valid_labels{ch}, b));

%% 2. Region of Interest (ROI) Channel Definitions
ch_theta_cp  = [10 9 19];      % Central-Parietal (CP: C4, C3, Pz)
ch_theta_ffp = [2 1 6];        % Frontal-Frontopolar (FFp: Fp2, Fp1, F8)
ch_alpha     = [11 10 13 19];   % Centro-Parietal/Posterior (T5, C4, P3, Pz)

%% 3. Average Epoch-Wise Power Computation per ROI
% Compute spatial average across selected channels for each epoch

% Theta CP
theta_cp_p = zeros(height(T), 1);
for ch = ch_theta_cp
    theta_cp_p = theta_cp_p + getBand(ch, 'Theta') / numel(ch_theta_cp);
end

% Theta FFp
theta_ffp_p = zeros(height(T), 1);
for ch = ch_theta_ffp
    theta_ffp_p = theta_ffp_p + getBand(ch, 'Theta') / numel(ch_theta_ffp);
end

% Alpha Posterior
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
med_R_cp  = zeros(n_subj, 1);  med_T_cp  = zeros(n_subj, 1);
med_R_ffp = zeros(n_subj, 1);  med_T_ffp = zeros(n_subj, 1);
med_R_alp = zeros(n_subj, 1);  med_T_alp = zeros(n_subj, 1);

p_intra_cp_all  = zeros(n_subj, 1);
p_intra_ffp_all = zeros(n_subj, 1);
p_intra_alp_all = zeros(n_subj, 1);

fprintf('========================================================================================\n');
fprintf('     1. INTRA-SUBJECT ANALYSIS (Mann-Whitney / Wilcoxon Rank-Sum per Epoch)\n');
fprintf('========================================================================================\n');

for s = 1:n_subj
    % Extract balanced epochs (1:n_min(s)) for Rest and Task
    r_cp  = theta_cp_p(subj_num == s & cond == "R");   r_cp  = r_cp(1:n_min(s));
    t_cp  = theta_cp_p(subj_num == s & cond == "T");   t_cp  = t_cp(1:n_min(s));
    
    r_ffp = theta_ffp_p(subj_num == s & cond == "R");  r_ffp = r_ffp(1:n_min(s));
    t_ffp = theta_ffp_p(subj_num == s & cond == "T");  t_ffp = t_ffp(1:n_min(s));
    
    r_alp = alpha_p(subj_num == s & cond == "R");      r_alp = r_alp(1:n_min(s));
    t_alp = alpha_p(subj_num == s & cond == "T");      t_alp = t_alp(1:n_min(s));
    
    % Store median band power values
    med_R_cp(s)  = median(r_cp);   med_T_cp(s)  = median(t_cp);
    med_R_ffp(s) = median(r_ffp);  med_T_ffp(s) = median(t_ffp);
    med_R_alp(s) = median(r_alp);  med_T_alp(s) = median(t_alp);
    
    % Execute intra-subject non-parametric rank-sum test across epochs
    p_intra_cp_all(s)  = ranksum(t_cp, r_cp);
    p_intra_ffp_all(s) = ranksum(t_ffp, r_ffp);
    p_intra_alp_all(s) = ranksum(t_alp, r_alp);
    
    fprintf('SUBJECT %d (%d epochs compared per condition):\n', s, n_min(s));
    fprintf('  - Theta CP  (Ch %s): Med R = %.4f | Med T = %.4f | p = %.4e\n', ...
        mat2str(ch_theta_cp), med_R_cp(s), med_T_cp(s), p_intra_cp_all(s));
    fprintf('  - Theta FFp (Ch %s): Med R = %.4f | Med T = %.4f | p = %.4e\n', ...
        mat2str(ch_theta_ffp), med_R_ffp(s), med_T_ffp(s), p_intra_ffp_all(s));
    fprintf('  - Alpha     (Ch %s): Med R = %.4f | Med T = %.4f | p = %.4e\n\n', ...
        mat2str(ch_alpha), med_R_alp(s), med_T_alp(s), p_intra_alp_all(s));
end

%% 6. Group-Level Analysis (Wilcoxon Signed-Rank Test, Paired N = 6)
fprintf('========================================================================================\n');
fprintf('     2. GROUP-LEVEL STATISTICAL ANALYSIS (Wilcoxon Signed-Rank Test, N = 6)\n');
fprintf('========================================================================================\n');

% Execute paired Wilcoxon signed-rank test on subject median values (Task vs Rest)
[p_grp_cp,  ~, stats_cp]  = signrank(med_T_cp,  med_R_cp);
[p_grp_ffp, ~, stats_ffp] = signrank(med_T_ffp, med_R_ffp);
[p_grp_alp, ~, stats_alp] = signrank(med_T_alp, med_R_alp);

diff_cp  = med_T_cp  - med_R_cp;
diff_ffp = med_T_ffp - med_R_ffp;
diff_alp = med_T_alp - med_R_alp;

% Display Theta CP results
fprintf('A. UNIFIED THETA — CENTRAL-PARIETAL ROI (Ch: %s):\n', mat2str(ch_theta_cp));
fprintf('   - Channels: %s\n', strjoin(string(valid_labels(ch_theta_cp)), ', '));
fprintf('   - Subjects with increase (T > R): %d out of %d\n', sum(diff_cp > 0), n_subj);
fprintf('   - Signed rank statistic (W)     : %d\n', stats_cp.signedrank);
fprintf('   - p-value                        : %.5f  %s\n\n', ...
    p_grp_cp, getSignifTag(p_grp_cp));

% Display Theta FFp results
fprintf('B. UNIFIED THETA — FRONTAL-FRONTOPOLAR ROI (Ch: %s):\n', mat2str(ch_theta_ffp));
fprintf('   - Channels: %s\n', strjoin(string(valid_labels(ch_theta_ffp)), ', '));
fprintf('   - Subjects with increase (T > R): %d out of %d\n', sum(diff_ffp > 0), n_subj);
fprintf('   - Signed rank statistic (W)     : %d\n', stats_ffp.signedrank);
fprintf('   - p-value                        : %.5f  %s\n\n', ...
    p_grp_ffp, getSignifTag(p_grp_ffp));

% Display Alpha results
fprintf('C. ALPHA — CENTRO-PARIETAL/POSTERIOR ROI (Ch: %s):\n', mat2str(ch_alpha));
fprintf('   - Channels: %s\n', strjoin(string(valid_labels(ch_alpha)), ', '));
fprintf('   - Subjects with desynchronization (T < R): %d out of %d\n', sum(diff_alp < 0), n_subj);
fprintf('   - Signed rank statistic (W)              : %d\n', stats_alp.signedrank);
fprintf('   - p-value                                : %.5f  %s\n', ...
    p_grp_alp, getSignifTag(p_grp_alp));
fprintf('========================================================================================\n');

%% 7. Summary Table Generation and Export
T_group_results = table(...
    {'Theta_CP'; 'Theta_FFp'; 'Alpha_Post'}, ...
    [sum(diff_cp > 0); sum(diff_ffp > 0); sum(diff_alp < 0)], ...
    [stats_cp.signedrank; stats_ffp.signedrank; stats_alp.signedrank], ...
    [p_grp_cp; p_grp_ffp; p_grp_alp], ...
    string({getSignifTag(p_grp_cp); getSignifTag(p_grp_ffp); getSignifTag(p_grp_alp)}), ...
    'VariableNames', {'ROI_Band', 'N_Subjects_ExpectedDirection', 'Statistic_W', 'p_value', 'Outcome'});

disp(T_group_results);
writetable(T_group_results, fullfile(output_dir, 'Wilcoxon_Group_Results.csv'));
fprintf('Summary table successfully exported to "%s/Wilcoxon_Group_Results.csv"\n', output_dir);

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
