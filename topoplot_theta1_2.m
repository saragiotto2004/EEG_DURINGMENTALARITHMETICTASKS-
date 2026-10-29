% Mappe scalpali (topoplot) della potenza di banda: R e poi T
output_dir = 'Output';
load(fullfile(output_dir, 'Extracted_Features.mat'), ...
'T_features', 'channel_labels', 'valid_labels', 'chanlocs');
T    = T_features;
n_ch = numel(channel_labels);
use_dB   = true;
bandList = {'Theta1', 'Theta2', 'Alpha'};


if istable(chanlocs), chanlocs = table2struct(chanlocs); end
if numel(chanlocs) ~= n_ch
    error('chanlocs ha %d canali, i dati %d.', numel(chanlocs), n_ch);
end
has_coords = (isfield(chanlocs, 'theta') && ~all(cellfun(@isempty, {chanlocs.theta}))) || ...
             (isfield(chanlocs, 'X')     && ~all(cellfun(@isempty, {chanlocs.X})));
%Condizione e soggetti
cond     = extractBefore(upper(string(T.Condition)), 2);  cond = cond(:);
isR      = cond == "R";
subj_num = cumsum([true; isR(2:end) & ~isR(1:end-1)]);
n_subj   = max(subj_num);
fprintf('Soggetti individuati: %d\n', n_subj);
%Potenza media per soggetto, condizione, banda, canale
% P(s, c, b, ch): c = 1 (R), 2 (T)
condLab = ["R","T"];
P = nan(n_subj, 2, numel(bandList), n_ch);
for s = 1:n_subj
for c = 1:2
        idx = (subj_num == s) & (cond == condLab(c));
for b = 1:numel(bandList)
for ch = 1:n_ch
                x = T.(sprintf('%s_%s', valid_labels{ch}, bandList{b}));
                P(s,c,b,ch) = mean(x(idx));
end
end
end
end
% Conversione in dB e etichetta dell'unità (nome diverso da "unit")
if use_dB
    P = 10*log10(P);
    unitLabel = 'dB (10 log_{10})';
else
    unitLabel = 'unit^2/Hz';
end
%Figure
condNames = {'R (rest)', 'T (task)'};
for s = 1:n_subj
    figure('Color','w','Name',sprintf('Soggetto %d',s),'Position',[100 50 800 900]);
for b = 1:numel(bandList)
        vals = squeeze(P(s,:,b,:));                  % [2 x n_ch]
        lim  = [min(vals(:)) max(vals(:))];          % scala comune a R e T
for c = 1:2
            subplot(numel(bandList), 2, (b-1)*2 + c);
            topoplot(squeeze(P(s,c,b,:))', chanlocs, 'maplimits', lim, ...
'electrodes', 'labels', 'style', 'both');
            colormap(gca, parula);
            cb = colorbar;
            ylabel(cb, unitLabel);
            title(sprintf('%s - %s', bandList{b}, condNames{c}));
end
end
    sgtitle(sprintf('Soggetto %d - potenza di banda: R e T', s));
end 