clear;
clc;
% L'obiettivo finale è processare i segnali EEG, estrarre delle feature significative
% (caratteristiche fisiologiche) e confrontare le due condizioni per identificare come
% cambia l'attività cerebrale durante il carico cognitivo.

% ANALISI STATISTICA E CONFRONTO TRA LE CONDIZIONI
% Variazione di Alpha, Theta1 e Theta2 tra Rest (R) e Task (T): rapporto R/T,
% per soggetto e per canale, con indicazione del canale che presenta la massima variazione.
% N.B. per alfa indico il valore più grande di R/T perché mi dice dove la potenza
% di T è diminuita di più, il contrario per theta1 e theta2.

output_dir='Output';
load(fullfile(output_dir, 'Extracted_Features.mat'), 'T_features', 'channel_labels', 'valid_labels');
T=T_features;

%Canali 
n_ch=numel(channel_labels);
getBand= @(ch, b) T.(sprintf('%s_%s', valid_labels{ch}, b));

%Condizione: R o T
cond=extractBefore(upper(string(T.Condition)), 2);
cond=cond(:);

%Divido i dati nei soggetti 
isR      =cond == "R";
newSubj  =[true; isR(2:end) & ~isR(1:end-1)];
subj_num =cumsum(newSubj);
n_subj   =max(subj_num);
fprintf('Soggetti individuati: %d\n', n_subj);

%Calcolo e visualizzazione
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

    % Alpha: R/T massimo | Theta1 e Theta2: R/T minimo
    [~, iA]   = max(Alpha_RT);
    [~, iTh1] = min(Theta1_RT);
    [~, iTh2] = min(Theta2_RT);

    Alpha_Max  = strings(n_ch, 1);  Alpha_Max(iA)    = "MAX";
    Theta1_Min = strings(n_ch, 1);  Theta1_Min(iTh1) = "MIN";
    Theta2_Min = strings(n_ch, 1);  Theta2_Min(iTh2) = "MIN";

    results{s} = table(Channel, Alpha_RT, Alpha_Max, ...
                       Theta1_RT, Theta1_Min, Theta2_RT, Theta2_Min);

    fprintf('\n---------- SOGGETTO %d ----------\n', s);
    disp(results{s});

    summary = [summary; table(s, Channel(iA), Alpha_RT(iA), ...
        Channel(iTh1), Theta1_RT(iTh1), Channel(iTh2), Theta2_RT(iTh2), ...
        'VariableNames', {'Subject','Alpha_CanaleMax','Alpha_RT_Max', ...
                          'Theta1_CanaleMin','Theta1_RT_Min', ...
                          'Theta2_CanaleMin','Theta2_RT_Min'})];
end

fprintf('\n---------- RIEPILOGO PER SOGGETTO ----------\n');
disp(summary);

% Salvataggio
xlsx_file =fullfile(output_dir, 'Alpha_Theta_RT.xlsx');
if exist(xlsx_file, 'file'), delete(xlsx_file); end

all_tbl=table();
for s = 1:n_subj
    writetable(results{s}, xlsx_file, 'Sheet', sprintf('Soggetto_%d', s));
    tmp=results{s};
    tmp.Subject=repmat(s, height(tmp), 1);
    all_tbl=[all_tbl; movevars(tmp, 'Subject', 'Before', 'Channel')]; 
end
writetable(summary, xlsx_file, 'Sheet', 'Riepilogo');
writetable(all_tbl, fullfile(output_dir, 'Alpha_Theta_RT.csv'));

fprintf('\nSalvato in "%s": Alpha_Theta_RT.xlsx e .csv\n', output_dir);



% Andamento della potenza di banda lungo le epoche, R vs T, per soggetto.
% Per ogni soggetto le serie R e T sono troncate al minimo numero di epoche
% tra le sue due condizioni.
% Canali scelti (da lettratura Seleznov)
%   Alpha  : Pz, P3, P4, O1, O2
%   Theta1 : Fz, F3, F4, Cz
%   Theta2 : gli stessi di alpha

names_alpha={'Pz','P3','P4','O1','O2'};
names_theta1={'Fz','F3','F4','Cz'};
names_theta2={'Pz','P3','P4','O1','O2'};
smooth_w=1; %finestra della media mobile=1 -> non applica smoothing

% Conversione nome -> indice, con controllo che tutti i canali esistano
nameToIdx= @(nm) cellfun(@(x) find(strcmpi(channel_labels, x), 1), nm, 'UniformOutput', false);
chk= @(nm) assert(all(~cellfun(@isempty, nameToIdx(nm))), ...
    'Canale non trovato tra: %s', strjoin(channel_labels, ', '));
chk(names_alpha); chk(names_theta1); chk(names_theta2);

ch_alpha=cell2mat(nameToIdx(names_alpha));
ch_theta1=cell2mat(nameToIdx(names_theta1));
ch_theta2=cell2mat(nameToIdx(names_theta2));

%Potenza per epoca, media sui canali selezionati: per ogni banda, per ogni epoca, fa la media aritmetica
% dei canali scelti
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
nomi   ={'Theta1', 'Theta2', 'Alpha'};%Escono tre vettori, theta1_p, theta2_p, alpha_p, lunghi quanto le righe della tabella.
canali ={names_theta1, names_theta2, names_alpha};
colR=[0.2 0.5 0.8];
colT=[0.85 0.35 0.2];

%Numero minimo di epoche per ogni soggetto (tra R e T) 
n_min=zeros(n_subj, 1);
for s= 1:n_subj
    nR=sum(subj_num == s & cond == "R");
    nT=sum(subj_num == s & cond == "T");
    n_min(s)=min(nR, nT);
    fprintf('Soggetto %d: R = %d epoche, T = %d epoche -> uso le prime %d\n', s, nR, nT, n_min(s));
end

%Figure 
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
        title(sprintf('Soggetto %d (%d epoche)', s, n_min(s)));
        xlabel('Epoca'); ylabel('Potenza (unit^2/Hz)');
        if s == 1, legend('Location', 'best'); end
        grid on; box off;
    end
    sgtitle(sprintf('%s - potenza per epoca, R vs T (canali %s)', ...
        nomi{k}, strjoin(canali{k}, ', ')));
end