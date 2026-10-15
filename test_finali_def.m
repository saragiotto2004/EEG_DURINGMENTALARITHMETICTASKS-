clear;
clc;
%L'obiettivo finale è processare i segnali EEG, estrarre delle feature significative
% (caratteristiche fisiologiche) e confrontare le due condizioni per identificare come 
% cambia l'attività cerebrale durante il carico cognitivo.

%ANALISI STATISTICA E CONFRONTO TRA LE CONDIZIONI


% Variazione di Alpha, Theta1 e Theta2 tra Rest (R) e Task (T): rapporto R/T,
% per soggetto e per canale, con indicazione del canale che presenta la massima variazione.
% N.B. per alfa indico il valore più grande di R/T perchè mi dice dove la potenza
% di T è diminuita di più, il contrario per theta1 e theta2.
output_dir='Output';
load(fullfile(output_dir, 'Extracted_Features.mat'), 'T_features');
T=T_features;

% Numero di canali 
vn=T.Properties.VariableNames;
n_ch=sum(startsWith(vn, 'Ch') & endsWith(vn, '_Delta'));

%Condizione: R o T 
cond=extractBefore(upper(string(T.Condition)), 2);
cond=cond(:);

%Divido i dati nei 6 soggetti
isR=cond == "R";
newSubj=[true; isR(2:end) & ~isR(1:end-1)];
subj_num=cumsum(newSubj);
n_subj=max(subj_num);

getBand=@(ch, b) T.(sprintf('Ch%d_%s', ch, b));

%Calcolo e visualizzazione
results=cell(n_subj, 1);
summary=table();

for s=1:n_subj
    iR=(subj_num == s) & (cond == "R");
    iT=(subj_num == s) & (cond == "T");

    Channel   = strings(n_ch, 1);
    Alpha_RT  = nan(n_ch, 1);
    Theta1_RT = nan(n_ch, 1);
    Theta2_RT = nan(n_ch, 1);

    for ch = 1:n_ch
        alpha  = getBand(ch, 'Alpha');
        theta1 = getBand(ch, 'Theta1');
        theta2 = getBand(ch, 'Theta2');

        Channel(ch)   = "Ch" + ch;
        Alpha_RT(ch)  = mean(alpha(iR))  / mean(alpha(iT));    
        Theta1_RT(ch) = mean(theta1(iR)) / mean(theta1(iT));  
        Theta2_RT(ch) = mean(theta2(iR)) / mean(theta2(iT));   
    end

    % Alpha: R/T massimo,Theta1 e Theta2: R/T minimo
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

    summary=[summary; table(s, Channel(iA), Alpha_RT(iA), ...
        Channel(iTh1), Theta1_RT(iTh1), Channel(iTh2), Theta2_RT(iTh2), ...
        'VariableNames', {'Subject','Alpha_CanaleMax','Alpha_RT_Max', ...
                          'Theta1_CanaleMin','Theta1_RT_Min', ...
                          'Theta2_CanaleMin','Theta2_RT_Min'})]; 
end

fprintf('\n---------- RIEPILOGO PER SOGGETTO ----------\n');
disp(summary);

%Salvataggio
xlsx_file=fullfile(output_dir, 'Alpha_Theta_RT.xlsx');
if exist(xlsx_file, 'file'), delete(xlsx_file); end

all_tbl=table();
for s=1:n_subj
    writetable(results{s}, xlsx_file, 'Sheet', sprintf('Soggetto_%d', s));
    tmp=results{s};
    tmp.Subject=repmat(s, height(tmp), 1);
    all_tbl=[all_tbl; movevars(tmp, 'Subject', 'Before', 'Channel')]; 
end
writetable(summary, xlsx_file, 'Sheet', 'Riepilogo');
writetable(all_tbl, fullfile(output_dir, 'Alpha_Theta_RT.csv'));

fprintf('\nSalvato in "%s": Alpha_Theta_RT.xlsx e .csv\n', output_dir);


%VISUALIZZAZIONE RISULTATI
% Andamento della potenza di banda lungo le epoche, R vs T, per soggetto
% Per ogni soggetto le serie R e T sono troncate al minimo numero di epoche
% tra le sue due condizioni.
% Alpha sui canali 11, 10, 13, 19; Theta1 e Theta2 sui canali 10, 9, 18, 19


ch_alpha=[11 10 13 19];
ch_theta1=[14 9 18 19];
ch_theta2=[10, 9, 19, 2, 1, 6];
smooth_w = 1;               

%Potenza per epoca, media sui canali selezionati
theta1_p=zeros(height(T), 1);
for ch=ch_theta1
    theta1_p=theta1_p + T.(sprintf('Ch%d_Theta1', ch)) / numel(ch_theta1)
end
theta2_p=zeros(height(T), 1);
for ch=ch_theta2
    theta2_p=theta2_p + T.(sprintf('Ch%d_Theta2', ch)) / numel(ch_theta2);
end
alpha_p=zeros(height(T), 1);
for ch=ch_alpha
    alpha_p=alpha_p + T.(sprintf('Ch%d_Alpha', ch)) / numel(ch_alpha);
end

pot={theta1_p, theta2_p, alpha_p};
nomi={'Theta1', 'Theta2', 'Alpha'};
canali={ch_theta1, ch_theta2, ch_alpha};
colR=[0.2 0.5 0.8];   
colT=[0.85 0.35 0.2];

%Numero minimo di epoche per ogni soggetto (tra R e T) 
n_min = zeros(n_subj, 1);
for s = 1:n_subj
    nR = sum(subj_num == s & cond == "R");
    nT = sum(subj_num == s & cond == "T");
    n_min(s) = min(nR, nT);
    fprintf('Soggetto %d: R = %d epoche, T = %d epoche -> uso le prime %d\n', s, nR, nT, n_min(s));
end

%Figure
for k=1:3
    figure('Color','w','Position',[100 100 1100 600], 'Name', nomi{k});
    for s=1:n_subj
        subplot(2,3,s); hold on;
        pR = pot{k}(subj_num==s & cond=="R");
        pT = pot{k}(subj_num==s & cond=="T");
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
        ch_str = ['canali ' char(strjoin(string(canali{k}), ', '))];
    end
    sgtitle(sprintf('%s - potenza per epoca, R vs T (%s)', nomi{k}, ch_str));
end



%PROVO AD USARE LA FUNZIONE topoplot PER LOCALIZZARE LE VARIAZIONI DI
%POTENZA NELLE FREQUENZE CONSIDERATE
