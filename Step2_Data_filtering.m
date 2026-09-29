
%% Metabolomics retained metabolites(%)
clear;clc;close all

%% Read Excel
T = readtable('E:\24h brain state\Metabolomics\20260928\Codes\Merged_Unique_Metabolites0914.xlsx');
T.Properties.VariableNames'
%% Metabolite names
met_names = T.MS2_name;
% qc_cv=T.QC_RSD;

%% Find QC columns

all_varnames = T.Properties.VariableNames;

qc_names = {'QC01','QC02','QC03','QC04','QC05','QC06'};

qc_idx = ismember(all_varnames, qc_names);

%% Quantification matrix
X = table2array(T(:,9:end));
sample_names = T.Properties.VariableNames(9:end);
%% QC matrix
qc_names = {'QC01','QC02','QC03','QC04','QC05','QC06'};
qc_idx = ismember(sample_names, qc_names);
X_QC = X(:,qc_idx);

%% Calculate QC CV

qc_mean = mean(X_QC,2,'omitnan');
qc_std  = std(X_QC,0,2,'omitnan');

qc_cv = qc_std ./ qc_mean * 100;

%% Filter

cv_threshold = 30;

keep_idx = qc_cv < cv_threshold;

fprintf('Total features = %d\n',length(qc_cv));
fprintf('Remain after CV filtering = %d\n',sum(keep_idx));
 %MS2 score cutoff
 MS2_score_threshold = 0.3;

% Read MS2 score
MS2_score = T.MS2_score;
MS2_score_tmp = MS2_score(keep_idx);
% MS2 score >= threshold
MS2_score_idx_tmp = MS2_score_tmp >= MS2_score_threshold;

fprintf('Remain after MS2-score filtering = %d\n', ...
    sum(MS2_score_idx_tmp));


%% Features with QC_RSD ≤ 30% and MS2 score  >= 0.3

MS2_filter_idx = false(size(keep_idx));

MS2_filter_idx(keep_idx) = MS2_score_idx_tmp;
%% Filtered matrix
X_tmp = X(MS2_filter_idx,:);

met_names_filt = met_names(MS2_filter_idx);

%% Remove metabolites with > 50% zeros

Ctr_Hyp_idx = startsWith(sample_names,'ctr') & contains(sample_names,'Hyp');
SD_Hyp_idx  = startsWith(sample_names,'SD')  & contains(sample_names,'Hyp');

Ctr_Mid_idx = startsWith(sample_names,'ctr') & contains(sample_names,'Mid');
SD_Mid_idx  = startsWith(sample_names,'SD')  & contains(sample_names,'Mid');

Ctr_PA_idx  = startsWith(sample_names,'ctr') & contains(sample_names,'PA');
SD_PA_idx   = startsWith(sample_names,'SD')  & contains(sample_names,'PA');

Ctr_CTX_idx = startsWith(sample_names,'ctr') & contains(sample_names,'CTX');
SD_CTX_idx  = startsWith(sample_names,'SD')  & contains(sample_names,'CTX');

%% Put all groups together
Ctr_Hyp_presence = mean(X_tmp (:,Ctr_Hyp_idx)>0,2);
SD_Hyp_presence  = mean(X_tmp (:,SD_Hyp_idx)>0,2);

Ctr_Mid_presence = mean(X_tmp (:,Ctr_Mid_idx)>0,2);
SD_Mid_presence  = mean(X_tmp (:,SD_Mid_idx)>0,2);

Ctr_PA_presence  = mean(X_tmp (:,Ctr_PA_idx)>0,2);
SD_PA_presence   = mean(X_tmp (:,SD_PA_idx)>0,2);

Ctr_CTX_presence = mean(X_tmp (:,Ctr_CTX_idx)>0,2);
SD_CTX_presence  = mean(X_tmp (:,SD_CTX_idx)>0,2);
%Metabolites with presence >= 0.5
Hyp_keep = ...
    (Ctr_Hyp_presence >= 0.5) & ...
    (SD_Hyp_presence  >= 0.5);
Mid_keep = ...
    (Ctr_Mid_presence >= 0.5) & ...
    (SD_Mid_presence  >= 0.5);

PA_keep = ...
    (Ctr_PA_presence >= 0.5) & ...
    (SD_PA_presence  >= 0.5);

CTX_keep = ...
    (Ctr_CTX_presence >= 0.5) & ...
    (SD_CTX_presence  >= 0.5);
% presence_idx_tmp = ...
%     Hyp_keep | ...
%     Mid_keep | ...
%     PA_keep  | ...
%     CTX_keep;
presence_idx_tmp = Hyp_keep & Mid_keep & PA_keep & CTX_keep;
presence_idx = false(size(MS2_filter_idx));

presence_idx(MS2_filter_idx) = presence_idx_tmp;


%%  Combine all filters
% final_idx = ...
%     keep_idx & ...
%     MS2_score_idx & ...
%     presence_idx;
final_idx = presence_idx;
T_filtered = T(final_idx,:);

fprintf('\n');
fprintf('CV filter retained      = %d\n',sum(keep_idx));
fprintf('Presence retained       = %d\n',sum(presence_idx));
fprintf('Final retained features = %d\n',sum(final_idx));


%% Generate filtered table

X_filt = X(final_idx,:);

met_names_filt = met_names(final_idx);

% ms2_score_filt = ms2_score(final_idx);

fprintf('\nFiltered matrix size:\n');
fprintf('%d metabolites × %d samples\n',...
    size(X_filt,1),...
    size(X_filt,2));

%%  Create output folder
output_folder = fullfile(pwd,'filtered metabolites0929');

if ~exist(output_folder,'dir')
    mkdir(output_folder);
end


%% Save filtered dataset

output_file = fullfile(output_folder,...
    'filtered_metabolites0929_0.3.xlsx');

writetable(T_filtered,output_file);

fprintf('\nFiltered dataset saved:\n%s\n',output_file);

%% Group-specific minimum/5 imputation
X_imp = X_filt;
Hyp_idx = contains(sample_names,'Hyp') & ...
          ~startsWith(sample_names,'QC');

Mid_idx = contains(sample_names,'Mid') & ...
          ~startsWith(sample_names,'QC');

PA_idx  = contains(sample_names,'PA') & ...
          ~startsWith(sample_names,'QC');

CTX_idx = contains(sample_names,'CTX') & ...
          ~startsWith(sample_names,'QC');

region_idx = {Hyp_idx,Mid_idx,PA_idx,CTX_idx};

for i = 1:size(X_imp,1)

    for r = 1:length(region_idx)

        idx = region_idx{r};

        vals = X_imp(i,idx);

        % find vals＞0
        pos_vals = vals(vals>0);

        if isempty(pos_vals)
            continue
        end

     
        fill_val = min(pos_vals)/5;

        vals(vals==0) = fill_val;

        X_imp(i,idx) = vals;

    end

end

%% Save imputed data
T_imputed = T_filtered;

T_imputed(:,9:end) = array2table( ...
    X_imp,...
    'VariableNames',...
    T_filtered.Properties.VariableNames(9:end));

%% Output file

output_file2 = fullfile(output_folder,...
    'filtered_metabolites_imputed0929.xlsx');

writetable(T_imputed,output_file2);

fprintf('\nImputed dataset saved:\n%s\n',output_file2);

%% abundance matrix after filtering

X_imp = table2array(T_imputed(:,9:end));

sample_names = T_imputed.Properties.VariableNames(9:end);
met_names_filt = T_imputed.MS2_name;
%% Data normalization
sample_median = median(X_imp,1);
X_norm = X_imp ./ sample_median;
T_norm = T_filtered;

T_norm(:,9:end) = array2table( ...
    X_norm,...
    'VariableNames',...
    T_filtered.Properties.VariableNames(9:end));

output_file_norm = fullfile(output_folder,...
    'filtered_metabolites_normalized0929.xlsx');

writetable(T_norm,output_file_norm);

fprintf('\nNormalized dataset saved:\n%s\n',...
    output_file_norm);
%%Log2转换数据
X_log = log2(X_norm + 1e-6);
%% Save log2 transformed data

T_log = T_filtered;

T_log(:,9:end) = array2table( ...
    X_log,...
    'VariableNames',...
    T_filtered.Properties.VariableNames(9:end));

output_file_log = fullfile(output_folder,...
    'filtered_metabolites_log2_0918.xlsx');

writetable(T_log,output_file_log);

fprintf('\nLog2 dataset saved:\n%s\n',...
    output_file_log);


%% PCA分析
% X_scaled = zscore(X_log,0,2);
% [coeff,score,latent,~,explained] = pca(X_scaled');
% sample_names = T_imputed.Properties.VariableNames(9:end);
% 
% nSample = length(sample_names);
qc_names = {'QC01','QC02','QC03','QC04','QC05','QC06'};

sample_names_all = T_imputed.Properties.VariableNames(9:end);

% Biological samples + QC01-QC06
pca_sample_idx = ~strcmp(sample_names_all,'QC_RSD');

X_PCA = X_log(:,pca_sample_idx);

sample_names = sample_names_all(pca_sample_idx);

X_scaled = zscore(X_PCA,0,2);

[coeff,score,latent,~,explained] = pca(X_scaled');

nSample = length(sample_names);
%%  PCA Plot

figure
hold on
for i = 1:nSample

    s = sample_names{i};

%     if startsWith(s,'QC')
   if ismember(s,qc_names)

        scatter(score(i,1),score(i,2),...
            150,...
            [0 0 0],...
            'filled',...
            'o');

    elseif contains(s,'Hyp')

        if contains(lower(s),'ctr')

            c = [1.0 0.6 0.6];      % light red

        else

            c = [0.8 0.0 0.0];      % dark red

        end

        scatter(score(i,1),score(i,2),...
            120,c,...
            'o','filled');

    elseif contains(s,'Mid')

        if contains(lower(s),'ctr')

            c = [0.6 0.8 1.0];      % light blue

        else

            c = [0.0 0.2 0.8];      % dark blue

        end

        scatter(score(i,1),score(i,2),...
            120,c,...
            's','filled');

    elseif contains(s,'PA')

        if contains(lower(s),'ctr')

            c = [0.6 1.0 0.6];      % light green

        else

            c = [0.0 0.6 0.0];      % dark green

        end

        scatter(score(i,1),score(i,2),...
            120,c,...
            '^','filled');


    elseif contains(s,'CTX')

        if contains(lower(s),'ctr')

            c = [0.8 0.6 1.0];      % light purple

        else

            c = [0.45 0.0 0.75];    % dark purple

        end

        scatter(score(i,1),score(i,2),...
            120,c,...
            'd','filled');

    end

end

xlabel(sprintf('PC1 (%.1f%%)',explained(1)))
ylabel(sprintf('PC2 (%.1f%%)',explained(2)))

title('Metabolomics PCA')

box on
grid on
set(gca,'FontSize',10)

% Legend
hold on
h1 = scatter(nan,nan,120,[1.0 0.6 0.6],'o','filled');
h2 = scatter(nan,nan,120,[0.8 0.0 0.0],'o','filled');

h3 = scatter(nan,nan,120,[0.6 0.8 1.0],'s','filled');
h4 = scatter(nan,nan,120,[0.0 0.2 0.8],'s','filled');

h5 = scatter(nan,nan,120,[0.6 1.0 0.6],'^','filled');
h6 = scatter(nan,nan,120,[0.0 0.6 0.0],'^','filled');

h7 = scatter(nan,nan,120,[0.8 0.6 1.0],'d','filled');
h8 = scatter(nan,nan,120,[0.45 0.0 0.75],'d','filled');

h9 = scatter(nan,nan,120,[0 0 0],'o','filled');

legend([h1 h2 h3 h4 h5 h6 h7 h8 h9],...
{'Ctr-Hyp','SD-Hyp',...
 'Ctr-Mid','SD-Mid',...
 'Ctr-PA','SD-PA',...
 'Ctr-CTX','SD-CTX',...
 'QC'},...
 'Location','eastoutside',...
 'FontSize',10,...
 'Box','off');

% Sample labels
for i = 1:length(sample_names)

    text(score(i,1),score(i,2),...
        sample_names{i},...
        'FontSize',8)

end



%% Pearson Correlation  of QC samples

idx_QC = find(startsWith(sample_names, 'QC'));

X_QC = X_log(:, idx_QC);

R_QC = corr(X_QC);

figure;

imagesc(R_QC);

axis square;
colormap(parula);
colorbar;
caxis([0 1]);

xticks(1:length(idx_QC));
yticks(1:length(idx_QC));

xticklabels(sample_names(idx_QC));
yticklabels(sample_names(idx_QC));

xtickangle(90);

title('QC Sample Cross-Correlation');

nQC = length(idx_QC);

for i = 1:nQC
    for j = 1:nQC
        
      
        text(j, i, sprintf('%.3f', R_QC(i,j)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 8, ...
            'Color', 'w');
        
    end
end

%% MS2 score cutoff threshold analysis


% MS2 scores after QC-CV filtering
MS2_score_plot = T.MS2_score(keep_idx);

% Remove NaN / Inf
MS2_score_plot = MS2_score_plot(isfinite(MS2_score_plot));
MS2_score_plot = MS2_score_plot(MS2_score_plot >= 0.30);

%% Parameters

bin_width = 0.02;
edges = 0:bin_width:1;

% Important cutoffs
cutoffs = [0.30 0.70];

%%  Calculate retained metabolites


n_total = length(MS2_score_plot);

n_retained = zeros(size(cutoffs));
percent_retained = zeros(size(cutoffs));

for i = 1:length(cutoffs)

    n_retained(i) = sum(MS2_score_plot >= cutoffs(i));

    percent_retained(i) = ...
        n_retained(i) / n_total * 100;

end

%% Print results

fprintf('\n========== MS2 Score Distribution ==========\n');

fprintf('Total after QC-CV filtering = %d\n', n_total);

for i = 1:length(cutoffs)

    fprintf('MS2 score >= %.2f : %d metabolites (%.1f%%)\n', ...
        cutoffs(i), ...
        n_retained(i), ...
        percent_retained(i));

end

fprintf('============================================\n');


%% Retention curve

cutoff_all = 0.30:0.01:1.00;

percent_retained_all = zeros(size(cutoff_all));

for i = 1:length(cutoff_all)

    cutoff = cutoff_all(i);

    n_retained_current = sum(MS2_score_plot >= cutoff);

    % Percentage relative to all metabolites entering this analysis
    percent_retained_all(i) = ...
        n_retained_current / n_total * 100;

end


figure('Color','w',...
       'Position',[200 200 1050 450]);


% Panel a: MS2 score distribution

subplot(1,2,1);

hold on;


histogram(MS2_score_plot,...
    'BinEdges',edges,...
    'FaceAlpha',0.65,...
    'EdgeColor','none');


%  Kernel density curve


[y_pdf,x_pdf] = ksdensity(MS2_score_plot);

% Convert probability density to histogram count scale
y_count = y_pdf * n_total * bin_width;

% Only show KDE within the observed score range
valid_kde = x_pdf >= 0.30 & x_pdf <= 1.00;

plot(x_pdf(valid_kde),...
     y_count(valid_kde),...
     'LineWidth',2.0);


%  Cutoff lines


xline(0.30,...
    ':',...
    'LineWidth',1.2);

xline(0.70,...
    ':',...
    'LineWidth',1.2);


%  Cutoff labels


text(0.30,118,...
    'Cutoff = 0.3 (All annotated)',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',10);

text(0.70,118,...
    'Cutoff = 0.7 (high-confidence)',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',10);


% Axes

xlabel('MS2 score',...
    'FontSize',12,...
    'FontWeight','normal');

ylabel('Number of metabolites',...
    'FontSize',12,...
    'FontWeight','normal');

xlim([0 1]);

ylim([0 120]);

xticks(0:0.2:1);

yticks(0:20:120);


%Panel label


text(-0.08,1.05,'a',...
    'Units','normalized',...
    'FontSize',14,...
    'FontWeight','bold');


%  Axis style

set(gca,...
    'FontSize',11,...
    'LineWidth',1.0,...
    'Box','off');

grid off;

hold off;


%%  Panel b: Retention curve


subplot(1,2,2);

hold on;



plot(cutoff_all,...
     percent_retained_all,...
     'LineWidth',1.5);


plot(cutoffs(1),...
     percent_retained(1),...
     'o',...
     'MarkerSize',7,...
     'MarkerFaceColor','w',...
     'LineWidth',1.2);

plot(cutoffs(2),...
     percent_retained(2),...
     'o',...
     'MarkerSize',7,...
     'MarkerFaceColor','w',...
     'LineWidth',1.2);


% Vertical cutoff line at 0.3


xline(0.30,...
    ':',...
    'LineWidth',1.0);


%  Vertical cutoff line at 0.7

xline(0.70,...
    ':',...
    'LineWidth',1.0);


% Horizontal line at retention percentage for cutoff 0.7


yline(percent_retained(2),...
    ':',...
    'LineWidth',1.0);


%  Percentage labels


text(0.30,...
     101,...
     sprintf('%.1f%%',percent_retained(1)),...
     'HorizontalAlignment','center',...
     'VerticalAlignment','bottom',...
     'FontSize',10);

text(0.70,...
     percent_retained(2)+2,...
     sprintf('%.1f%%',percent_retained(2)),...
     'HorizontalAlignment','left',...
     'VerticalAlignment','bottom',...
     'FontSize',10);



xlabel('MS2 score',...
    'FontSize',12,...
    'FontWeight','normal');

ylabel('Retained metabolites (%)',...
    'FontSize',12,...
    'FontWeight','normal');

xlim([0 1]);

ylim([0 100]);

xticks(0:0.1:1);

yticks(0:10:100);



text(-0.08,1.05,'b',...
    'Units','normalized',...
    'FontSize',14,...
    'FontWeight','bold');



set(gca,...
    'FontSize',11,...
    'LineWidth',1.0,...
    'Box','off');

grid off;

hold off;