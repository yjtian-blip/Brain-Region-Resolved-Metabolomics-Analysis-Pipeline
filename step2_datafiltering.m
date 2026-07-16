%%  Metabolomics PCA & Correlation Analysis

clear;clc;close all

%% Read Excel
T = readtable('G:\20260712\Merged_Unique_Metabolites0715.xlsx');
T.Properties.VariableNames'
%% Metabolite names
met_names = T.MS2_name;
% qc_cv=T.QC_RSD;

%% Find QC columns

all_varnames = T.Properties.VariableNames;

qc_idx = contains(all_varnames,'QC');

%% Quantification matrix
X = table2array(T(:,9:end));
sample_names = T.Properties.VariableNames(9:end);
%% QC matrix
qc_idx = startsWith(sample_names,'QC');
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

%% Filtered matrix

X_tmp = X(keep_idx,:);

met_names_filt = met_names(keep_idx);

%% Remove metabolites with >50% zeros

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
%代谢物在对照组和SD组保留大于0.5
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
% presence_idx_tmp = ...%% 4个脑区任一一个脑区满足条件则保留该代谢物
%     Hyp_keep | ...
%     Mid_keep | ...
%     PA_keep  | ...
%     CTX_keep;
presence_idx_tmp = Hyp_keep & Mid_keep & PA_keep & CTX_keep;%%4个脑区都存在
presence_idx = false(size(keep_idx));
presence_idx(keep_idx) = presence_idx_tmp;


%%  Combine all filters
final_idx = ...
    keep_idx & ...
    presence_idx;
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
output_folder = fullfile(pwd,'filtered metabolites0712');

if ~exist(output_folder,'dir')
    mkdir(output_folder);
end


%% Save filtered dataset

output_file = fullfile(output_folder,...
    'filtered_metabolites0715.xlsx');

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

        % 找出非零值
        pos_vals = vals(vals>0);

        if isempty(pos_vals)
            continue
        end

        % 最小非零值
        fill_val = min(pos_vals)/5;

        % 用最小值/5替换0
        vals(vals==0) = fill_val;

        X_imp(i,idx) = vals;

    end

end


%% Save imputed data
T_imputed = T_filtered;

% 第9列以后是定量数据
T_imputed(:,9:end) = array2table( ...
    X_imp,...
    'VariableNames',...
    T_filtered.Properties.VariableNames(9:end));

%% Output file

output_file2 = fullfile(output_folder,...
    'filtered_metabolites_imputed0715.xlsx');

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
    'filtered_metabolites_normalized0715.xlsx');

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
    'filtered_metabolites_log2_0715.xlsx');

writetable(T_log,output_file_log);

fprintf('\nLog2 dataset saved:\n%s\n',...
    output_file_log);


%% PCA分析
X_scaled = zscore(X_log,0,2);
[coeff,score,latent,~,explained] = pca(X_scaled');
sample_names = T_imputed.Properties.VariableNames(9:end);

nSample = length(sample_names);

%%  PCA Plot

figure
hold on
for i = 1:nSample

    s = sample_names{i};

    if startsWith(s,'QC')

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


%%  Pearson Correlation自定义 仅针对 QC 样本的相关性分析

% 找出QC样本的索引
idx_QC = find(startsWith(sample_names, 'QC'));

% 提取QC样本的数据（假设 X_log 已是对数转换/归一化后的矩阵，大小为 代谢物数 × 样本数）
X_QC = X_log(:, idx_QC);

% 计算QC样本间的Pearson相关系数矩阵
R_QC = corr(X_QC);

% 绘图
figure;
imagesc(R_QC);
axis square;
colormap(parula);
colorbar;
caxis([0.8 1]);          % 可根据实际数据范围调整

% 设置坐标轴标签为QC样本名称
xticks(1:length(idx_QC));
yticks(1:length(idx_QC));
xticklabels(sample_names(idx_QC));
yticklabels(sample_names(idx_QC));
xtickangle(90);

title('QC Sample Cross-Correlation');

%%  QC CV Distribution

figure

histogram(qc_cv,50)

xlabel('QC CV (%)')

ylabel('Number of Features')

title('QC CV Distribution')

xline(30,'--r','CV=30%')

box on

