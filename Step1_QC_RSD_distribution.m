%% Plot QC_RSD distribution under both POS and NEG mode

clear; clc; close all;

%%  Read data


NEG = readtable('E:\24h brain state\代谢组学结果\20260716\20260715\NEG_raw.xlsx', ...
    'VariableNamingRule','preserve');

POS = readtable('E:\24h brain state\代谢组学结果\20260716\20260715\POS_raw.xlsx', ...
    'VariableNamingRule','preserve');


%%  NEGATIVE MODE

QC_names = "QC" + compose('%02d',1:6);

QCcols_NEG = false(1,width(NEG));

for i = 1:length(QC_names)
    QCcols_NEG = QCcols_NEG | ...
        strcmp(NEG.Properties.VariableNames, QC_names(i));
end

% Check
fprintf('\n=============================================\n');
fprintf('NEGATIVE MODE\n');
fprintf('=============================================\n');

disp('QC columns detected:');
disp(NEG.Properties.VariableNames(QCcols_NEG));

if sum(QCcols_NEG) ~= 6
    error('NEG: Expected 6 QC samples (QC01-QC06), but %d were found.', ...
        sum(QCcols_NEG));
end


%%  Extract QC data

X_QC_NEG = NEG{:,QCcols_NEG};


%% Calculate QC CV / RSD for every feature


Mean_QC_NEG = mean(X_QC_NEG,2,'omitnan');

SD_QC_NEG = std(X_QC_NEG,0,2,'omitnan');

CV_NEG = SD_QC_NEG ./ Mean_QC_NEG * 100;


% Remove invalid values
validCV_NEG = isfinite(CV_NEG) & Mean_QC_NEG > 0;

CV_NEG_valid = CV_NEG(validCV_NEG);


%%  Print basic statistics

n_CV20_NEG = sum(CV_NEG_valid <= 20);
n_CV30_NEG = sum(CV_NEG_valid <= 30);

pct_CV20_NEG = n_CV20_NEG / length(CV_NEG_valid) * 100;
pct_CV30_NEG = n_CV30_NEG / length(CV_NEG_valid) * 100;

fprintf('\nNegative mode QC CV statistics:\n');
fprintf('Total features       : %d\n', height(NEG));
fprintf('Valid CV features    : %d\n', length(CV_NEG_valid));
fprintf('Mean CV              : %.2f %%\n', mean(CV_NEG_valid));
fprintf('Median CV            : %.2f %%\n', median(CV_NEG_valid));
fprintf('CV <= 20%%            : %d / %d (%.2f %%)\n', ...
    n_CV20_NEG, length(CV_NEG_valid), pct_CV20_NEG);

fprintf('CV <= 30%%            : %d / %d (%.2f %%)\n', ...
    n_CV30_NEG, length(CV_NEG_valid), pct_CV30_NEG);

%% Histogram of QC CV distribution

figure('Color','w', ...
       'Position',[200 150 850 600]);

binWidth = 5;

% Automatically determine reasonable upper limit
xMax_NEG = ceil(max(CV_NEG_valid)/10)*10;

% Bin edges
edges_NEG = 0:binWidth:xMax_NEG;

%% Calculate histogram counts

% CV <= 30%
CV_NEG_low = CV_NEG_valid(CV_NEG_valid <= 30);

% CV > 30%
CV_NEG_high = CV_NEG_valid(CV_NEG_valid > 30);

% Calculate counts using the same bin edges
counts_low = histcounts(CV_NEG_low, edges_NEG);
counts_high = histcounts(CV_NEG_high, edges_NEG);

% Bin centers
binCenters = edges_NEG(1:end-1) + binWidth/2;

%%  Plot CV <= 30%


bar(binCenters, counts_low, 1, ...
    'FaceColor',[192 0 0]/255, ...
    'FaceAlpha',0.5, ...
    'EdgeColor',[192 0 0]/255, ...
    'LineWidth',0.5);

hold on;

%%  Plot CV > 30%


bar(binCenters, counts_high, 1, ...
    'FaceColor',[242 242 242]/255, ...
    'FaceAlpha',1, ...
    'EdgeColor',[166 166 166]/255, ...
    'LineWidth',0.5);

%% 30% QC criterion


xline(30,'--', ...
    'Color',[0.85 0.33 0.10], ...
    'LineWidth',2);

xlabel('QC CV (%)', ...
    'FontSize',14, ...
    'FontWeight','bold');

ylabel('Number of features', ...
    'FontSize',14, ...
    'FontWeight','bold');

title('Negative Mode: QC CV Distribution', ...
    'FontSize',15, ...
    'FontWeight','bold');

set(gca, ...
    'FontSize',12, ...
    'LineWidth',1.2, ...
    'TickDir','in', ...
    'Box','off');

xlim([-5 xMax_NEG]);


text(0.97,0.93, ...
    sprintf('CV \\leq 30%%: %.1f%%', ...
    mean(CV_NEG_valid <= 30)*100), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'FontSize',11);



%% Save
exportgraphics(gcf, ...
    'NEG_QC_CV_Distribution.png', ...
    'Resolution',600);


%%  QC01-QC06 Cross-Correlation under NEG mode

R_NEG = corr(X_QC_NEG, ...
    'Rows','pairwise', ...
    'Type','Pearson');


%%  Plot correlation matrix

figure('Color','w', ...
       'Position',[300 150 700 650]);

imagesc(R_NEG);

axis square;

colormap(parula);
colorbar;

caxis([0 1]);

nQC = 6;

xticks(1:nQC);
yticks(1:nQC);

xticklabels(QC_names);
yticklabels(QC_names);

xtickangle(90);

title('Negative Mode: QC Sample Cross-Correlation', ...
    'FontSize',14, ...
    'FontWeight','bold');

hold on;


%% Grid lines

gridColor = [0.25 0.45 0.95];
gridWidth = 0.7;

for x = 0.5:1:(nQC+0.5)
    xline(x, ...
        'Color',gridColor, ...
        'LineWidth',gridWidth);
end

for y = 0.5:1:(nQC+0.5)
    yline(y, ...
        'Color',gridColor, ...
        'LineWidth',gridWidth);
end


%% Overlay Pearson r

for i = 1:nQC
    for j = 1:nQC
        
        text(j,i,sprintf('%.3f',R_NEG(i,j)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',9, ...
            'Color','w', ...
            'FontWeight','normal');
        
    end
end

xlim([0.5 nQC+0.5]);
ylim([0.5 nQC+0.5]);

set(gca, ...
    'FontSize',11, ...
    'LineWidth',1, ...
    'TickDir','in', ...
    'Box','off');

exportgraphics(gcf, ...
    'NEG_QC_CrossCorrelation.png', ...
    'Resolution',600);

%% ============================================================
% NEGATIVE MODE
% PCA preprocessing and PCA
%
% Workflow:
% QC CV <= 30%
% -> biological presence >= 3/6 per group
% -> group-specific min/5 imputation
% -> median normalization
% -> log2 transformation
% -> feature-wise autoscaling
% -> PCA
% PCA samples:
% 48 biological samples + 6 QC samples

fprintf('\n');
fprintf('====================================================\n');
fprintf('NEGATIVE MODE PCA ANALYSIS\n');
fprintf('====================================================\n');


%%  Define all sample columns


varNames_NEG = NEG.Properties.VariableNames;

SampleCols_NEG = ...
    startsWith(varNames_NEG,'ctr') | ...
    startsWith(varNames_NEG,'SD')  | ...
    startsWith(varNames_NEG,'QC');

sample_names_NEG = string(varNames_NEG(SampleCols_NEG));

X_NEG = NEG{:,SampleCols_NEG};

fprintf('Original feature number = %d\n',size(X_NEG,1));
fprintf('Total samples           = %d\n',size(X_NEG,2));


%% QC CV <= 30%


idx_CV30_NEG = validCV_NEG & (CV_NEG <= 30);

fprintf('Features with CV <= 30%%  = %d\n', ...
    sum(idx_CV30_NEG));

X_tmp_NEG = X_NEG(idx_CV30_NEG,:);


%%  Define biological groups

% IMPORTANT:
% QC samples are NOT used for presence filtering.

Ctr_Hyp_idx_NEG = ...
    startsWith(sample_names_NEG,'ctr') & ...
    contains(sample_names_NEG,'Hyp');

SD_Hyp_idx_NEG = ...
    startsWith(sample_names_NEG,'SD') & ...
    contains(sample_names_NEG,'Hyp');

Ctr_Mid_idx_NEG = ...
    startsWith(sample_names_NEG,'ctr') & ...
    contains(sample_names_NEG,'Mid');

SD_Mid_idx_NEG = ...
    startsWith(sample_names_NEG,'SD') & ...
    contains(sample_names_NEG,'Mid');

Ctr_PA_idx_NEG = ...
    startsWith(sample_names_NEG,'ctr') & ...
    contains(sample_names_NEG,'PA');

SD_PA_idx_NEG = ...
    startsWith(sample_names_NEG,'SD') & ...
    contains(sample_names_NEG,'PA');

Ctr_CTX_idx_NEG = ...
    startsWith(sample_names_NEG,'ctr') & ...
    contains(sample_names_NEG,'CTX');

SD_CTX_idx_NEG = ...
    startsWith(sample_names_NEG,'SD') & ...
    contains(sample_names_NEG,'CTX');


%% Check number of samples in each group

fprintf('\nBiological sample numbers:\n');

fprintf('Ctr-Hyp : %d\n',sum(Ctr_Hyp_idx_NEG));
fprintf('SD-Hyp  : %d\n',sum(SD_Hyp_idx_NEG));

fprintf('Ctr-Mid : %d\n',sum(Ctr_Mid_idx_NEG));
fprintf('SD-Mid  : %d\n',sum(SD_Mid_idx_NEG));

fprintf('Ctr-PA  : %d\n',sum(Ctr_PA_idx_NEG));
fprintf('SD-PA   : %d\n',sum(SD_PA_idx_NEG));

fprintf('Ctr-CTX : %d\n',sum(Ctr_CTX_idx_NEG));
fprintf('SD-CTX  : %d\n',sum(SD_CTX_idx_NEG));


%% Presence filtering

% At least 3 positive detections in EACH biological group


Ctr_Hyp_presence_NEG = ...
    sum(X_tmp_NEG(:,Ctr_Hyp_idx_NEG) > 0,2);

SD_Hyp_presence_NEG = ...
    sum(X_tmp_NEG(:,SD_Hyp_idx_NEG) > 0,2);

Ctr_Mid_presence_NEG = ...
    sum(X_tmp_NEG(:,Ctr_Mid_idx_NEG) > 0,2);

SD_Mid_presence_NEG = ...
    sum(X_tmp_NEG(:,SD_Mid_idx_NEG) > 0,2);

Ctr_PA_presence_NEG = ...
    sum(X_tmp_NEG(:,Ctr_PA_idx_NEG) > 0,2);

SD_PA_presence_NEG = ...
    sum(X_tmp_NEG(:,SD_PA_idx_NEG) > 0,2);

Ctr_CTX_presence_NEG = ...
    sum(X_tmp_NEG(:,Ctr_CTX_idx_NEG) > 0,2);

SD_CTX_presence_NEG = ...
    sum(X_tmp_NEG(:,SD_CTX_idx_NEG) > 0,2);



presence_keep_NEG = ...
    (Ctr_Hyp_presence_NEG >= 3) & ...
    (SD_Hyp_presence_NEG  >= 3) & ...
    (Ctr_Mid_presence_NEG >= 3) & ...
    (SD_Mid_presence_NEG  >= 3) & ...
    (Ctr_PA_presence_NEG  >= 3) & ...
    (SD_PA_presence_NEG   >= 3) & ...
    (Ctr_CTX_presence_NEG >= 3) & ...
    (SD_CTX_presence_NEG  >= 3);


fprintf('\nPresence filtering:\n');
fprintf('Retained features = %d\n', ...
    sum(presence_keep_NEG));


%% 6. Apply presence filter


X_filt_NEG = X_tmp_NEG(presence_keep_NEG,:);

fprintf('\nAfter CV + presence filtering:\n');
fprintf('%d features x %d samples\n', ...
    size(X_filt_NEG,1), ...
    size(X_filt_NEG,2));


%%  Missing-value imputation

% For each biological group:
% zero / NaN -> minimum positive value within that group / 5


X_imp_NEG = X_filt_NEG;

% Biological group indices
group_idx_NEG = { ...
    Ctr_Hyp_idx_NEG, ...
    SD_Hyp_idx_NEG, ...
    Ctr_Mid_idx_NEG, ...
    SD_Mid_idx_NEG, ...
    Ctr_PA_idx_NEG, ...
    SD_PA_idx_NEG, ...
    Ctr_CTX_idx_NEG, ...
    SD_CTX_idx_NEG};

group_names_NEG = { ...
    'Ctr-Hyp','SD-Hyp', ...
    'Ctr-Mid','SD-Mid', ...
    'Ctr-PA','SD-PA', ...
    'Ctr-CTX','SD-CTX'};

for r = 1:length(group_idx_NEG)

    idx = group_idx_NEG{r};

    for i = 1:size(X_imp_NEG,1)

        vals = X_imp_NEG(i,idx);

        % Positive observed values
        pos_vals = vals(isfinite(vals) & vals > 0);

        if isempty(pos_vals)
            continue
        end

        fill_val = min(pos_vals) / 5;

        % Replace zero or NaN/Inf
        missing_idx = ~isfinite(vals) | vals <= 0;

        vals(missing_idx) = fill_val;

        X_imp_NEG(i,idx) = vals;

    end

end


QC_idx_NEG = startsWith(sample_names_NEG,'QC');

for i = 1:size(X_imp_NEG,1)

    vals = X_imp_NEG(i,QC_idx_NEG);

    pos_vals = vals(isfinite(vals) & vals > 0);

    if isempty(pos_vals)
        % Fallback: use minimum positive value
        % among all samples for this feature
        all_vals = X_imp_NEG(i,:);
        all_pos = all_vals(isfinite(all_vals) & all_vals > 0);

        if isempty(all_pos)
            continue
        end

        fill_val = min(all_pos)/5;
    else
        fill_val = min(pos_vals)/5;
    end

    missing_idx = ~isfinite(vals) | vals <= 0;

    vals(missing_idx) = fill_val;

    X_imp_NEG(i,QC_idx_NEG) = vals;

end


%% Check missing values

fprintf('\nRemaining missing values after imputation:\n');

fprintf('NaN/Inf = %d\n', ...
    sum(~isfinite(X_imp_NEG),'all'));

fprintf('Zero    = %d\n', ...
    sum(X_imp_NEG <= 0,'all'));


%%  Median normalization


sample_median_NEG = median(X_imp_NEG,1);

X_norm_NEG = X_imp_NEG ./ sample_median_NEG;


%%  Log2 transformation


X_log_NEG = log2(X_norm_NEG + 1e-6);


%%  Feature-wise autoscaling


X_scaled_NEG = zscore(X_log_NEG,0,2);


%%  PCA


[coeff_NEG, ...
 score_NEG, ...
 latent_NEG, ...
 ~, ...
 explained_NEG, ...
 mu_NEG] = ...
    pca(X_scaled_NEG');




fprintf('\nPCA explained variance:\n');

for k = 1:min(10,length(explained_NEG))

    fprintf('PC%d = %.2f %%\n', ...
        k,explained_NEG(k));

end

fprintf('PC1 + PC2 = %.2f %%\n', ...
    explained_NEG(1)+explained_NEG(2));


%%  PCA plot


figure('Color','w', ...
    'Position',[250 150 950 700]);

hold on;

nSample_NEG = length(sample_names_NEG);

for i = 1:nSample_NEG

    s = sample_names_NEG(i);


    %% QC

    if startsWith(s,'QC')

        scatter(score_NEG(i,1), ...
                score_NEG(i,2), ...
                150, ...
                [0 0 0], ...
                'filled', ...
                'o');


    %% Hyp

    elseif contains(s,'Hyp')

        if startsWith(s,'ctr')

            c = [1.0 0.6 0.6];

        else

            c = [0.8 0.0 0.0];

        end

        scatter(score_NEG(i,1), ...
                score_NEG(i,2), ...
                120, ...
                c, ...
                'o', ...
                'filled');


    %% Mid

    elseif contains(s,'Mid')

        if startsWith(s,'ctr')

            c = [0.6 0.8 1.0];

        else

            c = [0.0 0.2 0.8];

        end

        scatter(score_NEG(i,1), ...
                score_NEG(i,2), ...
                120, ...
                c, ...
                's', ...
                'filled');


    %% PA

    elseif contains(s,'PA')

        if startsWith(s,'ctr')

            c = [0.6 1.0 0.6];

        else

            c = [0.0 0.6 0.0];

        end

        scatter(score_NEG(i,1), ...
                score_NEG(i,2), ...
                120, ...
                c, ...
                '^', ...
                'filled');


    %% CTX

    elseif contains(s,'CTX')

        if startsWith(s,'ctr')

            c = [0.8 0.6 1.0];

        else

            c = [0.45 0.0 0.75];

        end

        scatter(score_NEG(i,1), ...
                score_NEG(i,2), ...
                120, ...
                c, ...
                'd', ...
                'filled');

    end

end


%% Axis labels

xlabel(sprintf('PC1 (%.1f%%)', ...
    explained_NEG(1)), ...
    'FontSize',14);

ylabel(sprintf('PC2 (%.1f%%)', ...
    explained_NEG(2)), ...
    'FontSize',14);

title('Negative Mode: Metabolomics PCA', ...
    'FontSize',15);


%% Formatting

box on;
grid on;

set(gca, ...
    'FontSize',11, ...
    'LineWidth',1.2, ...
    'TickDir','in');



h1 = scatter(nan,nan,120,[1.0 0.6 0.6], ...
    'o','filled');

h2 = scatter(nan,nan,120,[0.8 0.0 0.0], ...
    'o','filled');

h3 = scatter(nan,nan,120,[0.6 0.8 1.0], ...
    's','filled');

h4 = scatter(nan,nan,120,[0.0 0.2 0.8], ...
    's','filled');

h5 = scatter(nan,nan,120,[0.6 1.0 0.6], ...
    '^','filled');

h6 = scatter(nan,nan,120,[0.0 0.6 0.0], ...
    '^','filled');

h7 = scatter(nan,nan,120,[0.8 0.6 1.0], ...
    'd','filled');

h8 = scatter(nan,nan,120,[0.45 0.0 0.75], ...
    'd','filled');

h9 = scatter(nan,nan,150,[0 0 0], ...
    'o','filled');


legend([h1 h2 h3 h4 h5 h6 h7 h8 h9], ...
    {'Ctr-Hyp','SD-Hyp', ...
     'Ctr-Mid','SD-Mid', ...
     'Ctr-PA','SD-PA', ...
     'Ctr-CTX','SD-CTX', ...
     'QC'}, ...
    'Location','eastoutside', ...
    'FontSize',10, ...
    'Box','off');




for i = 1:nSample_NEG

    text(score_NEG(i,1), ...
         score_NEG(i,2), ...
         sample_names_NEG(i), ...
         'FontSize',8, ...
         'Interpreter','none');

end


%%  Save PCA


exportgraphics(gcf, ...
    'NEG_PCA_CV30_Presence3.png', ...
    'Resolution',600);


fprintf('\nNEG PCA finished.\n');
%%  POSITIVE MODE

% Identify QC01-QC06


QCcols_POS = false(1,width(POS));

for i = 1:length(QC_names)
    QCcols_POS = QCcols_POS | ...
        strcmp(POS.Properties.VariableNames, QC_names(i));
end

fprintf('\n=============================================\n');
fprintf('POSITIVE MODE\n');
fprintf('=============================================\n');

disp('QC columns detected:');
disp(POS.Properties.VariableNames(QCcols_POS));

if sum(QCcols_POS) ~= 6
    error('POS: Expected 6 QC samples (QC01-QC06), but %d were found.', ...
        sum(QCcols_POS));
end


%Extract QC data


X_QC_POS = POS{:,QCcols_POS};


%% Calculate QC CV / RSD


Mean_QC_POS = mean(X_QC_POS,2,'omitnan');

SD_QC_POS = std(X_QC_POS,0,2,'omitnan');

CV_POS = SD_QC_POS ./ Mean_QC_POS * 100;

validCV_POS = isfinite(CV_POS) & Mean_QC_POS > 0;

CV_POS_valid = CV_POS(validCV_POS);


%%  Print statistics

n_CV20_POS = sum(CV_POS_valid <= 20);
n_CV30_POS = sum(CV_POS_valid <= 30);

pct_CV20_POS = n_CV20_POS / length(CV_POS_valid) * 100;
pct_CV30_POS = n_CV30_POS / length(CV_POS_valid) * 100;


fprintf('\nPositive mode QC CV statistics:\n');
fprintf('Total features       : %d\n', height(POS));
fprintf('Valid CV features    : %d\n', length(CV_POS_valid));
fprintf('Mean CV              : %.2f %%\n', mean(CV_POS_valid));
fprintf('Median CV            : %.2f %%\n', median(CV_POS_valid));
fprintf('CV <= 20%%            : %d / %d (%.2f %%)\n', ...
    n_CV20_POS, length(CV_POS_valid), pct_CV20_POS);

fprintf('CV <= 30%%            : %d / %d (%.2f %%)\n', ...
    n_CV30_POS, length(CV_POS_valid), pct_CV30_POS);


%% Histogram of QC CV distribution

figure('Color','w', ...
       'Position',[200 150 850 600]);

binWidth = 5;

% Automatically determine reasonable upper limit
xMax_POS = ceil(max(CV_POS_valid)/10)*10;

% Bin edges
edges_POS = 0:binWidth:xMax_POS;

%%  Calculate histogram counts


% CV <= 30%
CV_POS_low = CV_POS_valid(CV_POS_valid <= 30);

% CV > 30%
CV_POS_high = CV_POS_valid(CV_POS_valid > 30);

% Calculate counts using the same bin edges
counts_low = histcounts(CV_POS_low, edges_POS);
counts_high = histcounts(CV_POS_high, edges_POS);

% Bin centers
binCenters = edges_POS(1:end-1) + binWidth/2;

%% Plot CV <= 30%

bar(binCenters, counts_low, 1, ...
    'FaceColor',[192 0 0]/255, ...
    'FaceAlpha',0.5, ...
    'EdgeColor',[192 0 0]/255, ...
    'LineWidth',0.5);

hold on;

%%  Plot CV > 30%


bar(binCenters, counts_high, 1, ...
    'FaceColor',[242 242 242]/255, ...
    'FaceAlpha',1, ...
    'EdgeColor',[166 166 166]/255, ...
    'LineWidth',0.5);

%% 30% QC criterion


xline(30,'--', ...
    'Color',[0.85 0.33 0.10], ...
    'LineWidth',2);

%Labels


xlabel('QC CV (%)', ...
    'FontSize',14, ...
    'FontWeight','bold');

ylabel('Number of features', ...
    'FontSize',14, ...
    'FontWeight','bold');

title('Positive Mode: QC CV Distribution', ...
    'FontSize',15, ...
    'FontWeight','bold');

set(gca, ...
    'FontSize',12, ...
    'LineWidth',1.2, ...
    'TickDir','in', ...
    'Box','off');

xlim([-5 xMax_POS]);

%%  Annotation

text(0.97,0.93, ...
    sprintf('CV \\leq 30%%: %.1f%%', ...
    mean(CV_POS_valid <= 30)*100), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'FontSize',11);

%Export


exportgraphics(gcf, ...
    'POS_QC_CV_Distribution.png', ...
    'Resolution',600);

R_POS = corr(X_QC_POS, ...
    'Rows','pairwise', ...
    'Type','Pearson');


%%  Plot correlation matrix


figure('Color','w', ...
       'Position',[300 150 700 650]);

imagesc(R_POS);

axis square;

colormap(parula);
colorbar;

caxis([0 1]);

xticks(1:nQC);
yticks(1:nQC);

xticklabels(QC_names);
yticklabels(QC_names);

xtickangle(90);

title('Positive Mode: QC Sample Cross-Correlation', ...
    'FontSize',14, ...
    'FontWeight','bold');

hold on;


%% Grid lines

for x = 0.5:1:(nQC+0.5)
    xline(x, ...
        'Color',gridColor, ...
        'LineWidth',gridWidth);
end

for y = 0.5:1:(nQC+0.5)
    yline(y, ...
        'Color',gridColor, ...
        'LineWidth',gridWidth);
end


%% Overlay Pearson r

for i = 1:nQC
    for j = 1:nQC
        
        text(j,i,sprintf('%.3f',R_POS(i,j)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',9, ...
            'Color','w', ...
            'FontWeight','normal');
        
    end
end

xlim([0.5 nQC+0.5]);
ylim([0.5 nQC+0.5]);

set(gca, ...
    'FontSize',11, ...
    'LineWidth',1, ...
    'TickDir','in', ...
    'Box','off');

exportgraphics(gcf, ...
    'POS_QC_CrossCorrelation.png', ...
    'Resolution',600);


%% Summary

fprintf('\n=============================================\n');
fprintf('QC ANALYSIS FINISHED\n');
fprintf('=============================================\n');

fprintf('NEG features: %d\n',height(NEG));
fprintf('POS features: %d\n',height(POS));

fprintf('\nMean QC correlation:\n');

% Exclude diagonal
mask = ~eye(6);

fprintf('NEG mean r = %.4f\n',mean(R_NEG(mask)));
fprintf('POS mean r = %.4f\n',mean(R_POS(mask)));

fprintf('\nMedian QC correlation:\n');

fprintf('NEG median r = %.4f\n',median(R_NEG(mask)));
fprintf('POS median r = %.4f\n',median(R_POS(mask)));

fprintf('=============================================\n');
%% PCA preprocessing and PCA

fprintf('\n');
fprintf('====================================================\n');
fprintf('POSITIVE MODE PCA ANALYSIS\n');
fprintf('====================================================\n');


%%  Sample columns


varNames_POS = POS.Properties.VariableNames;

SampleCols_POS = ...
    startsWith(varNames_POS,'ctr') | ...
    startsWith(varNames_POS,'SD')  | ...
    startsWith(varNames_POS,'QC');

sample_names_POS = string(varNames_POS(SampleCols_POS));

X_POS = POS{:,SampleCols_POS};

fprintf('Original feature number = %d\n',size(X_POS,1));
fprintf('Total samples           = %d\n',size(X_POS,2));


%% CV <= 30%

idx_CV30_POS = validCV_POS & (CV_POS <= 30);

fprintf('Features with CV <= 30%% = %d\n', ...
    sum(idx_CV30_POS));

X_tmp_POS = X_POS(idx_CV30_POS,:);


%%  Biological groups

Ctr_Hyp_idx_POS = ...
    startsWith(sample_names_POS,'ctr') & ...
    contains(sample_names_POS,'Hyp');

SD_Hyp_idx_POS = ...
    startsWith(sample_names_POS,'SD') & ...
    contains(sample_names_POS,'Hyp');

Ctr_Mid_idx_POS = ...
    startsWith(sample_names_POS,'ctr') & ...
    contains(sample_names_POS,'Mid');

SD_Mid_idx_POS = ...
    startsWith(sample_names_POS,'SD') & ...
    contains(sample_names_POS,'Mid');

Ctr_PA_idx_POS = ...
    startsWith(sample_names_POS,'ctr') & ...
    contains(sample_names_POS,'PA');

SD_PA_idx_POS = ...
    startsWith(sample_names_POS,'SD') & ...
    contains(sample_names_POS,'PA');

Ctr_CTX_idx_POS = ...
    startsWith(sample_names_POS,'ctr') & ...
    contains(sample_names_POS,'CTX');

SD_CTX_idx_POS = ...
    startsWith(sample_names_POS,'SD') & ...
    contains(sample_names_POS,'CTX');

%% Presence filtering
Ctr_Hyp_presence_POS = ...
    sum(X_tmp_POS(:,Ctr_Hyp_idx_POS) > 0,2);

SD_Hyp_presence_POS = ...
    sum(X_tmp_POS(:,SD_Hyp_idx_POS) > 0,2);

Ctr_Mid_presence_POS = ...
    sum(X_tmp_POS(:,Ctr_Mid_idx_POS) > 0,2);

SD_Mid_presence_POS = ...
    sum(X_tmp_POS(:,SD_Mid_idx_POS) > 0,2);

Ctr_PA_presence_POS = ...
    sum(X_tmp_POS(:,Ctr_PA_idx_POS) > 0,2);

SD_PA_presence_POS = ...
    sum(X_tmp_POS(:,SD_PA_idx_POS) > 0,2);

Ctr_CTX_presence_POS = ...
    sum(X_tmp_POS(:,Ctr_CTX_idx_POS) > 0,2);

SD_CTX_presence_POS = ...
    sum(X_tmp_POS(:,SD_CTX_idx_POS) > 0,2);


%% Require >=3/6 in EVERY biological group


presence_keep_POS = ...
    (Ctr_Hyp_presence_POS >= 3) & ...
    (SD_Hyp_presence_POS  >= 3) & ...
    (Ctr_Mid_presence_POS >= 3) & ...
    (SD_Mid_presence_POS  >= 3) & ...
    (Ctr_PA_presence_POS  >= 3) & ...
    (SD_PA_presence_POS   >= 3) & ...
    (Ctr_CTX_presence_POS >= 3) & ...
    (SD_CTX_presence_POS  >= 3);


fprintf('\nPresence filtering:\n');
fprintf('Retained features = %d\n', ...
    sum(presence_keep_POS));


%%  Apply presence filter


X_filt_POS = X_tmp_POS(presence_keep_POS,:);

fprintf('\nAfter CV + presence filtering:\n');
fprintf('%d features x %d samples\n', ...
    size(X_filt_POS,1), ...
    size(X_filt_POS,2));


%%  Imputation

X_imp_POS = X_filt_POS;

group_idx_POS = { ...
    Ctr_Hyp_idx_POS, ...
    SD_Hyp_idx_POS, ...
    Ctr_Mid_idx_POS, ...
    SD_Mid_idx_POS, ...
    Ctr_PA_idx_POS, ...
    SD_PA_idx_POS, ...
    Ctr_CTX_idx_POS, ...
    SD_CTX_idx_POS};


%% Biological groups

for r = 1:length(group_idx_POS)

    idx = group_idx_POS{r};

    for i = 1:size(X_imp_POS,1)

        vals = X_imp_POS(i,idx);

        pos_vals = vals(isfinite(vals) & vals > 0);

        if isempty(pos_vals)
            continue
        end

        fill_val = min(pos_vals)/5;

        missing_idx = ...
            ~isfinite(vals) | vals <= 0;

        vals(missing_idx) = fill_val;

        X_imp_POS(i,idx) = vals;

    end

end


%% QC samples

QC_idx_POS = startsWith(sample_names_POS,'QC');

for i = 1:size(X_imp_POS,1)

    vals = X_imp_POS(i,QC_idx_POS);

    pos_vals = vals(isfinite(vals) & vals > 0);

    if isempty(pos_vals)

        all_vals = X_imp_POS(i,:);

        all_pos = ...
            all_vals(isfinite(all_vals) & all_vals > 0);

        if isempty(all_pos)
            continue
        end

        fill_val = min(all_pos)/5;

    else

        fill_val = min(pos_vals)/5;

    end

    missing_idx = ...
        ~isfinite(vals) | vals <= 0;

    vals(missing_idx) = fill_val;

    X_imp_POS(i,QC_idx_POS) = vals;

end



fprintf('\nRemaining NaN/Inf = %d\n', ...
    sum(~isfinite(X_imp_POS),'all'));

fprintf('Remaining zero    = %d\n', ...
    sum(X_imp_POS <= 0,'all'));


%%  Median normalization


sample_median_POS = median(X_imp_POS,1);

X_norm_POS = ...
    X_imp_POS ./ sample_median_POS;


%% Log2 transformation


X_log_POS = ...
    log2(X_norm_POS + 1e-6);


%% Autoscaling


X_scaled_POS = ...
    zscore(X_log_POS,0,2);


%%  PCA


[coeff_POS, ...
 score_POS, ...
 latent_POS, ...
 ~, ...
 explained_POS, ...
 mu_POS] = ...
    pca(X_scaled_POS');


%% Explained variance


fprintf('\nPCA explained variance:\n');

for k = 1:min(10,length(explained_POS))

    fprintf('PC%d = %.2f %%\n', ...
        k,explained_POS(k));

end

fprintf('PC1 + PC2 = %.2f %%\n', ...
    explained_POS(1)+explained_POS(2));


%% PCA plot


figure('Color','w', ...
    'Position',[250 150 950 700]);

hold on;

nSample_POS = length(sample_names_POS);

for i = 1:nSample_POS

    s = sample_names_POS(i);


    %% QC

    if startsWith(s,'QC')

        scatter(score_POS(i,1), ...
                score_POS(i,2), ...
                150,[0 0 0], ...
                'filled','o');


    %% Hyp

    elseif contains(s,'Hyp')

        if startsWith(s,'ctr')
            c = [1.0 0.6 0.6];
        else
            c = [0.8 0.0 0.0];
        end

        scatter(score_POS(i,1), ...
                score_POS(i,2), ...
                120,c,'o','filled');


    %% Mid

    elseif contains(s,'Mid')

        if startsWith(s,'ctr')
            c = [0.6 0.8 1.0];
        else
            c = [0.0 0.2 0.8];
        end

        scatter(score_POS(i,1), ...
                score_POS(i,2), ...
                120,c,'s','filled');


    %% PA

    elseif contains(s,'PA')

        if startsWith(s,'ctr')
            c = [0.6 1.0 0.6];
        else
            c = [0.0 0.6 0.0];
        end

        scatter(score_POS(i,1), ...
                score_POS(i,2), ...
                120,c,'^','filled');


    %% CTX

    elseif contains(s,'CTX')

        if startsWith(s,'ctr')
            c = [0.8 0.6 1.0];
        else
            c = [0.45 0.0 0.75];
        end

        scatter(score_POS(i,1), ...
                score_POS(i,2), ...
                120,c,'d','filled');

    end

end


%% Labels

xlabel(sprintf('PC1 (%.1f%%)', ...
    explained_POS(1)), ...
    'FontSize',14);

ylabel(sprintf('PC2 (%.1f%%)', ...
    explained_POS(2)), ...
    'FontSize',14);

title('Positive Mode: Metabolomics PCA', ...
    'FontSize',15);


box on;
grid on;

set(gca, ...
    'FontSize',11, ...
    'LineWidth',1.2, ...
    'TickDir','in');


%% Legend

h1 = scatter(nan,nan,120,[1.0 0.6 0.6], ...
    'o','filled');

h2 = scatter(nan,nan,120,[0.8 0.0 0.0], ...
    'o','filled');

h3 = scatter(nan,nan,120,[0.6 0.8 1.0], ...
    's','filled');

h4 = scatter(nan,nan,120,[0.0 0.2 0.8], ...
    's','filled');

h5 = scatter(nan,nan,120,[0.6 1.0 0.6], ...
    '^','filled');

h6 = scatter(nan,nan,120,[0.0 0.6 0.0], ...
    '^','filled');

h7 = scatter(nan,nan,120,[0.8 0.6 1.0], ...
    'd','filled');

h8 = scatter(nan,nan,120,[0.45 0.0 0.75], ...
    'd','filled');

h9 = scatter(nan,nan,150,[0 0 0], ...
    'o','filled');


legend([h1 h2 h3 h4 h5 h6 h7 h8 h9], ...
    {'Ctr-Hyp','SD-Hyp', ...
     'Ctr-Mid','SD-Mid', ...
     'Ctr-PA','SD-PA', ...
     'Ctr-CTX','SD-CTX', ...
     'QC'}, ...
    'Location','eastoutside', ...
    'FontSize',10, ...
    'Box','off');


%% Sample labels

for i = 1:nSample_POS

    text(score_POS(i,1), ...
         score_POS(i,2), ...
         sample_names_POS(i), ...
         'FontSize',8, ...
         'Interpreter','none');

end


%% Save

exportgraphics(gcf, ...
    'POS_PCA_CV30_Presence3.png', ...
    'Resolution',600);

fprintf('\nPOS PCA finished.\n');