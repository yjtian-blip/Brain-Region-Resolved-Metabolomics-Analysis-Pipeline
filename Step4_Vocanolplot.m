
clear; clc; close all;

%% Read table
filename = 'E:\24h brain state\Metabolomics\filtered metabolites0712\MS2 score 0.7\CTX\Metabolite_with_statistics0920.xlsx';   % 替换为实际文件名
T = readtable(filename);

col_Met  = 'Metabolites';      % Metabolites
col_P    = 'P_value';         % P value
col_FC   = 'Log2FC';         % Log2FC
col_VIP  = 'VIP1';             % VIP value


metabolites = T.(col_Met);
pval        = T.(col_P);
log2FC      = T.(col_FC);
vip         = T.(col_VIP);
neglogP = -log10(pval);


valid = ~(isnan(pval) | isnan(log2FC) | isnan(vip));
% valid = ~(isnan(pval) | isnan(log2FC));
metabolites = metabolites(valid);
pval        = pval(valid);
log2FC      = log2FC(valid);
vip         = vip(valid);


sig_up   = (pval <= 0.05) & (log2FC > 0) & (vip   > 1);
sig_down = (pval <= 0.05) & (log2FC < 0) & (vip   > 1);
nonsig   = ~(sig_up | sig_down);

total_metabolites = length(metabolites);
up_regulated      = sum(sig_up);
down_regulated    = sum(sig_down);

fprintf('Total metabolites: %d\n', total_metabolites);
fprintf('Up-regulated : %d\n', up_regulated);
fprintf('Down-regulated : %d\n', down_regulated);

figure('Position', [100 100 800 600]);
hold on;

% Non significant（Grey）
scatter(log2FC(nonsig), -log10(pval(nonsig)), ...
    40, [0.75 0.75 0.75], 'filled', 'DisplayName', 'Non-significant');

% Up-regulated（Red）
scatter(log2FC(sig_up), -log10(pval(sig_up)), ...
    70, [192 0 0]/255, 'filled','DisplayName', 'Up-regulated');

% Down-regulated（Blue）
scatter(log2FC(sig_down), -log10(pval(sig_down)), ...
    70, [0 0 255]/255, 'filled', 'DisplayName', 'Down-regulated');

% Threshold
xline(0, 'k--', 'LineWidth', 1);
yline(-log10(0.05), 'k--', 'LineWidth', 1);
% Highloght top5

% Axis
xlabel('Log_2 Fold Change');
ylabel('-Log_{10}(P-value)');
title('Volcano Plot');

legend('Location', 'best');
grid on;
box on;
set(gca, 'FontSize', 10);
grid off
xlim([-4 4]);
ylim([0 6]);

%% Save results
save('volcano_stats2.mat', 'total_metabolites', 'up_regulated', 'down_regulated');
fprintf('\nStatistics saved to volcano_stats.mat\n');