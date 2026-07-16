
%% 清空环境
clear; clc; close all;

%% ================= 1. 读取 Excel 文件 =================
% 请修改为你的文件路径和文件名
filename = 'G:\20260715\CTX\Metabolite_with_statistics.xlsx';   % 替换为实际文件名
T = readtable(filename);

% 根据你的列名调整以下变量名（请检查 T.Properties.VariableNames）
col_Met  = 'Metabolites';      % 代谢物名称列
col_P    = 'P_value';         % P值列
col_FC   = 'Log2FC';         % Log2 Fold Change 列
col_VIP  = 'VIP1';             % VIP 值列

% 提取所需数据（若列名不同，请修改）
metabolites = T.(col_Met);
pval        = T.(col_P);
log2FC      = T.(col_FC);
vip         = T.(col_VIP);
neglogP = -log10(pval);

% 移除缺失值（如果有）
valid = ~(isnan(pval) | isnan(log2FC) | isnan(vip));
% valid = ~(isnan(pval) | isnan(log2FC));
metabolites = metabolites(valid);
pval        = pval(valid);
log2FC      = log2FC(valid);
vip         = vip(valid);

%% ================= 2. 显著性分类 =================
sig_up   = (pval <= 0.05) & (log2FC > 0) & (vip   > 1);
sig_down = (pval <= 0.05) & (log2FC < 0) & (vip   > 1);
nonsig   = ~(sig_up | sig_down);

% 统计数量
total_metabolites = length(metabolites);
up_regulated      = sum(sig_up);
down_regulated    = sum(sig_down);

fprintf('Total metabolites: %d\n', total_metabolites);
fprintf('Up-regulated : %d\n', up_regulated);
fprintf('Down-regulated : %d\n', down_regulated);

%% ================= 3. 绘制火山图 =================
figure('Position', [100 100 800 600]);
hold on;

% 非显著点（浅灰色）
scatter(log2FC(nonsig), -log10(pval(nonsig)), ...
    40, [0.75 0.75 0.75], 'filled', 'DisplayName', 'Non-significant');

% 显著上调（红色）
scatter(log2FC(sig_up), -log10(pval(sig_up)), ...
    70, [192 0 0]/255, 'filled','DisplayName', 'Up-regulated');

% 显著下调（蓝色）
scatter(log2FC(sig_down), -log10(pval(sig_down)), ...
    70, [0 0 255]/255, 'filled', 'DisplayName', 'Down-regulated');

% 阈值线
xline(0, 'k--', 'LineWidth', 1);
yline(-log10(0.05), 'k--', 'LineWidth', 1);
% Highloght top5

% 坐标轴标签
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

%% ================= 4. 保存统计结果到 .mat 文件 =================
save('volcano_stats2.mat', 'total_metabolites', 'up_regulated', 'down_regulated');
fprintf('\nStatistics saved to volcano_stats.mat\n');