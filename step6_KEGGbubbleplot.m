clc;
clear;
% close all;

%% =========================
% 读取Excel
%% =========================
T = readtable('E:\24h brain state\代谢组学结果\20260704\20260712\Pons\pathway_results_PA.csv');
Pathway = string(T.Pathway);
Impact  = T.Impact;
AdjustP = T.AdjustP;

%% =========================
% 按 AdjustP 显著性排序 (升序)
%% =========================
% 一般情况下 P 值越小越显著，所以按升序排列，取前5名
[AdjustP, idx] = sort(AdjustP, 'descend');
Pathway = Pathway(idx);
Impact  = Impact(idx);

%% =========================
% 提取 Top 5
%% =========================
topN = min(7, length(Pathway)); % 防止数据不足5行报错
Pathway_top = Pathway(1:topN);
Impact_top  = Impact(1:topN);
AdjustP_top = AdjustP(1:topN);

% 将数据翻转，使得排名第 1 (AdjustP 最小) 的通路显示在Y轴最上方
Pathway_top = flip(Pathway_top);
Impact_top  = flip(Impact_top);
AdjustP_top = flip(AdjustP_top);

%% =========================
% Y轴位置计算
%% =========================
Y = 1:topN;

%% =========================
% 气泡大小 (对应 Impact 数值)
%% =========================
% 基础大小为 50，系数可以根据实际 Impact 的范围自行微调
BubbleSize = 50 + 800 * Impact_top; 

%% =========================
% 绘图
%% =========================
figure('Position',[200 100 800 450]);

% X轴: Impact, Y轴: 通路, Size: Impact, Color: AdjustP
% scatter(Impact_top,...
%         Y,...
%         BubbleSize,...
%         AdjustP_top,...
%         'filled',...
%         'MarkerEdgeColor','k');
scatter(AdjustP_top,...
        Y,...
        BubbleSize,...
        AdjustP_top,...
        'filled',...
        'MarkerEdgeColor','k');
%% =========================
% 坐标轴设置
%% =========================
yticks(Y);
yticklabels(Pathway_top);
xlabel('AdjustP'); % 修正了原代码中错误的X轴标签
ylabel('');
set(gca,...
    'FontSize',11,...
    'LineWidth',1.2,...
    'TickLength',[0 0],...
    'YLim',[0.5 topN+0.5]); % 上下增加留白，让气泡显示更完整
% xticks([0, 1, 2, 3, 4]);   % 例如，设置 6 个常用刻度
% xticklabels({'0','1','2','3','4'}); 
%% =========================
% 颜色与图例 (数值越大，颜色越深)
%% =========================
% 自定义浅蓝到深蓝渐变色：
% 第一行是浅色 (对应最小的AdjustP)，最后一行是深色 (对应最大的AdjustP)
% custom_cmap = [linspace(0.85, 0.1, 256)', linspace(0.9, 0.3, 256)', linspace(1, 0.6, 256)'];
% 起点：浅灰 [200 200 200]/255
% 终点：红色   [200  12  12]/255
R = linspace(200/255, 200/255, 256)';   % 恒为 200/255
G = linspace(200/255,  12/255, 256)';   % 从 200 降到 12
B = linspace(200/255,  12/255, 256)';   % 从 200 降到 12
custom_cmap = [R, G, B];
colormap(custom_cmap);

cb = colorbar;
cb.Label.String = 'Adjusted P';

%% =========================
% 颜色范围设置
%% =========================
% 如果Top5的AdjustP有差异，则应用动态范围
if max(AdjustP_top) > min(AdjustP_top)
    caxis([min(AdjustP_top) max(AdjustP_top)]); 
%     注意: 新版 MATLAB 推荐使用 clim([min max]) 替代 caxis
end

%% =========================
% 网格与背景
%% =========================
grid on;
box off;
set(gca,'Color',[0.95 0.95 0.95]);
title('Top 5 KEGG Pathway Enrichment');
