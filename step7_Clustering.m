clc; clear; close all;

%% 生成commonly shared代谢物矩阵
regions = {'Hyp','Mid','Pons','CTX'};

AllMet = {};

for i = 1:length(regions)

    filename = sprintf('%s_sigdata.xlsx',regions{i});

    T = readtable(filename,'VariableNamingRule','preserve');

    AllMet = [AllMet; string(T.Metabolites)];

end

% 所有代谢物
[UniqueMet,~,idx] = unique(AllMet);

% 每个代谢物出现次数
Count = accumarray(idx,1);

Result = table(UniqueMet,Count,...
    'VariableNames',{'Metabolites','RegionNumber'});

% 4个脑区共有
Shared4 = Result(Result.RegionNumber==4,:);

% 至少3个脑区共有
Shared3or4 = Result(Result.RegionNumber>=3,:);

writetable(Shared4,...
'Shared_Significant_Metabolites_4Regions1.xlsx');

writetable(Shared3or4,...
'Shared_Significant_Metabolites_in_3or4_Regions1.xlsx');

%%  读取上一步生成的“4个脑区显著代谢物”名单 
% 读取刚刚保存的结果文件
shared_file = 'G:\20260715\sigdata\Shared_Significant_Metabolites_4Regions.xlsx';
SharedT = readtable(shared_file, 'VariableNamingRule', 'preserve');

% 提取代谢物名称，并强制转换为 string 类型以确保匹配精准
SharedMetNames = string(SharedT.Metabolites); 

%% 2. 读取包含 Peak Area 的原始数据总表
 
% 请确认这个路径和文件名是否仍然准确
raw_file = 'G:\20260715\filtered_metabolites_imputed.xlsx';
Raw = readtable(raw_file, 'VariableNamingRule', 'preserve');

% 提取原始表中的代谢物名称（注意：根据你之前的代码，这里原名叫 MS2_name）
Raw_Metabolites = string(Raw. MS2_name);

%% 3. 匹配代谢物并提取对应的 Peak Area 数据
% 在原始表中寻找 SharedMetNames 对应的位置
[tf, loc] = ismember(SharedMetNames, Raw_Metabolites);

% 提取成功匹配的代谢物名称和行索引
matched_loc = loc(tf);
Final_Met_Names = SharedMetNames(tf);

% 提取 Peak Area 矩阵

HeatData = Raw{matched_loc, 9:end};

% 提取对应的样本名称（列名）
SampleNames = Raw.Properties.VariableNames(9:end);

%% 4. 数据标准化处理 (Log2 + Z-score)


HeatData_log = log2(HeatData + 1);


HeatData_z = zscore(HeatData_log, 0, 2); 

%% 5. 结果整理与保存

% 将标准化后的 Z-score 矩阵转换回表格格式，贴上对应的列名
Result_Zscore = array2table(HeatData_z, 'VariableNames', SampleNames);

% 把代谢物名字插回到表格的最前面 (作为第 1 列)
Result_Zscore = addvars(Result_Zscore, Final_Met_Names, 'Before', 1, 'NewVariableNames', 'MS2Name');

% 导出为新的 Excel 表格
save_zscore_filename = 'Shared_Metabolites_Zscore_Data2.xlsx';
writetable(Result_Zscore, save_zscore_filename);

fprintf('==== 处理完成 ====\n');
fprintf('成功提取并标准化了 %d 个共享代谢物的 Peak Area。\n', length(Final_Met_Names));
fprintf('Z-score 矩阵已保存至: %s\n', save_zscore_filename);

%% 6. Clustering Heatmap
SampleNames = string(SampleNames);

CTX = contains(SampleNames,"CTX");
Hyp = contains(SampleNames,"Hyp");
Mid = contains(SampleNames,"Mid");
PA  = contains(SampleNames,"PA");

% newOrder = [find(CTX),find(Hyp),find(Mid),find(PA)];
CTX_Ctrl = find(contains(SampleNames,"CTX") & contains(lower(SampleNames),"ctr"));
CTX_SD   = find(contains(SampleNames,"CTX") & contains(lower(SampleNames),"sd"));

Hyp_Ctrl = find(contains(SampleNames,"Hyp") & contains(lower(SampleNames),"ctr"));
Hyp_SD   = find(contains(SampleNames,"Hyp") & contains(lower(SampleNames),"sd"));

Mid_Ctrl = find(contains(SampleNames,"Mid") & contains(lower(SampleNames),"ctr"));
Mid_SD   = find(contains(SampleNames,"Mid") & contains(lower(SampleNames),"sd"));

PA_Ctrl  = find(contains(SampleNames,"PA") & contains(lower(SampleNames),"ctr"));
PA_SD    = find(contains(SampleNames,"PA") & contains(lower(SampleNames),"sd"));

newOrder = [...
    CTX_Ctrl CTX_SD ...
    Hyp_Ctrl Hyp_SD ...
    Mid_Ctrl Mid_SD ...
    PA_Ctrl  PA_SD];


HeatData_z = HeatData_z(:,newOrder);

SampleNames = SampleNames(newOrder);

%行聚类

D = pdist(HeatData_z,'correlation');

Z = linkage(D,'ward');

leafOrder = optimalleaforder(Z,D);

HeatData_z = HeatData_z(leafOrder,:);

Final_Met_Names = Final_Met_Names(leafOrder);
figure('Position', [150, 150, 1000, 600]);
% imagesc(HeatData_z);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')
nexttile

dendrogram(Z,...
    0,...
    'Orientation','left',...
    'Reorder',leafOrder,...
    'ColorThreshold','default');

set(gca,'XTick',[])
set(gca,'YTick',[])
box off
nexttile

imagesc(HeatData_z)

colormap(redbluecmap)

caxis([-2 2])

colorbar
yticks(1:length(Final_Met_Names))

yticklabels(Final_Met_Names)

xticks(1:length(SampleNames))

xticklabels(SampleNames)

xtickangle(90)

xline(12.5,'k','LineWidth',2)
% 
xline(24.5,'k','LineWidth',2)
% 
xline(36.5,'k','LineWidth',2)
title('Z-score Normalized Peak Area of Shared Metabolites', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Metabolites');
xlabel('Samples');
%% 7. No clustering Heatmap
SampleNames = string(SampleNames);

CTX = contains(SampleNames,"CTX");
Hyp = contains(SampleNames,"Hyp");
Mid = contains(SampleNames,"Mid");
PA  = contains(SampleNames,"PA");

% newOrder = [find(CTX),find(Hyp),find(Mid),find(PA)];
CTX_Ctrl = find(contains(SampleNames,"CTX") & contains(lower(SampleNames),"ctr"));
CTX_SD   = find(contains(SampleNames,"CTX") & contains(lower(SampleNames),"sd"));

Hyp_Ctrl = find(contains(SampleNames,"Hyp") & contains(lower(SampleNames),"ctr"));
Hyp_SD   = find(contains(SampleNames,"Hyp") & contains(lower(SampleNames),"sd"));

Mid_Ctrl = find(contains(SampleNames,"Mid") & contains(lower(SampleNames),"ctr"));
Mid_SD   = find(contains(SampleNames,"Mid") & contains(lower(SampleNames),"sd"));

PA_Ctrl  = find(contains(SampleNames,"PA") & contains(lower(SampleNames),"ctr"));
PA_SD    = find(contains(SampleNames,"PA") & contains(lower(SampleNames),"sd"));

newOrder = [...
    CTX_Ctrl CTX_SD ...
    Hyp_Ctrl Hyp_SD ...
    Mid_Ctrl Mid_SD ...
    PA_Ctrl  PA_SD];


HeatData_z = HeatData_z(:,newOrder);

SampleNames = SampleNames(newOrder);
figure('Position',[150 150 1200 800]);

imagesc(HeatData_z);

colormap(redbluecmap);
caxis([-2 2]);
colorbar;

yticks(1:length(Final_Met_Names));
yticklabels(Final_Met_Names);

xticks(1:length(SampleNames));

xticklabels(SampleNames);
xtickangle(90);

ylabel('Metabolites');
xlabel('Samples');
title('Z-score normalized metabolite abundance');