clc; clear;
%% Merge metabolites with VIP score
% 读取表格
T1 = readtable('G:\20260715\CTX\CTX_imputed.xlsx');
T2 = readtable('G:\20260715\CTX\\oplsda_vip.csv');

% 匹配代谢物名称
[tf,loc] = ismember(T1.Metabolites,T2.Metabolites);

% 新建VIP1列
T1.VIP1 = nan(height(T1),1);

% 将VIP1赋值给对应代谢物
T1.VIP1(tf) = T2.V1(loc(tf));

% 输出未匹配成功的代谢物
disp('未匹配到VIP1的代谢物数量：')
disp(sum(~tf))

% 保存结果
writetable(T1,'Metabolites_with_VIP.xlsx');

%% 提取数据
% D:I 为Control
Ctrl = table2array(T1(:,4:9));

% J:O 为SD
SD = table2array(T1(:,10:15));

%% 初始化P值
nMetabolites = height(T1);
Pvalue = nan(nMetabolites,1);
MeanCtrl = nan(nMetabolites,1);
MeanSD   = nan(nMetabolites,1);
FC       = nan(nMetabolites,1);
Log2FC   = nan(nMetabolites,1);
%% Student t-test
for i = 1:nMetabolites

    x = Ctrl(i,:);
    y = SD(i,:);

    % 去掉缺失值（如果存在）
    x = x(~isnan(x));
    y = y(~isnan(y));
    % 计算均值
    MeanCtrl(i) = mean(x);
    MeanSD(i)   = mean(y);
    % Fold Change
    if MeanCtrl(i) ~= 0
        FC(i) = MeanSD(i) / MeanCtrl(i);
        Log2FC(i) = log2(FC(i));
    end

 % Student t-test
    if numel(x)>=2 && numel(y)>=2

        [~,p] = ttest2(x,y,'Vartype','equal');

        % 若想使用Welch t test
        % [~,p] = ttest2(x,y,'Vartype','unequal');

        Pvalue(i) = p;

    end

end
%% 添加到表格
T1.P_value  = Pvalue;
T1.MeanCtrl = MeanCtrl;
T1.MeanSD   = MeanSD;
T1.FC        = FC;
T1.Log2FC    = Log2FC;

%% 保存
 writetable(T1,...
    'Metabolite_with_statistics.xlsx');

%%  筛选显著代谢物

sig_idx = (T1.P_value <= 0.05) & ...
          (T1.VIP1 > 1);

SigTable = T1(sig_idx,...
    {'Metabolites','SuperClass','VIP1','P_value','Log2FC'});

%% 保存显著代谢物


writetable(SigTable,...
    'Significant_Metabolites.xlsx');

fprintf('显著代谢物数量：%d\n',height(SigTable));
%% SuperClass统计分析（含上调/下调细分）

% 全部代谢物所属Superclass
AllClass = string(T1.SuperClass);

% 显著变化代谢物所属Superclass及变化方向
SigClass = string(SigTable.SuperClass);
% 假设log2FC > 0 为上调，log2FC < 0 为下调（请根据实际列名调整）
if ismember('Log2FC', SigTable.Properties.VariableNames)
    isUp = SigTable.Log2FC > 0;
else
    error('SigTable 中缺少 log2FC 列，无法判断上/下调。');
end

% 获取所有Superclass名称（基于全部代谢物）
ClassName = unique(AllClass);
nClass = length(ClassName);

% 初始化统计数组
TotalCount   = zeros(nClass,1);   % 该类别代谢物总数
SigCount     = zeros(nClass,1);   % 显著变化总数
SigUpCount   = zeros(nClass,1);   % 显著上调数
SigDownCount = zeros(nClass,1);   % 显著下调数

for i = 1:nClass
    % 当前类别在所有代谢物中的逻辑索引
    inAll = (AllClass == ClassName(i));
    TotalCount(i) = sum(inAll);
    
    % 当前类别在显著代谢物中的逻辑索引
    inSig = (SigClass == ClassName(i));
    SigCount(i) = sum(inSig);
    
    % 细分为上调和下调
    SigUpCount(i)   = sum(inSig & isUp);
    SigDownCount(i) = sum(inSig & ~isUp);
end

% 各类别在显著代谢物中的占比（上调、下调、合计）
Composition_Up   = SigUpCount   / sum(SigCount) * 100;   % 上调占所有显著代谢物的比例
Composition_Down = SigDownCount / sum(SigCount) * 100;   % 下调占所有显著代谢物的比例
Composition_Total = SigCount    / sum(SigCount) * 100;   % 合计（应与 Composition_Up+Composition_Down 一致）

% 各类别在自身总数中的命中率（上调、下调、合计）
HitRatio_Up   = SigUpCount   ./ TotalCount * 100;
HitRatio_Down = SigDownCount ./ TotalCount * 100;
HitRatio_Total = SigCount    ./ TotalCount * 100;        % 即原 HitRatio

%% 汇总结果表格
Result = table(...
    ClassName,...
    TotalCount,...
    SigCount, SigUpCount, SigDownCount,...
    Composition_Total, Composition_Up, Composition_Down,...
    HitRatio_Total, HitRatio_Up, HitRatio_Down,...
    'VariableNames', {...
    'SuperClass',...
    'TotalCount',...
    'SigTotal', 'SigUp', 'SigDown',...
    'Comp_Total_%', 'Comp_Up_%', 'Comp_Down_%',...
    'HitRatio_Total_%', 'HitRatio_Up_%', 'HitRatio_Down_%'});

%% 按 Composition_Total 降序排列（也可按 SigTotal 等排序）
Result = sortrows(Result, 'Comp_Total_%', 'descend');

%% 显示结果
disp(Result)

%% 保存为 Excel 文件
writetable(Result, 'Superclass_Statistics.xlsx');