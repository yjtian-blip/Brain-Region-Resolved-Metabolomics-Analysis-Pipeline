
%%按照RSD&Peak area进行合并
%读取数据
NEG = readtable('G:\20260715\NEG_MS2score_cutoff.xlsx');
POS = readtable('G:\20260715\POS_MS2score_cutoff.xlsx');

NEG.MS2_name = string(NEG.("MS2_name"));
POS.MS2_name = string(POS.("MS2_name"));
%计算QC-RSD
QCcols_NEG = startsWith(NEG.Properties.VariableNames,'QC');
QCcols_POS = startsWith(POS.Properties.VariableNames,'QC');

NEG.QC_RSD = ...
    std(NEG{:,QCcols_NEG},0,2,'omitnan') ./ ...
    mean(NEG{:,QCcols_NEG},2,'omitnan') *100;

POS.QC_RSD = ...
    std(POS{:,QCcols_POS},0,2,'omitnan') ./ ...
    mean(POS{:,QCcols_POS},2,'omitnan') *100;
%计算peak area
SampleCols_NEG = ...
    startsWith(NEG.Properties.VariableNames,'Ctr_') | ...
    startsWith(NEG.Properties.VariableNames,'SD_');

SampleCols_POS = ...
    startsWith(POS.Properties.VariableNames,'Ctr_') | ...
    startsWith(POS.Properties.VariableNames,'SD_');

NEG.MeanArea = ...
    mean(NEG{:,SampleCols_NEG},2,'omitnan');

POS.MeanArea = ...
    mean(POS{:,SampleCols_POS},2,'omitnan');


% 保留同时具有 MS2_name 和 MS2_score 的 feature

NEG.MS2_name = strtrim(string(NEG.MS2_name));
POS.MS2_name = strtrim(string(POS.MS2_name));

idxNEG = ...
    NEG.MS2_name ~= "" & ...
    ~ismissing(NEG.MS2_name) & ...
    ~isnan(NEG.MS2_score);

idxPOS = ...
    POS.MS2_name ~= "" & ...
    ~ismissing(POS.MS2_name) & ...
    ~isnan(POS.MS2_score);

NEG = NEG(idxNEG,:);
POS = POS(idxPOS,:);


%找交集
CommonMet = intersect(NEG.MS2_name,...
                      POS.MS2_name);
%保留最佳模式
Merged = NEG([],:);
ScoreThreshold = 0.10;

for i = 1:length(CommonMet)

    met = CommonMet(i);

    idxN = find(NEG.MS2_name==met,1);
    idxP = find(POS.MS2_name==met,1);

    rsdN = NEG.QC_RSD(idxN);
    rsdP = POS.QC_RSD(idxP);

    scoreN = NEG.MS2_score(idxN);
    scoreP = POS.MS2_score(idxP);

    areaN = NEG.MeanArea(idxN);
    areaP = POS.MeanArea(idxP);

   
    %% 优先比较MS2 score
 
    if abs(scoreN-scoreP) >  ScoreThreshold

        if scoreN > scoreP
            Merged = [Merged; NEG(idxN,:)];
        else
            Merged = [Merged; POS(idxP,:)];
        end

    else
%% 比较RSD
    if rsdN < rsdP

        Merged = [Merged; NEG(idxN,:)];

    elseif rsdP < rsdN

        Merged = [Merged; POS(idxP,:)];

    else

        %% RSD相同比较PeakArea

        if areaN >= areaP

            Merged = [Merged; NEG(idxN,:)];

        else

            Merged = [Merged; POS(idxP,:)];

        end

    end

    end
end
%加入非交集代谢物
NEG_unique = NEG(~ismember(NEG.MS2_name,...
                           CommonMet),:);

POS_unique = POS(~ismember(POS.MS2_name,...
                           CommonMet),:);

Merged = [Merged;
          NEG_unique;
          POS_unique];



%% 排序并保存
Merged = sortrows(Merged,'MS2_name');

%% 删除辅助变量（如果存在）
varsToRemove = {'QC_RSD','MeanArea'};
varsToRemove = intersect(varsToRemove, Merged.Properties.VariableNames);

Merged_out = removevars(Merged, varsToRemove);


writetable(Merged,...
          'Merged_Unique_Metabolites0715.xlsx');