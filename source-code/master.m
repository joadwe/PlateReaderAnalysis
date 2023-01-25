%% plate reader analysis script
% Author: Joshua A. Welsh
% Written: 2023-01-25
% Updated: 2023-01-25

% clear previous run
clear; close all; clc

% allow user to select the spreadsheet
[file,path]=uigetfile({'*.xlsx','*.csv'},'File Selector');

filename = fullfile(path,file);
data = readmatrix(filename,'Sheet','Microplate End point','Range','B15:M22'); % import raw sheets
map = readcell(filename,'Sheet','Plate Map'); % import plate map to link samples

if sum(cell2mat(cellfun(@isnumeric, map(1,2:end),'UniformOutput',false))) > 0
    map(1,2:end) = cellfun(@num2str, map(1,2:end),'UniformOutput',false);
end

% define expected plate labels
RowID = {'A','B','C','D','E','F','G','H'};
ColID = strsplit(num2str(1:12));

% create empty array for sample name aggregation
SampleNames = [];
SamplePos = [];

% create empty array for standard name aggregation
StdNames =[];
StdPos = [];

RowPos = [];
ColPos = [];

% go through each name in the platemap to identify where samples are
for i = 1:size(map, 1)

    for ii = 1:size(map,2)

        testName = map{i,ii};

        if ismissing(testName) % check for missing value
        elseif max(strcmp(testName, RowID)) == 1 % check if row header
            RowPos = [RowPos; i ii];
        elseif max(strcmp(testName, ColID)) == 1 % check if column header
            ColPos = [ColPos; i ii];
        elseif contains(testName, {'Std','Blank'}) % check if a standard
            StdNames = [StdNames; {testName}];
            StdPos = [StdPos; i ii];
        else
            SampleNames = [SampleNames; {testName}];
            SamplePos = [SamplePos; i ii];
        end
    end
end

% check row headers are labelled
if ~sum(RowPos(:,2)) == 8
    error('Plate map row header not correctly labelled e.g. A, B, C,...')
end

% check column headers are labelled
if ~sum(ColPos(:,1)) == 12
    error('Plate map column header not correctly labelled e.g. 1, 2, 3,...')
end

% once labelling confirmed reduce plate map labels to samples
PlateMap = map(2:end,2:end);

% extract unique sample names
UniqueNames.Sample = unique(SampleNames);
UniqueNames.Std = unique(StdNames);

% extract standard curve information
StdPos = StdPos-1; % remove row and column headers from index
StdReadings = data(sub2ind(size(data),StdPos(:,1),StdPos(:,2))); % obtain absorbance readings for standard
StdAverage = nan(size(UniqueNames.Std));
StdConc = nan(size(UniqueNames.Std));

% obtain the average absorbance reading for standard
for i = 1:numel(UniqueNames.Std)
    ind = strcmp(UniqueNames.Std{i},StdNames);
    StdAverage(i) = mean(StdReadings(ind));

    if strcmpi(UniqueNames.Std{i}, 'Blank')
        StdConc(i) = 0;
    else
        StdConcStr = regexp(UniqueNames.Std{i},'\d+(\.\d+)?','match');
        StdConc(i) = str2num(StdConcStr{1});
    end
end

% extract sample information
SamplePos = SamplePos-1;
SampleReadings = data(sub2ind(size(data),SamplePos(:,1),SamplePos(:,2))); % obtain absorbance readings for standard
SampleAverage = nan(size(UniqueNames.Sample));

% obtain the average absorbance reading for samples
for i = 1:numel(UniqueNames.Sample)
    ind = strcmp(UniqueNames.Sample{i}, SampleNames);
    SampleAverage(i) = mean(SampleReadings(ind));
end

% linear region of standards
LinRegion = 1:numel(StdConc);

% perform analysis
[Derived_Conc,fig] = BCA_Analysis(StdConc, StdAverage, SampleAverage, LinRegion);

% create pass/fail array
Outcome = repmat({'Passed'},size(Derived_Conc,1),1);
Outcome(Derived_Conc(:,2)==2) = {'Failed'};

% construct export cell array
Export = [UniqueNames.Sample, num2cell(Derived_Conc(:,1)),Outcome];

% sort output into alphabetical order of samples
[~,sortInd] = natsort(UniqueNames.Sample);

% write data to spreadsheet
writecell(Export(sortInd,:), filename, 'Sheet','Concentration')

% export QC data
print(fig, 'Standard Curve.jpeg','-djpeg','-r300')
close(fig)

function [Export_Sample_Conc_Extr, fig] = BCA_Analysis(Conc, Abs, Samples, LinRegion)

% create placeholder for sample concentration export
Export_Sample_Conc_Extr = [];

% set threshold for deviation of linear curve from raw data for
% extrapolation of concentration
ThresholdPercentage = 10;

% set default graphical properties
set(groot,'defaultAxesFontSize',15, 'defaultAxesBox','on',...
    'DefaultLineLineWidth',2,'DefaultAxesLineWidth',2,...
    'DefaultAxesTickDir','out',  'defaultAxesTickDirMode', 'manual')

% sort the concentrations and absorbance readings in order
[Conc, sortInd] = sort(Conc);
Abs = Abs(sortInd);

% obtain absorbance limit of detection from blank
LoD.Abs = [Abs(Conc==0) Abs(Conc==max(Conc))];

% isolate concentrations and absorbance readings above the blank
ConcAboveZero = Conc(Conc>0);
AbsAboveZero = Abs(Conc>0);
Grad = gradient(log10(ConcAboveZero), log10(AbsAboveZero));

% obtained the gradient of the standard curve
for i = 1:numel(ConcAboveZero)
    GradConsInd(i,:) = Grad > Grad(i)*0.9 & Grad < Grad(i)*1.1;
end
GradConsNum = sum(GradConsInd,2);

% check the gradient a region of consistency for linear regression
if max(GradConsNum) < 3
    error('Standard curve has poor linearity')
end

% isolate the standard curve data in the linear region
LinearGradRegion = find(GradConsNum==max(GradConsNum));
ConcLinearRegion = ConcAboveZero(GradConsInd(LinearGradRegion(end),:));
AbsLinearRegion = AbsAboveZero(GradConsInd(LinearGradRegion(end),:));

% perform linear regression
reg1 = robustfit(log10(ConcLinearRegion),log10(AbsLinearRegion));

% compare expected vs acquired concentration
[AcquiredConc] = getConc(reg1, AbsAboveZero);

% obtain LoDs for concentrations and absorbance
LowestLinearConc = 10^(interp1((100.*(AcquiredConc./ConcAboveZero))-100,log10(ConcAboveZero),ThresholdPercentage, 'linear'));
LoD.Conc = [LowestLinearConc max(ConcLinearRegion)];
LoD.Abs(1) = getAbs(reg1, LoD.Conc(1));

% obtain concentations of samples
for i = 1:numel(Samples)
    [Sample_Conc_Extr1(i)] = getConc(reg1, Samples(i));
end

col = [0 0.5 0;
    0.5 0 0];

fig = figure('units','centimeters','position',[0 0 21 29.7],'visible','off');
tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
ax = nexttile;

pt(1) = plot(ax, Conc, Abs, '-o', 'MarkerFaceColor','k','MarkerEdgeColor','k','Color','k');
hold on
xlim(ax, [10^min(floor(log10(Conc(Conc>0)))) 10^max(ceil(log10(Conc(Conc>0))))])
ylim(ax, [10^min(floor(log10(Abs(Abs>0)))) 10^max(ceil(log10(Abs(Abs>0))))])

pt(2) = plot(ax, ConcLinearRegion, AbsLinearRegion, 'or','MarkerFaceColor','r');
pt(3) = plot(ax, getConc(reg1, ax.YLim), ax.YLim,'--b','linewidth',2);

pt(4) = fill([repmat(LoD.Conc(1),1,2), repmat(LoD.Conc(2), 1, 2)],[ax.YLim(1) repmat(LoD.Abs(1), 1, 2), ax.YLim(1)],[0.8 0 0],'FaceAlpha',0.15,'LineStyle','none');
fill([repmat(LoD.Conc(1),1,2), repmat(LoD.Conc(2), 1, 2)],[ax.YLim(2) repmat(LoD.Abs(2), 1, 2), ax.YLim(2)],[0.8 0 0],'FaceAlpha',0.15,'LineStyle','none');
fill([repmat(ax.XLim(1),1,2), repmat(LoD.Conc(1), 1, 2)],[ax.YLim(1) repmat(ax.YLim(2), 1, 2), ax.YLim(1)],[0.8 0 0],'FaceAlpha',0.15,'LineStyle','none');
fill([repmat(ax.XLim(2),1,2), repmat(LoD.Conc(2), 1, 2)],[ax.YLim(1) repmat(ax.YLim(2), 1, 2), ax.YLim(1)],[0.8 0 0],'FaceAlpha',0.15,'LineStyle','none');

pt(5) = fill([repmat(LoD.Conc(1),1,2), repmat(LoD.Conc(2), 1, 2)],[LoD.Abs(1) repmat(LoD.Abs(2), 1, 2), LoD.Abs(1)],[0 0.5 0],'FaceAlpha',0.15,'LineStyle','none');

for i = 1:numel(Samples)
    if Sample_Conc_Extr1(i) > LoD.Conc(1) & Sample_Conc_Extr1(i) < LoD.Conc(2)
        pt(6) = plot([1e-4 Sample_Conc_Extr1(i)],[Samples(i) Samples(i)],'-','Color',col(1,:),'LineStyle','-');
        pt(7) = plot([Sample_Conc_Extr1(i) Sample_Conc_Extr1(i)],[1e-1 Samples(i)],'-','Color',col(1,:),'LineStyle','--');
        Export_Sample_Conc_Extr = [Export_Sample_Conc_Extr; Sample_Conc_Extr1(i),1];
    else
        pt(8) = plot([1e-4 Sample_Conc_Extr1(i)],[Samples(i) Samples(i)],'-','Color',col(2,:),'LineStyle','-');
        pt(9) = plot([Sample_Conc_Extr1(i) Sample_Conc_Extr1(i)],[1e-1 Samples(i)],'-','Color',col(2,:),'LineStyle','--');
        Export_Sample_Conc_Extr = [Export_Sample_Conc_Extr; Sample_Conc_Extr1(i),2];
    end
end

xlabel(ax,'Concentration [µg mL^{-1}]')
ylabel(ax,'Absorbance [a.u.]')
set(ax,'XScale','log','YScale','log')

axis square
grid on
names = {'Raw standards data', 'Raw linear region data',...
    'Standard curve','Non-linear region','Linear Region',...
    'Sample Abs (Passed)','Sample Conc (Passed)',...
    'Sample Abs (Failed)','Sample Conc (Failed)',...
    };
legend(pt, names, 'NumColumns',1,'location','northeastoutside','box','off');

end

function [Conc_Extr1]=getConc(regInfo, Abs)

Conc_Extr1 = 10.^((log10(Abs)-regInfo(1))/regInfo(2));

end

function [Conc_Extr1]=getAbs(regInfo, Abs)

Conc_Extr1 = 10.^(log10(Abs)*regInfo(2)+regInfo(1));

end