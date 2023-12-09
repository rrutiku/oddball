%% Clear the workspace
% ft_defaults
clear all
clc

%% 

PATH    = 'D:\_W\itcf-lab\oddball_comparison\data\data_4_github\\single_trial_data_LearningOddball\\';

AllSubs_P3b = struct();
trNr_P3b = zeros(15,2);
for s = 1:15

    load([PATH, 'pp_', num2str(s), '\\LearningOddball']) 

    % Define outlier trials and remove them
    [idx,d] = nt_find_outlier_trials( nt_trial2mat( good_data.trial ), 2 , 0);  % 2SD
    cfg = [];
    cfg.trials = idx;
    good_data  = ft_redefinetrial(cfg, good_data);
    CNDs = CNDs(idx,:);
    
    % Conditions and averaging 
    cfg = [];
    cfg.demean         = 'yes';                         % baseline correction
    cfg.baselinewindow = [-0.1 0.0];

    cfg.trials     = find(CNDs(:,1) == 1);
    standards      = ft_preprocessing(cfg, good_data);
    trNr_P3b(s,1)  = length(find(CNDs(:,1) == 1));
    
    cfg.trials     = find(CNDs(:,1) == 2); 
    deviants       = ft_preprocessing(cfg, good_data);
    trNr_P3b(s,2)  = length(find(CNDs(:,1) == 2));
    
    cfg = [];
    cfg.keeptrials = 'no';
    AllSubs_P3b.(['s',num2str(s)]).standards = ft_timelockanalysis(cfg, standards);
    AllSubs_P3b.(['s',num2str(s)]).all_dev   = ft_timelockanalysis(cfg, deviants);

end
save([PATH, 'AllSubs_P3b'], 'AllSubs_P3b', 'trNr_P3b');
%%

mean(trNr_P3b(:,1))
std(trNr_P3b(:,1))
[min(trNr_P3b(:,1)) max(trNr_P3b(:,1))]

mean(trNr_P3b(:,2))
std(trNr_P3b(:,2))
[min(trNr_P3b(:,2)) max(trNr_P3b(:,2))]

%% Difference waves

% PATH    = '.........................\\single_trial_data_LearningOddball\\';
% load([PATH, 'AllSubs_P3b'], 'AllSubs_P3b', 'trNr_P3b');

for s = 1:15
    AllSubs_P3b.(['s',num2str(s)]).allDiff     = AllSubs_P3b.(['s',num2str(s)]).all_dev;
    AllSubs_P3b.(['s',num2str(s)]).allDiff.avg = AllSubs_P3b.(['s',num2str(s)]).allDiff.avg - AllSubs_P3b.(['s',num2str(s)]).standards.avg;
end

%%

cfg = [];
cfg.keepindividual = 'no';

GrandP3b = ft_timelockgrandaverage(cfg, AllSubs_P3b.s1.allDiff, AllSubs_P3b.s2.allDiff, AllSubs_P3b.s3.allDiff, ...
    AllSubs_P3b.s4.allDiff, AllSubs_P3b.s5.allDiff, AllSubs_P3b.s6.allDiff, AllSubs_P3b.s7.allDiff, AllSubs_P3b.s8.allDiff, ...
    AllSubs_P3b.s9.allDiff, AllSubs_P3b.s10.allDiff, AllSubs_P3b.s11.allDiff, AllSubs_P3b.s12.allDiff, ...
    AllSubs_P3b.s13.allDiff,AllSubs_P3b.s14.allDiff, AllSubs_P3b.s15.allDiff);
    
GrandStandard = ft_timelockgrandaverage(cfg, AllSubs_P3b.s1.standards, AllSubs_P3b.s2.standards, AllSubs_P3b.s3.standards, ...
    AllSubs_P3b.s4.standards, AllSubs_P3b.s5.standards, AllSubs_P3b.s6.standards, AllSubs_P3b.s7.standards, AllSubs_P3b.s8.standards, ...
    AllSubs_P3b.s9.standards, AllSubs_P3b.s10.standards, AllSubs_P3b.s11.standards, AllSubs_P3b.s12.standards, ...
    AllSubs_P3b.s13.standards, AllSubs_P3b.s14.standards, AllSubs_P3b.s15.standards);

GrandDeviant = ft_timelockgrandaverage(cfg, AllSubs_P3b.s1.all_dev, AllSubs_P3b.s2.all_dev, AllSubs_P3b.s3.all_dev, ...
    AllSubs_P3b.s4.all_dev, AllSubs_P3b.s5.all_dev, AllSubs_P3b.s6.all_dev, AllSubs_P3b.s7.all_dev, AllSubs_P3b.s8.all_dev, ...
    AllSubs_P3b.s9.all_dev, AllSubs_P3b.s10.all_dev, AllSubs_P3b.s11.all_dev, AllSubs_P3b.s12.all_dev, ...
    AllSubs_P3b.s13.all_dev, AllSubs_P3b.s14.all_dev, AllSubs_P3b.s15.all_dev);

%% Plots

% Determine the min-max range
MinMax = [ceil(max(abs([min(min(GrandP3b.avg)) max(max(GrandP3b.avg))])))*-1 ...
          ceil(max(abs([min(min(GrandP3b.avg)) max(max(GrandP3b.avg))])))];
%MinMax = [-1.1 1.4];
    
figure
set(gcf,'Color','w')
subplot_tight(2,1,1);
for c = 1:length(GrandP3b.label)
    plot(GrandP3b.time, GrandP3b.avg(c,:), 'k', 'LineWidth', 1.0); hold all;
end
axis([-0.1,0.7, MinMax])
line([-0.0 -0.0],   get(gca, 'ylim'), 'color', [1.0 0.0 0.0], 'LineWidth', 2.5); hold all;
run('figure_prop.m')
print([PATH, 'P3b_diffWave'],'-dpng','-r300');

%%

cfg = [];
cfg.keepindividual = 'yes';

GrandStandard_wt = ft_timelockgrandaverage(cfg, AllSubs_P3b.s1.standards, AllSubs_P3b.s2.standards, AllSubs_P3b.s3.standards, ...
    AllSubs_P3b.s4.standards, AllSubs_P3b.s5.standards, AllSubs_P3b.s6.standards, AllSubs_P3b.s7.standards, AllSubs_P3b.s8.standards, ...
    AllSubs_P3b.s9.standards, AllSubs_P3b.s10.standards, AllSubs_P3b.s11.standards, AllSubs_P3b.s12.standards, ...
    AllSubs_P3b.s13.standards, AllSubs_P3b.s14.standards, AllSubs_P3b.s15.standards);

GrandDeviant_wt = ft_timelockgrandaverage(cfg, AllSubs_P3b.s1.all_dev, AllSubs_P3b.s2.all_dev, AllSubs_P3b.s3.all_dev, ...
    AllSubs_P3b.s4.all_dev, AllSubs_P3b.s5.all_dev, AllSubs_P3b.s6.all_dev, AllSubs_P3b.s7.all_dev, AllSubs_P3b.s8.all_dev, ...
    AllSubs_P3b.s9.all_dev, AllSubs_P3b.s10.all_dev, AllSubs_P3b.s11.all_dev, AllSubs_P3b.s12.all_dev, ...
    AllSubs_P3b.s13.all_dev, AllSubs_P3b.s14.all_dev, AllSubs_P3b.s15.all_dev);

% Prepare the neighborhood structure
cfg = [];
cfg.method = 'distance';
cfg.neighbourdist = 0.12;
cfg.layout = 'M1_XYZ3_new.sfp';
neighbours = ft_prepare_neighbours(cfg, standards);

% stats for overall MMN
cfg = [];
cfg.channel          = 'all';
cfg.latency          = [0.0 0.7];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.neighbours       = neighbours;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.numrandomization = 10000;

design = [repmat(1,1,15), repmat(2,1,15); 1:15, 1:15];

cfg.design = design;
cfg.ivar   = 1;
cfg.uvar   = 2;

overall_aP3b_stat = ft_timelockstatistics(cfg, GrandDeviant_wt, GrandStandard_wt);
save([PATH, 'overall_aP3b_stat'], 'overall_aP3b_stat', 'cfg');
contourf(overall_aP3b_stat.mask)

%%
load ([PATH, 'overall_aP3b_stat'], 'overall_aP3b_stat');

[overall_aP3b_stat.posclusters.prob]
[overall_aP3b_stat.negclusters.prob]

overall_aP3b_stat.posclusterslabelmat = overall_aP3b_stat.mask .* overall_aP3b_stat.posclusterslabelmat;
overall_aP3b_stat.negclusterslabelmat = overall_aP3b_stat.mask .* overall_aP3b_stat.negclusterslabelmat;

overall_aP3b_stat.negclusterslabelmat(overall_aP3b_stat.negclusterslabelmat == 1) = 0;
overall_aP3b_stat.posclusterslabelmat(1:30,1:260) = 0;

overall_aP3b_stat.mask = overall_aP3b_stat.negclusterslabelmat + overall_aP3b_stat.posclusterslabelmat;
contourf(overall_aP3b_stat.mask)

%% Single channel effects

overall_aP3b_1st_neg_chan = {'Fz'};
overall_aP3b_1st_pos_chan = {'Pz'};

% Fz
fuu = subplot(2,1,1); hold all;
cfg = [];
cfg.channel    = overall_aP3b_1st_neg_chan;
cfg.linestyle  = {'-','-'};
cfg.graphcolor = [0.5 0.5 0.5; 0 0 0];   % [1 0.498 0.314; 0 0.98 0.604];
cfg.linewidth  = 2.5;
cfg.figure     = fuu;
ft_singleplotER(cfg, GrandStandard, GrandDeviant); hold all;

for i = find(overall_aP3b_stat.posclusterslabelmat(find(strcmp(overall_aP3b_stat.label, ...
        overall_aP3b_1st_neg_chan{1,1})),:) > 0)
    plot(overall_aP3b_stat.time(i), -4.3, 'r.', 'MarkerSize',20); hold all;
end
for i = find(overall_aP3b_stat.negclusterslabelmat(find(strcmp(overall_aP3b_stat.label, ...
        overall_aP3b_1st_neg_chan{1,1})),:) > 0)
    plot(overall_aP3b_stat.time(i), -4.3, 'b.', 'MarkerSize',20); hold all;
end
line([-0.0 -0.0],   get(gca, 'ylim'), 'color', [0.0 0.0 0.0], 'LineWidth', 1.0); hold all;
axis([-0.1,0.7, -4.3,1])
title('')
run('figure_prop.m')

% Pz
fuu = subplot(2,1,2); hold all; 
cfg = [];
cfg.channel    = overall_aP3b_1st_pos_chan;
cfg.linestyle  = {'-','-'};
cfg.graphcolor = [0.5 0.5 0.5; 0 0 0];   % [1 0.498 0.314; 0 0.98 0.604];
cfg.linewidth  = 2.5;
cfg.figure     = fuu;
ft_singleplotER(cfg, GrandStandard, GrandDeviant); hold all;

ynow = get(gca, 'ylim');
for i = find(overall_aP3b_stat.posclusterslabelmat(find(strcmp(overall_aP3b_stat.label, ...
        overall_aP3b_1st_pos_chan{1,1})),:) > 0)
    plot(overall_aP3b_stat.time(i), -1, 'r.', 'MarkerSize',20); hold all;
end
for i = find(overall_aP3b_stat.negclusterslabelmat(find(strcmp(overall_aP3b_stat.label, ...
        overall_aP3b_1st_pos_chan{1,1})),:) > 0)
    plot(overall_aP3b_stat.time(i), -1, 'b.', 'MarkerSize',20); hold all;
end
line([-0.0 -0.0],   get(gca, 'ylim'), 'color', [0.0 0.0 0.0], 'LineWidth', 1.0); hold all;
axis([-0.1,0.7, -1,3])
title('')
run('figure_prop.m')

print([PATH, 'aP3b_PosNeg'],'-dpng','-r300');

%% Topoplots post-stimulus

figure
set(gcf,'Color','w')

count = 1;
j = [0.0:0.1:0.4, 0.7];
for k = 1:5
    
    fuu = subplot_tight(3,5,count); count = count+1; hold all;
    
    cfg = [];
    cfg.xlim = [j(k) j(k+1)];
    cfg.zlim = [-3 3];
    cfg.comment = 'no';
    cfg.colorbar = 'no';
    cfg.style = 'straight';
    cfg.layout = 'M1_XYZ3_new.sfp';
    cfg.marker = 'on';
    cfg.markersymbol = '.';
    cfg.figure = fuu;
    
    ft_topoplotER(cfg, GrandP3b); hold all;
end

j = [0.0:0.1:0.4, 0.7];
i = [1:100:401, 701];
for k = 1:5
    
    fuu = subplot_tight(3,5,count); count = count+1; % 4 Rows so it would look like the GlobLoc ones
    
    cfg = [];
    cfg.xlim = [j(k) j(k+1)];
    cfg.comment = 'no';
    cfg.colorbar = 'no';
    cfg.style = 'blank';
    cfg.highlight = {'on', 'on'};
    cfg.highlightchannel = {overall_aP3b_stat.label(sum(overall_aP3b_stat.negclusterslabelmat(:, i(k):i(k+1)) > 0, 2) > 33), ...
                            overall_aP3b_stat.label(sum(overall_aP3b_stat.posclusterslabelmat(:, i(k):i(k+1)) > 0, 2) > 33)};
    cfg.highlightsymbol = {'.', '.'};
    cfg.highlightsize = {45, 45};
    cfg.highlightcolor = {[0 0 1], [1 0 0]};
    cfg.layout = 'M1_XYZ3_new.sfp';
    cfg.marker = 'off';
    cfg.figure = fuu;
   
    ft_topoplotER(cfg, GrandP3b); hold all;
    
end

print([PATH, 'aP3b_diffTopoStat'],'-dpng','-r300');

%% STATS 4 PAPER

chan = 14; %find( sum(overall_aMMN_stat.negclusterslabelmat,2));
linz = mean(GrandP3b.avg(chan,:),1);

mean(linz(402:801),2)
mean(linz(602:801),2)

%%
% FOR PEAKS, MEAN AMPLITUDE, ETC.
%%

cfg = [];
cfg.keepindividual = 'yes';
GrandP3b_wt = ft_timelockgrandaverage(cfg, AllSubs_P3b.s1.allDiff, AllSubs_P3b.s2.allDiff, AllSubs_P3b.s3.allDiff, ...
    AllSubs_P3b.s4.allDiff, AllSubs_P3b.s5.allDiff, AllSubs_P3b.s6.allDiff, AllSubs_P3b.s7.allDiff, AllSubs_P3b.s8.allDiff, ...
    AllSubs_P3b.s9.allDiff, AllSubs_P3b.s10.allDiff, AllSubs_P3b.s11.allDiff, AllSubs_P3b.s12.allDiff, ...
    AllSubs_P3b.s13.allDiff, AllSubs_P3b.s14.allDiff, AllSubs_P3b.s15.allDiff);

fuu = zeros(15,801);
fuuMin = zeros(15,2);
subplot(3,1,1:2); hold all;
for s = 1:15
    load([PATH, 'pp_', num2str(s), '\\aP3b_stat']);
    
    GrandPosCluster = 14;
    
    fuu(s,:) = squeeze(mean(GrandP3b_wt.individual(s,GrandPosCluster,:),2))';
    fuuMin(s,1) = mean(fuu(s,402:601),2);
    fuuMin(s,2) = mean(fuu(s,602:801),2);
    
    if s == 13
        plot(GrandP3b_wt.time, fuu(s,:), 'color', [0.3 0.3 0.3], 'LineWidth', 2.0); hold all;
    else
        plot(GrandP3b_wt.time, fuu(s,:), 'color', [0.5 0.5 0.5], 'LineWidth', 2.0); hold all;
        aP3b_stat.posclusterslabelmat = aP3b_stat.mask .* aP3b_stat.posclusterslabelmat;
        bar = find(sum(aP3b_stat.posclusterslabelmat(:,1:end),1) > 0);
        plot(GrandP3b_wt.time(bar+100), fuu(s,bar+100), 'color', [1.0 0.0 0.0], 'LineWidth', 2.0); hold all;
    end
    
end
axis([-0.1,0.7, -2,9])
line([-0.0 -0.0],   get(gca, 'ylim'), 'color', [1.0 0.0 0.0], 'LineWidth', 0.5); hold off;
run('figure_prop.m')
print([PATH, 'aP3b_15subClusters'],'-dpng','-r300');

%%

mean(fuuMin(:,1))
std(fuuMin(:,1))
[min(fuuMin(:,1)) max(fuuMin(:,1))]

mean(fuuMin(:,2))
std(fuuMin(:,2))
[min(fuuMin(:,2)) max(fuuMin(:,2))]

subplot(2,1,1);
scatter(1:15, fuuMin(:,1)); hold all; 
scatter(1:15, fuuMin(:,2)); hold all; 
axis([0,16, -1,10]); hold all; 
for s = 1:15
    line([s s], [fuuMin(s,2),0], 'color', [0.5 0.5 0.5], 'LineWidth', 1.0); hold all;
    line([s s], [fuuMin(s,1),0], 'color', [0.5 0.5 0.5], 'LineWidth', 1.0); hold all;
end
%run('figure_prop.m')
print([PATH, 'aP3b_indVals'],'-dpng','-r300');

%%

% fuuStat = zeros(15,62,701);
% fuuFlat = zeros(15,701);
% for s = 10
%     load([PATH, 'CF_S', subNr{s}, '_CTR/', 'aP3b_stat']);
%     contourf(aP3b_stat.mask,2)
%     
%     mask = zeros(size(aP3b_stat.mask));
%     mask( aP3b_stat.posclusterslabelmat == 1 ) = 1;
%     %mask( aP3b_stat.posclusterslabelmat == 3 ) = 1;
%     %mask( aP3b_stat.negclusterslabelmat == 2 ) = -1;
%     contourf(mask,2)
% 
%     fuuStat(s,:,:) = mask;
%     fuuFlat(s,:) = sum(mask,1);
%     
% end
% save([PATH, 'indiv_aP3b_contourf'], 'fuuStat');
load([PATH, 'indiv_aP3b_contourf'], 'fuuStat');

mask = zeros(size(overall_aP3b_stat.mask));
mask( overall_aP3b_stat.posclusterslabelmat == 1 ) = 1;
mask( overall_aP3b_stat.negclusterslabelmat == 1 ) = -1;
fuuStat(16,:,:) = mask;
% fuuFlat(find(fuuFlat > 0)) = 1;
% fuuFlat(find(fuuFlat < 0)) = -1;
fuuFlat = squeeze(sum(fuuStat,2));
% contourf(fuuFlat,2);
subplot(3,1,1:2); hold all;
for s = 1:15
    
    Cneg = [find(fuuFlat(s,:) <= -3, 1, 'first') ... % onset & offset
            find(fuuFlat(s,:) <= -3, 1, 'last')];
    if ~isempty(Cneg)
        line([Cneg(1) Cneg(2)], [s,s], 'color', [0 0 1], 'LineWidth', 10.0); hold all;
    end

    Cpos = [find(fuuFlat(s,:) >= 3, 1, 'first') ... % onset & offset
            find(fuuFlat(s,:) >= 3, 1, 'last')];
    if ~isempty(Cpos)
        line([Cpos(1) Cpos(2)], [s,s], 'color', [1 0 0], 'LineWidth', 10.0); hold all;
    end

end
axis([-100,700, -1,16])
line([-0.0 -0.0],   get(gca, 'ylim'), 'color', [0.0 0.0 0.0], 'LineWidth', 1); hold off;
%run('figure_prop.m')
print([PATH, 'aP3b_indStats'],'-dpng','-r300');
