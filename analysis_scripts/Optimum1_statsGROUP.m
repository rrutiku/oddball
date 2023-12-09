%% Clear the workspace
% ft_defaults
clear all
clc

%% 

PATH    = '................................\\single_trial_data_Optimum1\\';

AllSubs_MMN = struct();
trNr_MMN = zeros(15,2);
for s = 1:15

    load([PATH, 'pp_', num2str(s), '\\Optimum1']) 

    % Define outlier trials and remove them
    [idx,d] = nt_find_outlier_trials( nt_trial2mat( good_data.trial ), 2 , 1);  % 2SD
    cfg = [];
    cfg.trials = idx;
    good_data  = ft_redefinetrial(cfg, good_data);
    CNDs = CNDs(idx);
    
    % Conditions and averaging 
    cfg = [];
    cfg.demean         = 'yes';                         % baseline correction
    cfg.baselinewindow = [-0.1 0.0];
    
    cfg.trials     = find(CNDs == 1);
    standards      = ft_preprocessing(cfg, good_data);
    trNr_MMN(s,1)  = length(find(CNDs == 1));

    cfg.trials     = find(ismember(CNDs, [2 3 4 5])); % 
    deviants       = ft_preprocessing(cfg, good_data);
    trNr_MMN(s,2)  = length(find(ismember(CNDs, [2 3 4 5]))); % 
    
    cfg = [];
    cfg.keeptrials = 'no';
    AllSubs_MMN.(['s',num2str(s)]).standards = ft_timelockanalysis(cfg, standards);
    AllSubs_MMN.(['s',num2str(s)]).all_dev   = ft_timelockanalysis(cfg, deviants);

end
save([PATH, 'AllSubs_MMN'], 'AllSubs_MMN', 'trNr_MMN');
%%

mean(trNr_MMN(:,1))
std(trNr_MMN(:,1))
[min(trNr_MMN(:,1)) max(trNr_MMN(:,1))]

mean(trNr_MMN(:,2))
std(trNr_MMN(:,2))
[min(trNr_MMN(:,2)) max(trNr_MMN(:,2))]

fuu = trNr_MMN(:,1) - trNr_MMN(:,2);
mean(fuu)
std(fuu)
[min(fuu) max(fuu)]

%% Difference waves

% PATH    = '..............................\\single_trial_data_Optimum1\\';
% load([PATH, 'AllSubs_MMN']);

for s = 1:15
    AllSubs_MMN.(['s',num2str(s)]).allDiff     = AllSubs_MMN.(['s',num2str(s)]).all_dev;
    AllSubs_MMN.(['s',num2str(s)]).allDiff.avg = AllSubs_MMN.(['s',num2str(s)]).allDiff.avg - AllSubs_MMN.(['s',num2str(s)]).standards.avg;
end

%%

cfg = [];
cfg.keepindividual = 'no';

GrandMMN = ft_timelockgrandaverage(cfg, AllSubs_MMN.s1.allDiff, AllSubs_MMN.s2.allDiff, AllSubs_MMN.s3.allDiff, ...
    AllSubs_MMN.s4.allDiff, AllSubs_MMN.s5.allDiff, AllSubs_MMN.s6.allDiff, AllSubs_MMN.s7.allDiff, AllSubs_MMN.s8.allDiff, ...
    AllSubs_MMN.s9.allDiff, AllSubs_MMN.s10.allDiff, AllSubs_MMN.s11.allDiff, AllSubs_MMN.s12.allDiff, ...
    AllSubs_MMN.s13.allDiff,AllSubs_MMN.s14.allDiff, AllSubs_MMN.s15.allDiff);
    
GrandStandard = ft_timelockgrandaverage(cfg, AllSubs_MMN.s1.standards, AllSubs_MMN.s2.standards, AllSubs_MMN.s3.standards, ...
    AllSubs_MMN.s4.standards, AllSubs_MMN.s5.standards, AllSubs_MMN.s6.standards, AllSubs_MMN.s7.standards, AllSubs_MMN.s8.standards, ...
    AllSubs_MMN.s9.standards, AllSubs_MMN.s10.standards, AllSubs_MMN.s11.standards, AllSubs_MMN.s12.standards, ...
    AllSubs_MMN.s13.standards, AllSubs_MMN.s14.standards, AllSubs_MMN.s15.standards);

GrandDeviant = ft_timelockgrandaverage(cfg, AllSubs_MMN.s1.all_dev, AllSubs_MMN.s2.all_dev, AllSubs_MMN.s3.all_dev, ...
    AllSubs_MMN.s4.all_dev, AllSubs_MMN.s5.all_dev, AllSubs_MMN.s6.all_dev, AllSubs_MMN.s7.all_dev, AllSubs_MMN.s8.all_dev, ...
    AllSubs_MMN.s9.all_dev, AllSubs_MMN.s10.all_dev, AllSubs_MMN.s11.all_dev, AllSubs_MMN.s12.all_dev, ...
    AllSubs_MMN.s13.all_dev, AllSubs_MMN.s14.all_dev, AllSubs_MMN.s15.all_dev);

%% Plots

% Determine the min-max range
MinMax = [-1.1 1.4];
    
figure
set(gcf,'Color','w')
subplot_tight(2,1,1);
for c = 1:length(GrandMMN.label)
    plot(GrandMMN.time, GrandMMN.avg(c,:), 'k', 'LineWidth', 1.0); hold all;
end
plot(GrandMMN.time, GrandMMN.avg([19 29],:), 'y', 'LineWidth', 2.0); hold all; % TP10 and TP9
plot(GrandMMN.time, GrandMMN.avg(52,:), 'b', 'LineWidth', 2.0); hold all; % Fz

axis([-0.1,0.4, MinMax])
line([-0.0 -0.0],   get(gca, 'ylim'), 'color', [1.0 0.0 0.0], 'LineWidth', 2.5); hold all;
run('figure_prop.m')
print([PATH, 'aMMN_diffWave'],'-dpng','-r300');

%%

cfg = [];
cfg.keepindividual = 'yes';

GrandStandard_wt = ft_timelockgrandaverage(cfg, AllSubs_MMN.s1.standards, AllSubs_MMN.s2.standards, AllSubs_MMN.s3.standards, ...
    AllSubs_MMN.s4.standards, AllSubs_MMN.s5.standards, AllSubs_MMN.s6.standards, AllSubs_MMN.s7.standards, AllSubs_MMN.s8.standards, ...
    AllSubs_MMN.s9.standards, AllSubs_MMN.s10.standards, AllSubs_MMN.s11.standards, AllSubs_MMN.s12.standards, ...
    AllSubs_MMN.s13.standards, AllSubs_MMN.s14.standards, AllSubs_MMN.s15.standards);

GrandDeviant_wt = ft_timelockgrandaverage(cfg, AllSubs_MMN.s1.all_dev, AllSubs_MMN.s2.all_dev, AllSubs_MMN.s3.all_dev, ...
    AllSubs_MMN.s4.all_dev, AllSubs_MMN.s5.all_dev, AllSubs_MMN.s6.all_dev, AllSubs_MMN.s7.all_dev, AllSubs_MMN.s8.all_dev, ...
    AllSubs_MMN.s9.all_dev, AllSubs_MMN.s10.all_dev, AllSubs_MMN.s11.all_dev, AllSubs_MMN.s12.all_dev, ...
    AllSubs_MMN.s13.all_dev, AllSubs_MMN.s14.all_dev, AllSubs_MMN.s15.all_dev);

% Prepare the neighborhood structure
cfg = [];
cfg.method = 'distance';
cfg.neighbourdist = 0.12;
cfg.layout = 'M1_XYZ3_new.sfp';
neighbours = ft_prepare_neighbours(cfg, AllSubs_MMN.s1.standards);

% stats for overall MMN
cfg = [];
cfg.channel          = 'all';
cfg.latency          = [-0.0 0.4]; % [-0.1 0.0];
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

overall_aMMN_stat = ft_timelockstatistics(cfg, GrandDeviant_wt, GrandStandard_wt);
save([PATH, 'overall_aMMN_stat'], 'overall_aMMN_stat', 'cfg'); % 'aMMN_stat_baseline'
%     contourf(overall_aMMN_stat.mask)
%%
load ([PATH, 'overall_aMMN_stat'], 'overall_aMMN_stat');

[overall_aMMN_stat.posclusters.prob]
[overall_aMMN_stat.negclusters.prob]

overall_aMMN_stat.posclusterslabelmat = overall_aMMN_stat.mask .* overall_aMMN_stat.posclusterslabelmat;
overall_aMMN_stat.negclusterslabelmat = overall_aMMN_stat.mask .* overall_aMMN_stat.negclusterslabelmat;

overall_aMMN_stat.posclusterslabelmat(overall_aMMN_stat.posclusterslabelmat == 1) = 0;

overall_aMMN_stat.mask = overall_aMMN_stat.negclusterslabelmat + overall_aMMN_stat.posclusterslabelmat;
contourf(overall_aMMN_stat.mask)

%% Single channel effects

overall_aMMN_1st_neg_chan = {'Fz'};
overall_aMMN_1st_pos_chan = {'Pz'};

% Fz
fuu = subplot(2,1,1); hold all;
cfg = [];
cfg.channel    = 'Fz';
cfg.linestyle  = {'-','-'};
cfg.graphcolor = [0.5 0.5 0.5; 0 0 0];
cfg.linewidth  = 2.5;
cfg.figure     = fuu;
ft_singleplotER(cfg, GrandStandard, GrandDeviant); hold all;

for i = find(overall_aMMN_stat.posclusterslabelmat(find(strcmp(overall_aMMN_stat.label, overall_aMMN_1st_neg_chan{1,1})),:) > 0)
    plot(overall_aMMN_stat.time(i), -1.3, 'r.', 'MarkerSize',20); hold all;
end
for i = find(overall_aMMN_stat.negclusterslabelmat(find(strcmp(overall_aMMN_stat.label, overall_aMMN_1st_neg_chan{1,1})),:) > 0)
    plot(overall_aMMN_stat.time(i), -1.3, 'b.', 'MarkerSize',20); hold all;
end
line([-0.0 -0.0],   get(gca, 'ylim'), 'color', [0.0 0.0 0.0], 'LineWidth', 1.0); hold all;
axis([-0.1,0.4, -1.3,1.3])
title('')
run('figure_prop.m')

print([PATH, 'aMMN_PosNeg'],'-dpng','-r300');

%% Topoplots post-stimulus

figure
set(gcf,'Color','w')

count = 1;
j = [0.0:0.1:0.4];
for k = 1:4
    
    fuu = subplot_tight(3,4,count); count = count+1; hold all;
    
    cfg = [];
    cfg.xlim = [j(k) j(k+1)];
    cfg.zlim = [-1 1];
    cfg.comment = 'no';
    cfg.colorbar = 'no';
    cfg.style = 'straight';
    cfg.layout = 'M1_XYZ3_new.sfp';
    cfg.marker = 'on';
    cfg.markersymbol = '.';
    %cfg.markersize = 1;
    cfg.figure = fuu;
    
    ft_topoplotER(cfg, GrandMMN); hold all;
end

j = [0.0:0.1:0.4];
i = [1:100:401];
for k = 1:4
    
    fuu = subplot_tight(3,4,count); count = count+1; % 4 Rows so it would look like the GlobLoc ones
    
    cfg = [];
    cfg.xlim = [j(k) j(k+1)];
    cfg.comment = 'no';
    cfg.colorbar = 'no';
    cfg.style = 'blank';
    cfg.highlight = {'on', 'on'};
    cfg.highlightchannel = {overall_aMMN_stat.label(sum(overall_aMMN_stat.negclusterslabelmat(:, i(k):i(k+1)) > 0, 2) > 33), ...
                            overall_aMMN_stat.label(sum(overall_aMMN_stat.posclusterslabelmat(:, i(k):i(k+1)) > 0, 2) > 33)};
    cfg.highlightsymbol = {'.', '.'};
    cfg.highlightsize = {45, 45};
    cfg.highlightcolor = {[0 0 1], [1 0 0]};
    cfg.layout = 'M1_XYZ3_new.sfp';
    cfg.marker = 'off';
    cfg.figure = fuu;
   
    ft_topoplotER(cfg, GrandMMN); hold all;
    
end

print([PATH, 'aMMN_diffTopoStat'],'-dpng','-r300');

%% STATS 4 PAPER

chan = 52; %find( sum(overall_aMMN_stat.negclusterslabelmat,2));
linz = mean(GrandMMN.avg(chan,:),1);

GrandMMN.time(find(linz == min(linz),1,'first'))
linz(find(linz == min(linz),1,'first'))

%% Location of the electrodes with max. statistical results

cfg = [];
cfg.xlim = [j(k) j(k+1)];
cfg.comment = 'no';
cfg.colorbar = 'no';
cfg.style = 'blank';
cfg.highlight = {'on', 'on'};
cfg.highlightchannel = {overall_aMMN_1st_pos_chan, overall_aMMN_1st_neg_chan};
cfg.highlightsymbol = { '.', '.'};
cfg.highlightsize = {50, 50};
cfg.highlightcolor = {[0 0 0], [0 0 0]};
cfg.layout = 'M1_XYZ3_new.sfp';
cfg.marker = 'off';

ft_topoplotER(cfg, GrandMMN);
    
print([PATH, 'aMMN_PosNegTopo'],'-dpng','-r300');
%%
% FOR PEAKS, MEAN AMPLITUDE, ETC.
%%

cfg = [];
cfg.keepindividual = 'yes';
GrandMMN_wt = ft_timelockgrandaverage(cfg, AllSubs_MMN.s1.allDiff, AllSubs_MMN.s2.allDiff, AllSubs_MMN.s3.allDiff, ...
    AllSubs_MMN.s4.allDiff, AllSubs_MMN.s5.allDiff, AllSubs_MMN.s6.allDiff, AllSubs_MMN.s7.allDiff, AllSubs_MMN.s8.allDiff, ...
    AllSubs_MMN.s9.allDiff, AllSubs_MMN.s10.allDiff, AllSubs_MMN.s11.allDiff, AllSubs_MMN.s12.allDiff, ...
    AllSubs_MMN.s13.allDiff, AllSubs_MMN.s14.allDiff, AllSubs_MMN.s15.allDiff);

fuu = zeros(15,501);
fuuMin = zeros(15,2);
CLUSTER = zeros(62,1);
subplot(3,1,1:2); hold all;
for s = [1:15 1]
    load([PATH, 'pp_', num2str(s), '\\aMMN_stat']);
    
    GrandNegCluster = 52;
%     GrandNegCluster = find( sum(aMMN_stat.negclusterslabelmat == 1,2) );
    CLUSTER(GrandNegCluster) = CLUSTER(GrandNegCluster) + 1; 

    fuu(s,:) = squeeze(mean(GrandMMN_wt.individual(s,GrandNegCluster,:),2))';
    fuuMin(s,1) = GrandMMN_wt.time(find(fuu(s,:) == min(fuu(s,:)), 1, 'first'));
    fuuMin(s,2) = fuu(s,find(fuu(s,:) == min(fuu(s,:)), 1, 'first'));
    
    plot(GrandMMN_wt.time, fuu(s,:), 'color', [0.5 0.5 0.5], 'LineWidth', 2.0); hold all;
    
    if s == 1
        plot(GrandMMN_wt.time, fuu(s,:), 'color', [0.3 0.3 0.3], 'LineWidth', 2.0); hold all;
        aMMN_stat.negclusterslabelmat(aMMN_stat.negclusterslabelmat ~= 1) = 0;
        bar = find(sum(aMMN_stat.negclusterslabelmat(52,1:210),1) > 0);
        plot(GrandMMN_wt.time(bar+100), fuu(s,bar+100), 'color', [1.0 0.0 0.0], 'LineWidth', 2.0); hold all;
    else
        aMMN_stat.negclusterslabelmat = aMMN_stat.mask .* aMMN_stat.negclusterslabelmat;
        bar = find(sum(aMMN_stat.negclusterslabelmat(52,1:210),1) > 0);
        plot(GrandMMN_wt.time(bar+100), fuu(s,bar+100), 'color', [0.0 0.0 1.0], 'LineWidth', 2.0); hold all;
    end
end
axis([-0.1,0.4, -2,2])
line([-0.0 -0.0],   get(gca, 'ylim'), 'color', [1.0 0.0 0.0], 'LineWidth', 0.5); hold off;
run('figure_prop.m')
print([PATH, 'aMMN_15subClusters'],'-dpng','-r300');

%%

mean(fuuMin(:,1))
std(fuuMin(:,1))
[min(fuuMin(:,1)) max(fuuMin(:,1))]

mean(fuuMin(:,2))
std(fuuMin(:,2))
[min(fuuMin(:,2)) max(fuuMin(:,2))]

subplot(2,1,1);
scatter(fuuMin(:,1), fuuMin(:,2)); hold all; 
axis([0.07,0.2, -2,0]); hold all; 
for s = 1:15
    line([fuuMin(s,1) fuuMin(s,1)], [fuuMin(s,2),0], 'color', [0.5 0.5 0.5], 'LineWidth', 1.0); hold all;
end
%run('figure_prop.m')
print([PATH, 'aMMN_indVals'],'-dpng','-r300');

%%

% fuuStat = zeros(15,62,401);
% fuuFlat = zeros(15,401);
% for s = 15
%     load([PATH, 'CF_S', subNr{s}, '_CTR/', 'aMMN_stat']);
%     contourf(aMMN_stat.mask,2)
%     
%     mask = zeros(size(aMMN_stat.mask));
%     %mask( aMMN_stat.posclusterslabelmat == 1 ) = 1;
%     %mask( aMMN_stat.posclusterslabelmat == 2 ) = 1;
%     mask( aMMN_stat.negclusterslabelmat == 1 ) = -1;
%     contourf(mask,2)
% 
%     fuuStat(s,:,:) = mask;
%     fuuFlat(s,:) = sum(mask,1);
%     
% end
% save([PATH, 'indiv_aMMN_contourf'], 'fuuStat');
load([PATH, 'indiv_aMMN_contourf'], 'fuuStat');

mask = zeros(size(overall_aMMN_stat.mask));
mask( overall_aMMN_stat.posclusterslabelmat == 2 ) = 1;
mask( overall_aMMN_stat.negclusterslabelmat == 1 ) = -1;
fuuStat(16,:,:) = mask;
% fuuFlat(find(fuuFlat > 0)) = 1;
% fuuFlat(find(fuuFlat < 0)) = -1;
fuuFlat = squeeze(sum(fuuStat,2));
% contourf(fuuFlat,2);
subplot(3,1,1:2); hold all;
for s = 1:15
    
    Cneg = [find(fuuFlat(s,:) <= -3, 1, 'first') ... % onset & offset
            find(fuuFlat(s,:) <= -3, 1, 'last')];
        
    line([Cneg(1) Cneg(2)], [s,s], 'color', [0 0 1], 'LineWidth', 10.0); hold all;
    

    Cpos = [find(fuuFlat(s,:) >= 3, 1, 'first') ... % onset & offset
            find(fuuFlat(s,:) >= 3, 1, 'last')];
    if ~isempty(Cpos)
        line([Cpos(1) Cpos(2)], [s,s], 'color', [1 0 0], 'LineWidth', 10.0); hold all;
    end

end
axis([-100,400, -1,16])
line([-0.0 -0.0],   get(gca, 'ylim'), 'color', [0.0 0.0 0.0], 'LineWidth', 1); hold off;
% run('figure_prop.m')
print([PATH, 'aMMN_indStats'],'-dpng','-r300');
