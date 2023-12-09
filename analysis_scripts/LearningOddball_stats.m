%% Clear the workspace
% ft_defaults
clear all
clc

%% 

PATH    = '.........................\\single_trial_data_LearningOddball\\';

for s = 1:15
    
    load([PATH, 'pp_', num2str(s), '\\LearningOddball']) 

    % Define outlier trials and remove them
    [idx,d] = nt_find_outlier_trials( nt_trial2mat( good_data.trial ), 2 , 0);  % 2SD
    cfg = [];
    cfg.trials = idx;
    good_data  = ft_redefinetrial(cfg, good_data);
    CNDs = CNDs(idx,:);

    
    cfg = [];
    cfg.demean         = 'yes';                         % baseline correction
    cfg.baselinewindow = [-0.1 0.0]; 
    cfg.trials     = find(CNDs(:,1) == 1);
    standards      = ft_preprocessing(cfg, good_data);
    cfg.trials     = find(CNDs(:,1) == 2); 
    deviants       = ft_preprocessing(cfg, good_data);

%     cfg = [];
%     cfg.keeptrials = 'no';
%     all_stand_avg  = ft_timelockanalysis(cfg, standards);
%     all_dev_avg    = ft_timelockanalysis(cfg, deviants);
% 
%     cfg = [];
%     cfg.layout = 'M1_XYZ3_new.sfp';
%     ft_multiplotER(cfg, all_stand_avg, all_dev_avg);

    clear good_data

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfg = [];
    cfg.method = 'distance';
    cfg.neighbourdist = 0.12;
    cfg.layout = 'M1_XYZ3_new.sfp';
    neighbours = ft_prepare_neighbours(cfg, standards);

    cfg = [];
    cfg.keeptrials = 'yes';
    standards_wt = ft_timelockanalysis(cfg, standards);
    deviants_wt  = ft_timelockanalysis(cfg, deviants);

    if size(standards_wt.trial,1) > size(deviants_wt.trial,1)
        ntrl = randsample(size(standards_wt.trial,1),size(deviants_wt.trial,1));
        standards_wt.trial = standards_wt.trial(ntrl,:,:);
    elseif size(standards_wt.trial,1) < size(deviants_wt.trial,1)
        ntrl = randsample(size(deviants_wt.trial,1),size(standards_wt.trial,1));
        deviants_wt.trial = deviants_wt.trial(ntrl,:,:);
    end

    
    cfg = [];
    cfg.channel          = 'all';
    cfg.latency          = [0.0 0.7];
    cfg.method           = 'montecarlo';
    cfg.parameter        = 'trial';
    cfg.statistic        = 'indepsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan        = 2;
    cfg.neighbours       = neighbours;
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 10000;


    design = [repmat(1,1,size(standards_wt.trial,1)), repmat(2,1,size(deviants_wt.trial,1))];

    cfg.design = design;
    cfg.ivar   = 1;

    aP3b_stat = ft_timelockstatistics(cfg, deviants_wt, standards_wt);

    save([PATH, 'pp_', num2str(s), '\\aP3b_stat'], 'aMMN_stat', 'cfg'); 
    %load([PATH, SUBJECT, 'aP3b_stat']);
    
%     figure
%     contourf(aP3b_stat.mask)
%     figure
% %     contourf(aP3b_stat.negclusterslabelmat==1)
%     contourf(aP3b_stat.posclusterslabelmat==1)
end
