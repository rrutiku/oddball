%% Clear the workspace
ft_defaults                           
clear all
clc
%% Define where the data is
PATH    = 'D:\\_W\\itcf-lab\\oddball_comparison\\data\\';	% Project folder
SUBJECT = '.........\\';                                    % Subject folder
DATSETS = {};                                               % Datasets for this particular paradigm
for d = 2:3 % Check these numbers!!!
    DATSETS{end+1} = ['......................', num2str(d)];
end
% load([PATH, SUBJECT, 'GlobLoc_trl'])
%% Define trials

blocks = {1:4, 5:8};                                                        % Which blocks are in each dataset?
run('GlobLoc_prep2500_xtra')                                                % Execute the extra script for trigger management

% Check the diffs for irrelevant triggers at the beginning
% Sanity checks:
% fuu = diff(events)/2500;
% mean(fuu(fuu < 0.7)); [min(fuu(fuu < 0.7)) max(fuu(fuu < 0.7))]; hist(fuu(fuu < 0.7))
% [min(diffs{1}), min(diffs{2})]
% bar = diff(tz(1:find(tz,1,'last'))); hist(bar(bar < 0.3))
%%
CFGs    = {};
CFGs{1} = cfgs;
% If not yet present, start a '_trl' file
if exist([PATH, SUBJECT, 'GlobLoc_trl.mat'], 'file') == 0
    save([PATH, SUBJECT, 'GlobLoc_trl'], 'CFGs')
else
	save([PATH, SUBJECT, 'GlobLoc_trl'], 'CFGs', '-append')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in the data

dat = {};
for ds = 1:length(blocks)
    
    CFGs{end+1} = [];
    CFGs{end}.dataset    = [PATH, SUBJECT, DATSETS{ds}, '.vhdr'];
    CFGs{end}.continuous = 'yes';   
    dat{ds} = ft_preprocessing(CFGs{end});  % load the raw data
    
%     dat{ds}.trial{1,1}([19,29],:) = dat{ds}.trial{1,1}([29,19],:); 

    % add a hp-filter here if necessary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Browse through the raw data
    cfg = [];
    cfg.channel = {'all' '-VEOG' '-HEOG'};
    cfg.continuous = 'yes';
    cfg.layout = 'M1_XYZ3_new.sfp';
    cfg.blocksize = 10;
    cfg.preproc.demean = 'yes';
    cfg = ft_databrowser(cfg, dat{ds});

end

%% Are TP9 and TP10 switched?

R = corrcoef(dat{1}.trial{1,1}(19,:),dat{1}.trial{1,1}(29,:))

R = corrcoef(dat{1}.trial{1,1}(19,:),dat{1}.trial{1,1}(20,:))
R = corrcoef(dat{1}.trial{1,1}(19,:),dat{1}.trial{1,1}(28,:))


R = corrcoef(dat{1}.trial{1,1}(29,:),dat{1}.trial{1,1}(28,:))
R = corrcoef(dat{1}.trial{1,1}(29,:),dat{1}.trial{1,1}(20,:))
%%
for ds = 1:length(blocks)
    
    CFGs{end+1}   = [];
    CFGs{end}.trl = cfgs{ds}.trl;
    dat{ds} = ft_redefinetrial(CFGs{end}, dat{ds});
    
    if ds == 1
        raw_data  = dat{ds};
    else
        cfg = [];
        raw_data  = ft_appenddata(cfg, raw_data, dat{ds});
    end
    
end; clear dat
%FORDS = cell2mat({cat(1, CHORDS{:})});
%FORDS   = FORDS(FORDS(:,2)==1,1); 
%% downsample the data and add bipolar EOG + electrode info

% downsample the data
CFGs{end+1} = [];
CFGs{end}.resamplefs  = 1000;
CFGs{end}.detrend     = 'no';
CFGs{end}.demean      = 'yes'; 
raw_data      = ft_resampledata(CFGs{end}, raw_data);

% Remove edges because there are artifacts and one should not filter over them
CFGs{end+1} = [];
CFGs{end}.toilim = [-0.9 0.8];
raw_data  = ft_redefinetrial(CFGs{end}, raw_data);

% filter the data below muscle activity
CFGs{end+1} = [];
CFGs{end}.lpfilter  = 'yes';
CFGs{end}.lpfreq    = 30;
CFGs{end}.hpfilter  = 'no';
CFGs{end}.detrend   = 'no';
raw_data = ft_preprocessing(CFGs{end}, raw_data);

% Remove edges again due to possible edge artifacts
CFGs{end+1} = [];
CFGs{end}.toilim = [-0.5 0.7];
raw_data  = ft_redefinetrial(CFGs{end}, raw_data);

% Change sampleinfo to avoid repeating sample numbers
raw_data.sampleinfo = [1000:2000:2000*length(raw_data.trial)]';
raw_data.sampleinfo(:,2) = raw_data.sampleinfo(:,1)+length(raw_data.time{1})-1;
fake_trl_raw = raw_data.sampleinfo;
raw_data.cfg.event = [];

% Add electrode locations and bipolar EOG
raw_data = bipolar_EOG(raw_data);
save([PATH, SUBJECT, 'GlobLoc_trl'], 'fake_trl_raw', '-append')

%% Look at the data
global badchan allbadchan
badchan = {}; allbadchan = {};

% artifact rejection 1
cfg = [];
cfg.channel    = {'all' '-VEOG' '-HEOG'};
cfg.continuous = 'no';
cfg.ylim = [-250 250];
cfg.layout = 'M1_XYZ3_new.sfp';              % specify the layout file that should be used for optional topoplotting
cfg.viewmode = 'butterfly';  % 'vertical'; % 
cfg = ft_databrowser(cfg,raw_data);

artifacts_1  = cfg.artfctdef.visual.artifact;
badchan_1    = badchan;
allbadchan_1 = allbadchan{end};
save([PATH, SUBJECT, 'GlobLoc_trl'], 'artifacts_1', 'badchan_1', 'allbadchan_1', '-append')
%% Repair the bad channels and reject artifacts

% Interpolate faulty channels
raw_data = repair_channels(raw_data, badchan_1, allbadchan_1); % % Note: deletes the EOG channels

cfg = [];
cfg.artfctdef.eog.artifact = artifacts_1;
cfg = ft_rejectartifact(cfg,raw_data);
pre_finals = cfg.cfg.trl;

CFGs{end+1} = [];
CFGs{end}.trl = pre_finals;
raw_data = ft_redefinetrial(CFGs{end},raw_data);

save([PATH, SUBJECT, 'GlobLoc_trl'], 'CFGs', 'pre_finals', '-append')

%% Perform ICA 1

tmpdata   = cat(3, raw_data.trial{1:10});
tmpdata   = reshape(tmpdata, [size(raw_data.trial{1},1) size(raw_data.trial{1},2)*10]);
tempCpca  = svd(tmpdata);
th = 0.01;
figure,semilogy(tempCpca,'.-')
a = ylim;
hold on
Cpca = length(raw_data.label) - length(allbadchan_1) +2;
line([Cpca Cpca],[a(1) a(2)],'Color','k','LineStyle','--');
hold off

%%
cfg = [];
cfg.channel = {'all' '-VEOG' '-HEOG'};
cfg.demean = 'yes';              
cfg.method = 'runica';          % FieldTrip supports multiple ways to perform ICA, 'runica' is one of them.
% cfg.runica.extended = 1;
cfg.runica.pca = Cpca;
cfg.updatesens = 'no';

comp = ft_componentanalysis(cfg, raw_data);
% save([PATH, SUBJECT, 'GlobLoc_trl'],'comp'); % load([PATH, SUBJECT, 'GlobLoc_trl'],'comp');

%% 

EEG = fieldtrip2eeglab(raw_data);

EEG.compvars    = repmat(100,1,Cpca);
EEG.icasphere   = eye(62);
EEG.icaweights  = comp.unmixing;
EEG.icawinv     = comp.topo;
EEG.icaact      = cat(3,comp.trial{:});
EEG.comp2remove = [];        

%%
% second trials rejection
TMPREJ=[];
eegplot(EEG.icaact,'winlength',5,'command','pippo');
%%
if ~isempty(TMPREJ)
    [trialrej elecrej]=eegplot2trial(TMPREJ,size(EEG.data,2),size(EEG.data,3));
else
    trialrej=[];
end
tr2reject=find(trialrej==1);
%%
pre_finals_2 = pre_finals;
pre_finals_2(tr2reject,:) = [];

cfg = [];
cfg.trl = pre_finals_2;
raw_data = ft_redefinetrial(cfg, raw_data);
% If this change was made rerun ICA
save([PATH, SUBJECT, 'GlobLoc_trl'], 'tr2reject', 'pre_finals_2', '-append')

%% Perform ICA 2

cfg = [];
cfg.channel = {'all' '-VEOG' '-HEOG'};
cfg.demean = 'yes';              
cfg.method = 'runica';          % FieldTrip supports multiple ways to perform ICA, 'runica' is one of them.
% cfg.runica.extended = 1;
cfg.runica.pca = Cpca;
cfg.updatesens = 'no';

comp = ft_componentanalysis(cfg, raw_data);
save([PATH, SUBJECT, 'GlobLoc_trl'],'comp'); % load([PATH, SUBJECT, DATSET, '_comp'],'comp');

%% 

EEG = fieldtrip2eeglab(raw_data);

EEG.compvars    = repmat(100,1,Cpca);
EEG.icasphere   = eye(62);
EEG.icaweights  = comp.unmixing;
EEG.icawinv     = comp.topo;
EEG.icaact      = cat(3,comp.trial{:});
EEG.comp2remove = [];        

% cd Preprocessing_BRAINAMP/
ICA_selectandremove(EEG)
% save([PATH, SUBJECT, DATSET, '_EEGstruct'], 'EEG', 'data_GUI');

%%
comp_out = data_GUI.comp2remove;
save([PATH, SUBJECT, 'GlobLoc_trl'], 'comp_out', '-append')

%% Be sure that comps_out contains only those components that should be cast out

% remove the bad components and backproject the data
CFGs{end+1}   = [];
CFGs{end}.component = comp_out; % to be removed component(s)
CFGs{end}.updatesens = 'no';
data = ft_rejectcomponent(CFGs{end}, comp);

%% Preprocess the data with the options for ERP analysis

CFGs{end+1}   = [];
CFGs{end}.reref       = 'yes';
CFGs{end}.refchannel  = 'all'; % 

data = ft_preprocessing(CFGs{end}, data);

%% Look at filtered and repaired data
badchan = {}; allbadchan = {};

cfg = [];
cfg.channel    =  {'all','-VEOG','-HEOG'};
cfg.continuous = 'no';
cfg.ylim       = [-100 100];
cfg.viewmode = 'butterfly';  % 'vertical'; % 
cfg = ft_databrowser(cfg, data);

artifacts_2  = cfg.artfctdef.visual.artifact;
badchan_2    = badchan;
allbadchan_2 = {}; allbadchan_2 = allbadchan{end}; % 
%% Repair the bad channels and reject artifacts

save([PATH, SUBJECT, 'GlobLoc_trl'], 'artifacts_2', 'badchan_2', 'allbadchan_2', '-append')

% Interpolate faulty channels
data.elec = raw_data.elec;
data = repair_channels(data, badchan_2, allbadchan_2); % % Note: deletes the EOG channels

cfg = [];
cfg.artfctdef.eog.artifact = artifacts_2;
cfg = ft_rejectartifact(cfg,data);
finals = cfg.cfg.trl;

CFGs{end+1} = [];
CFGs{end}.trl   = finals;
good_data = ft_redefinetrial(CFGs{end},data);

save([PATH, SUBJECT, 'GlobLoc_trl'], 'CFGs', 'finals', '-append');
%% Final save

save([PATH, SUBJECT, 'GlobLoc_trl'], 'CHORDS', 'CFGs', 'missTrig', 'fake_trl_raw', ...
     'artifacts_1', 'badchan_1', 'allbadchan_1', 'pre_finals', 'tr2reject', 'pre_finals_2', 'comp_out', ...
     'artifacts_2', 'badchan_2', 'allbadchan_2', 'finals');
 
save([PATH, SUBJECT, 'GlobLoc'], 'good_data');

%%
clear all
clc