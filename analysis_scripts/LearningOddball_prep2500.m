%% Clear the workspace
ft_defaults                            
clear all
clc
%% Define where the data is
PATH    = 'D:\\_W\\itcf-lab\\oddball_comparison\\data\\';	% Project folder
SUBJECT = '..........\\';                                   % Subject folder
DATSET  = '.........................';                      % Dataset for this particular paradigm
% load([PATH, SUBJECT, DATSET, '_trl'])
%% Define the trials 
CFGs    = {};
CFGs{1} = [];
CFGs{end}.dataset             = [PATH, SUBJECT, DATSET, '.vhdr'];
CFGs{end}.trialdef.eventtype  = 'Response';
CFGs{end}.trialdef.eventvalue = 'R128';
CFGs{end}.continuous          = 'yes';   
CFGs{end}.trialdef.prestim    = 1.0;          % in seconds; includes space for edge artifacts (x2 due to resampling ect)
CFGs{end}.trialdef.poststim   = 1.0;          % in seconds; includes space for edge artifacts

CFGs{end} = ft_definetrial(CFGs{end});
%% If not yet present, start a '_trl' file
if exist([PATH, SUBJECT, DATSET, '_trl.mat'], 'file') == 0
    save([PATH, SUBJECT, DATSET, '_trl'], 'CFGs')
else
	save([PATH, SUBJECT, DATSET, '_trl'], 'CFGs', '-append')
end
%% Make sure that all the above used triggers were actually stimulus triggers
% (and not test triggers at the beginning of the experiment, for example)

difts = diff(CFGs{end}.trl(:,1)) / 2500; % should be around 1.05
% min(difts); hist(difts(difts < 1.7))
% Remove the irrelevant triggers from the definition
% cfg.trl([1],:) = [];
%% Interpolate the missing triggers based on psychtoolbox output

% This is the psychtoolbox output for each stimulus
d = dir([PATH, SUBJECT, 'logs/aP3b_critical_events*']);
load([PATH, SUBJECT, 'logs/', d.name])

tz = [t0z tz];

TZ_line = zeros(1,round((tz(end)-tz(1))*2500)+100);
TZ_line( round((tz-tz(1))*2500.15) +1 ) = 2;            % .15 correction of sampling rate (2500 Hz)
TZ = find(TZ_line)' + CFGs{end}.trl(1,1)-1;

events   = CFGs{end}.trl(:,1);

missTrig = [];
for i = 1:size(TZ,1)
    if TZ(i,1)>events(i,1)-50 && TZ(i,1)<events(i,1)+50 % if the sample number is about the same
        continue 
    else
        missTrig = [missTrig; i];
        vai2ins = TZ(i,1) + events(i-1,1) - TZ(i-1,1);  % estimate the sample number
        events = [events(1:i-1);   vai2ins;    events(i:end)];
        CFGs{end}.trl = [CFGs{end}.trl(1:i-1,:);   [vai2ins vai2ins+(CFGs{end}.trl(1,2)-CFGs{end}.trl(1,1)) ...
                    CFGs{end}.trl(1,3) CFGs{end}.trl(1,4)];    CFGs{end}.trl(i:end,:)];
    end
end

trl_raw = CFGs{end}.trl;
% TZ(:,2) = events; TZ(:,3) = TZ(:,1) - TZ(:,2);
% fuu = diff(events)/2500; hist(fuu)
save([PATH, SUBJECT, DATSET, '_trl'], 'trl_raw', 'missTrig', '-append')
%% Load and downsample the data; add bipolar EOG + electrode info 

CFGs{end+1} = [];
CFGs{end}.dataset    = [PATH, SUBJECT, DATSET, '.vhdr'];
CFGs{end}.continuous = 'yes';   
% add a hp-filter here if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFGs{end}.hpfilter  = 'yes';
% CFGs{end}.hpfreq    = 0.5;
% CFGs{end}.hpfiltord = 4;
raw_data = ft_preprocessing(CFGs{end});

% Browse through the raw data
cfg = [];
cfg.channel = {'all' '-VEOG' '-HEOG'};
cfg.continuous = 'yes';
cfg.layout = 'M1_XYZ3_new.sfp';
cfg.blocksize = 10;
cfg.preproc.demean = 'yes';
cfg = ft_databrowser(cfg, raw_data);

%% Are TP9 and TP10 switched?

R = corrcoef(raw_data.trial{1,1}(19,:),raw_data.trial{1,1}(29,:))

R = corrcoef(raw_data.trial{1,1}(19,:),raw_data.trial{1,1}(20,:))
R = corrcoef(raw_data.trial{1,1}(19,:),raw_data.trial{1,1}(28,:))


R = corrcoef(raw_data.trial{1,1}(29,:),raw_data.trial{1,1}(28,:))
R = corrcoef(raw_data.trial{1,1}(29,:),raw_data.trial{1,1}(20,:))

% raw_data.trial{1,1}([19,29],:) = raw_data.trial{1,1}([29,19],:); 
%%

% Cut the data into trials
CFGs{end+1}   = [];
CFGs{end}.trl = trl_raw;
raw_data = ft_redefinetrial(CFGs{end}, raw_data);

% downsample the data
CFGs{end+1} = [];
CFGs{end}.resamplefs  = 1000;
CFGs{end}.detrend     = 'no';
CFGs{end}.demean      = 'yes'; 
raw_data      = ft_resampledata(CFGs{end}, raw_data);

% Remove edges because there are artifacts and one should not filter over them
CFGs{end+1} = [];
CFGs{end}.toilim = [-0.9 0.9];
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
save([PATH, SUBJECT, DATSET, '_trl'], 'fake_trl_raw', '-append')

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
save([PATH, SUBJECT, DATSET, '_trl'], 'artifacts_1', 'badchan_1', 'allbadchan_1', '-append')
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

save([PATH, SUBJECT, DATSET, '_trl'], 'CFGs', 'pre_finals', '-append')

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
% save([PATH, SUBJECT, DATSET, '_comp'],'comp'); % load([PATH, SUBJECT, DATSET, '_comp'],'comp');

%% 

EEG = fieldtrip2eeglab(raw_data);

EEG.compvars    = repmat(100,1,Cpca);
EEG.icasphere   = eye(62);
EEG.icaweights  = comp.unmixing;
EEG.icawinv     = comp.topo;
EEG.icaact      = cat(3,comp.trial{:});
EEG.comp2remove = [];        

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
save([PATH, SUBJECT, DATSET, '_trl'], 'tr2reject', 'pre_finals_2', '-append')

%% Perform ICA 2

cfg = [];
cfg.channel = {'all' '-VEOG' '-HEOG'};
cfg.demean = 'yes';              
cfg.method = 'runica';          % FieldTrip supports multiple ways to perform ICA, 'runica' is one of them.
% cfg.runica.extended = 1;
cfg.runica.pca = Cpca;
cfg.updatesens = 'no';

comp = ft_componentanalysis(cfg, raw_data);
save([PATH, SUBJECT, DATSET, '_comp'],'comp'); % load([PATH, SUBJECT, DATSET, '_compS'],'comp');

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
%% Be sure that comps_out contains only those components that should be cast out

comp_out = data_GUI.comp2remove;
save([PATH, SUBJECT, DATSET, '_trl'], 'comp_out', '-append') 

%%
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

save([PATH, SUBJECT, DATSET, '_trl'], 'artifacts_2', 'badchan_2', 'allbadchan_2', '-append')

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

save([PATH, SUBJECT, DATSET, '_trl'], 'CFGs', 'finals', '-append');
%% Final save

save([PATH, SUBJECT, DATSET, '_trl'], 'CFGs', 'missTrig', 'trl_raw',  'fake_trl_raw', ...
     'artifacts_1', 'badchan_1', 'allbadchan_1', 'pre_finals', 'tr2reject', 'pre_finals_2', 'comp_out', ...
     'artifacts_2', 'badchan_2', 'allbadchan_2', 'finals');
 
save([PATH, SUBJECT, DATSET], 'good_data');

%%
clear all
clc
