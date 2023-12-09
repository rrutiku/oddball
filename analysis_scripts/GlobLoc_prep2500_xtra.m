%% Load the stimulation sequence and recreate the conditions

d = dir([PATH, SUBJECT, 'logs/GlobLoc_stimulation_SEQUENCE*']);
load([PATH, SUBJECT, 'logs/', d.name], 'SEQUENCE')

CHORDS = {};
for ds = 1:length(blocks)
    
    CHORDS{ds} = [];
    
    for i = blocks{ds}
    
        standard = SEQUENCE{1,i}(1,:);
    
        if isequal(standard, [1 1 1 1 1])
           block = 1;
       elseif isequal(standard, [2 2 2 2 2])
           block = 2;
       elseif isequal(standard, [1 1 1 1 2])
           block = 3;
       elseif isequal(standard, [2 2 2 2 1])
           block = 4;     
        end
   
        thc = repmat(block,size(SEQUENCE{1,i},1),5);
        thc(SEQUENCE{1,i}(:,5) ~= standard(5),5) = block*-1;
        thc = reshape(thc.',1,[]);
    
        % for the 5th triggers afterwards
        trg = repmat([0 0 0 0 1],size(SEQUENCE{1,i},1),1);
        trg = reshape(trg.',1,[]);
    
        CHORDS{ds} = [CHORDS{ds}; [thc' trg']];
    end

end

%% Identify triggers and check them

cfg = [];
cfg.trialdef.eventtype  = 'Response';
cfg.trialdef.eventvalue = 'R128';
cfg.continuous          = 'yes';   
cfg.trialdef.prestim    = 1.0; % in seconds
cfg.trialdef.poststim   = 0.9; % in seconds

cfgs = {}; diffs = {};
for ds = 1:length(blocks)

    cfg.dataset = [PATH, SUBJECT, DATSETS{ds}, '.vhdr'];
    cfgs{ds}    = ft_definetrial(cfg);
    
    diffs{ds} = diff(cfgs{ds}.trl(:,1)) / 2500;
    
end

% cfgs{1}.trl = cfgs{1}.trl(1:end-10,:);    % For S05 (with -1 to 0.9 sec)

%% Look at the diffs and decide whether some triggers at the beginning do not belong to the paradigm

% should be around a sequence of aprox. 0.15 0.15 0.15 0.15 >1.0

% 1st is the test trigger?
% cfgs{...}.trl(1,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate the missing triggers based on psychtoolbox output

% This is the psychtoolbox output for each stimulus
d = dir([PATH, SUBJECT, 'logs/GlobLoc_critical_events*']);
load([PATH, SUBJECT, 'logs/', d.name]) % diff(tz)';

count = 1;
missTrig = {};

% % For S02, S04, S06
% CHORDS{2} = CHORDS{2}(2:end,:);
% tz = tz( [ 1:length(CHORDS{1})  length(CHORDS{1})+2:find(tz,1,'last') ] );

% % For S05 (with -1 to 0.9 sec)
% CHORDS{1} = CHORDS{1}(6:end,:);
% CHORDS{2} = CHORDS{2}(84:end,:);
% CHORDS{3} = CHORDS{3}(2:end,:);
% tz = tz( [ 6:2050   2134:2775  2777:find(tz,1,'last') ] );

for ds = 1:length(blocks)
    
    tzz = tz(count:count+size(CHORDS{ds},1)-1); % diff(tz); diff(tzz);
    
    TZ_line = zeros(1,round((tzz(end)-tzz(1))*2500)+100);
    TZ_line( round((tzz-tzz(1))*2500.15) +1 ) = 2;                  % .15 correction of sampling rate (2500 Hz)
    TZ = find(TZ_line)' + cfgs{ds}.trl(1,1)-1;
    
    events   = cfgs{ds}.trl(:,1);
    
    missTrig{ds} = [];
    for i = 1:size(TZ,1)
        if i > numel(events)
            missTrig{ds} = [missTrig{ds}; i];
            vai2ins = TZ(i,1) + events(i-1,1) - TZ(i-1,1); % estimate the sample number
            events = [events(1:i-1);   vai2ins];
            cfgs{ds}.trl = [cfgs{ds}.trl(1:i-1,:);   [vai2ins vai2ins+(cfgs{ds}.trl(1,2)-cfgs{ds}.trl(1,1)) ...
                cfgs{ds}.trl(1,3) cfgs{ds}.trl(1,4)]];
        end
        if TZ(i,1)>events(i,1)-50 && TZ(i,1)<events(i,1)+50 % if the sample number is about the same % 50 for 2500 Hz
            continue % Do nothing
        else
            missTrig{ds} = [missTrig{ds}; i];
            vai2ins = TZ(i,1) + events(i-1,1) - TZ(i-1,1); % estimate the sample number
            events = [events(1:i-1);   vai2ins;    events(i:end)];
            cfgs{ds}.trl = [cfgs{ds}.trl(1:i-1,:);   [vai2ins vai2ins+(cfgs{ds}.trl(1,2)-cfgs{ds}.trl(1,1)) ...
                cfgs{ds}.trl(1,3) cfgs{ds}.trl(1,4)];    cfgs{ds}.trl(i:end,:)];
        end
    end

    count = count + size(CHORDS{ds},1);
end

%% Define only the 5th tone triggers as trials

for ds = 1:length(blocks)
    cfgs{ds}.trl = cfgs{ds}.trl(CHORDS{ds}(:,2)==1,:); 
end


