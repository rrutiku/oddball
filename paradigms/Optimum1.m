
%{
    Created July 2016
	@author: Renate Rutiku
	contact: renate.rutiku@gmail.com

    ---------------------------------------------------------------------
    ------ Optimum-1 paradigm for the multi-feature MMN assessment ------
    ---------------------------------------------------------------------

    This scrips runs the auditory Optimum-1 MMN paradigm (N‰‰t‰nen et al., 
    2004; Pakarinen et al., 2007) with minor adjustments.
    Most parameter values are replicated from Pakarinen et al. (2007), 
    correspond to the highest of the 6 tested deviation levels.
    
    In order for the script to work:
    1. Psychtoolbox (http://psychtoolbox.org/) has to be installed.
    2. The 'GenerateTone.m' and 'AttenuateSound.m' files from the 
       PSYCHOACOUSTICS toolbox (Soranzo & Grassi, 2014) have to be in the 
       same folder as this script.  

    The most important output file has the prefix 'aMMN_stimulation_SEQUENCE_'.
    It contains a definition of each sound in the 1450-sounds sequence. 

    In order to use this script for an EEG experiment, functioning trigger
    code has to be added. Set up whatever trigger system you are using by
    replacing the current code on lines 166-176. The triggers for the sounds 
    have to be sent on lines 204/205, 217/218, and 226/227. 

    ---------------------------------------------------------------------

    If you use this script in your own work, please cite: 

    Rutiku, R., Fiscone, C., Massimini, M., & Sarasso, S. 
    (currently under 2nd review). 
    Assessing MMN and P3b sensitivity within-individual ñ a comparison 
    between the local-global paradigm and two specialized oddball sequences.

    AND

    Soranzo, A., & Grassi, M. (2014). PSYCHOACOUSTICS: a comprehensive MATLAB 
    toolbox for auditory testing. Frontiers in psychology, 5, 712.

    AND 

    The recommended publications to cite for psychtoolbox: 
    http://psychtoolbox.org/credits

%}

%% Set up the session

clear; clc;                                 % clear MATLAB environment

timestamp = datestr(now);                   % timestamp for data files
timestamp = regexprep(timestamp, '-','');
timestamp = regexprep(timestamp, ' ','');
timestamp = regexprep(timestamp, ':','');

diary(['aMMN_diary_', timestamp])           % a record of the Command Window output

%% Generate the sounds

sf        = 48000;                          % sampling frequency (Hz); make sure that's what you sound card is running
sDuration = 50;                             % standard duration (ms)
sFreqs    = [523, 1046, 1569];              % standard chord frequencies (Hz)
sAmps     = [1, 1/2, 1/4];                  % normalized chord amplitude modulations (dB)

% Cosine ramp; from http://www.h6.dion.ne.jp/~fff/old/technique/auditory/matlab.html
% works much better than GenerateEnvelope 
dr = 0.005;                                 % 5 ms ramps
nr = floor(sf * dr);
CSramp = sin(linspace(0, pi/2, nr));
CSramp = [CSramp, ones(1, floor(sf * sDuration/1000) - nr * 2), fliplr(CSramp)];


% The STANDARD TONE
sTone = GenerateTone(sf, sDuration, sFreqs, sAmps);
sTone = sTone .* CSramp';                   % on-off ramps 
sTone = [sTone'; sTone'];                   % two identical sound channels for both ears
%sound(sTone, sf); % 16-bit default
%plot(sTone')

% The FREQUENCY DEVIANTS
FHdFreqs = [609, 609*2, 609*3];             % higher than the standard sound
FHdTone  = GenerateTone(sf, sDuration, FHdFreqs, sAmps);
FHdTone  = FHdTone .* CSramp'; 
FHdTone = [FHdTone'; FHdTone'];
%sound(FHdTone, sf); %plot(FHdTone')
FLdFreqs = [450, 450*2, 450*3];             % lower than the standard sound
FLdTone  = GenerateTone(sf, sDuration, FLdFreqs, sAmps);
FLdTone  = FLdTone .* CSramp';
FLdTone = [FLdTone'; FLdTone'];
%sound(FLdTone, sf); %plot(FLdTone')

% The INTENSITY DEVIANT
IdTone = GenerateTone(sf, sDuration, sFreqs, sAmps);
IdTone = IdTone .* CSramp';
IdTone = AttenuateSound(IdTone, - 15);      % softer by 15 dB
IdTone = [IdTone'; IdTone'];
%sound(IdTone, sf); %plot(IdTone')

% The DURATION DEVIANT
DdDuration = sDuration-27;         
DdTone     = GenerateTone(sf, DdDuration, sFreqs, sAmps);
Dddr       = dr;                            % new ramps for shorter sound
Ddnr       = floor(sf * Dddr);
DdCSramp   = sin(linspace(0, pi/2, Ddnr));
DdCSramp   = [DdCSramp, ones(1, floor(sf * DdDuration/1000) - Ddnr * 2), fliplr(DdCSramp)];
DdTone     = DdTone .* DdCSramp';
DdTone = [DdTone'; DdTone'];
%sound(DdTone, sf); %plot(DdTone')

% The LOCATION DEVIANTS
LLdTone = GenerateTone(sf, sDuration, sFreqs, sAmps,[0 0 0],-0.5,-1);
LLdTone = LLdTone' .* [CSramp;CSramp];      % coming from the left
%sound(LLdTone, sf); %plot(LLdTone')
LRdTone = GenerateTone(sf, sDuration, sFreqs, sAmps,[0 0 0],0.5,1);
LRdTone = LRdTone' .* [CSramp;CSramp];      % coming from the right
%sound(LRdTone, sf); %plot(LRdTone')

% sound(LRdTone, sf); pause(0.25); sound(LLdTone, sf); 

%% Generate the stimulation sequence

stims = {};
stims{1,1} = sTone;     % 1. index == standard
stims{1,2} = FHdTone;   % 2. index == frequency deviant
stims{2,2} = FLdTone; 
stims{1,3} = IdTone;    % 3. index == intensity deviant
stims{1,4} = DdTone;    % 4. index == duration deviant
stims{1,5} = LLdTone;   % 5. index == location deviant
stims{2,5} = LRdTone;


% Sequence rules from Pakarinen et al. (2007):
% 1. Every other tone is a standard (50 %)
% 2. Each deviant type occurs once in an array of four successive deviants
% 3. Two successive deviants are never the same

SEQUENCE = zeros(1,720);
SEQUENCE(1:4) = randperm(4)+1;
for ar = 5:4:720
    ar_order = randperm(4)+1;
    if SEQUENCE(ar-1) ~= ar_order(1)
        SEQUENCE(ar:ar+3) = ar_order;
    else
        SEQUENCE(ar:ar+3) = fliplr(ar_order);
    end
end
% length(find(SEQUENCE==2)); find(diff(SEQUENCE) ==0) % Sanity checks

ROWS = ones(1,720);
% high/low frequency deviants
HighLow = [ones(1,90), ones(1,90)*2];
HighLow = HighLow(randperm(numel(HighLow)));
ROWS(find(SEQUENCE==2)) = HighLow;
% left/right location deviants
LeftRight = [ones(1,90), ones(1,90)*2];
LeftRight = LeftRight(randperm(numel(LeftRight)));
ROWS(find(SEQUENCE==5)) = LeftRight;

% SAVE!!!!
save(['aMMN_stimulation_SEQUENCE_', timestamp], 'SEQUENCE', 'ROWS')

%% Initiate triggers 

% %Open Device
% UseTriggerBox('Open','C:\Users\Neuro\Desktop\auditory_paradigms\TriggerBoxNET_Net3.dll','GeMSTR13007');
% disp('Open device...');
% 
% % Test trigger
% UseTriggerBox('Send',1,0,0,0,0,0,0,0);
% pause(0.5); % due to 10 ms output latency
% UseTriggerBox('Send',0,0,0,0,0,0,0,0);

%% Run the Optimum-1 sequence

% Running on PTB-3? Abort otherwise.
AssertOpenGL;
% Perform basic initialization of the sound driver:
InitializePsychSound(1);

% Open the default audio device [], with default mode [] (==Only playback),
% and a required latencyclass of 1 == friendly low-latency mode
% a frequency of sf and 2 sound channels
% This returns a handle to the audio device:
pahandle = PsychPortAudio('Open', [], [], 1, sf, 2);


pause(2.0) % short delay before the sequence begins 

t0z = zeros(1,10);
t1z = zeros(1,720);
t2z = zeros(1,720);
tic;
for first = 1:10
    
    PsychPortAudio('FillBuffer', pahandle, stims{1,1});
    
    % start it immediately (0) and wait for the playback to start, return onset timestamp.
    t0z(first) = PsychPortAudio('Start', pahandle, [], 0, 1);
    
    %UseTriggerBox('Send',1,0,0,0,0,0,0,0);
    %UseTriggerBox('Send',0,0,0,0,0,0,0,0);
    
    % Pause for 500 ms ISI + 50 ms stimulus time (-10 ms to account for output latency) 
    % NB! Check output latency and adjust ISI if necessary
    java.lang.Thread.sleep(540);  % in ms
    
end

for trl = 1:720
    
    PsychPortAudio('FillBuffer', pahandle, stims{ROWS(trl), SEQUENCE(trl)});
    t1z(trl) = PsychPortAudio('Start', pahandle, [], 0, 1);
    
    %UseTriggerBox('Send',1,0,0,0,0,0,0,0);
    %UseTriggerBox('Send',0,0,0,0,0,0,0,0);
    
    java.lang.Thread.sleep(540);     % NB! Check output latency and adjust ISI if necessary
    
    
    PsychPortAudio('FillBuffer', pahandle, stims{1,1});
    t2z(trl) = PsychPortAudio('Start', pahandle, [], 0, 1);
    
    %UseTriggerBox('Send',1,0,0,0,0,0,0,0);
    %UseTriggerBox('Send',0,0,0,0,0,0,0,0);
    
    java.lang.Thread.sleep(540); 
    
end
toc; % Takes 802.8 seconds
PsychPortAudio('Close')
%UseTriggerBox('Close')
diary off

% SAVE!!!!
save(['aMMN_critical_events', timestamp], 't0z', 't1z', 't2z')

%% Save a backup of the entire Workspace 
%(can be removed after verifying that the procedure works properly)

save(['aMMN_backup_dump_', timestamp])

