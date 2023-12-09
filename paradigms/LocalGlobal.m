
%{
    Created July 2016
	@author: Renate Rutiku
	contact: renate.rutiku@gmail.com

    ---------------------------------------------------------------------
    --------- Local-global paradigm for MMN and P3b assessment ----------
    ---------------------------------------------------------------------

    This scrips runs the auditory local-global paradigm (Bekinschtein et al., 
    2009) with minor adjustments.
    
    In order for the script to work:
    1. Psychtoolbox (http://psychtoolbox.org/) has to be installed.
    2. The 'GenerateTone.m' file from the PSYCHOACOUSTICS toolbox (Soranzo 
       & Grassi, 2014) has to be in the same folder as this script. 

    NB! Note that this is an active oddball sequence. Therefore, participants 
    should be instructed to pay attention to the rare global deviant quintlets 
    and count them for each block separately, for example.

    The most important output file has the prefix 'GlobLoc_stimulation_SEQUENCE_'.
    It contains a definition of each sound in the entire sequence. 

    In order to use this script for an EEG experiment, functioning trigger
    code has to be added. Set up whatever trigger system you are using by
    replacing the current code on lines 137-147. The triggers for the sounds 
    have to be sent on lines 174/175. 

    ---------------------------------------------------------------------

    If you use this script in your own work, please cite: 

    Rutiku, R., Fiscone, C., Massimini, M., & Sarasso, S. 
    (currently under 2nd review). 
    Assessing MMN and P3b sensitivity within-individual – a comparison 
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

timestamp = datestr(now);
timestamp = regexprep(timestamp, '-','');
timestamp = regexprep(timestamp, ' ','');
timestamp = regexprep(timestamp, ':','');

diary(['GlobLoc_diary_', timestamp])        % a record of the Command Window output

%% Generate the sounds

sf        = 48000;                          % sampling frequency (Hz); make sure that's what you sound card is running
sDuration = 50;                             % standard duration (ms)
sFreqs    = [350,700,1400];                 % standard chord frequencies (Hz)
sAmps     = [1, 1/2, 1/4];                  % normalized chord amplitude modulations (dB)

% Cosine ramp; from http://www.h6.dion.ne.jp/~fff/old/technique/auditory/matlab.html
% works much better than GenerateEnvelope 
dr = 0.005; % 5 ms ramps
nr = floor(sf * dr);
CSramp = sin(linspace(0, pi/2, nr));
CSramp = [CSramp, ones(1, floor(sf * sDuration/1000) - nr * 2), fliplr(CSramp)];


% The STANDARD TONE
sTone = GenerateTone(sf, sDuration, sFreqs, sAmps);
sTone = sTone .* CSramp';                   % on-off ramps 
sTone = [sTone'; sTone'];                   % two identical sound channels for both ears
%sound(sTone, sf); % 16-bit default
%plot(sTone')

% The FREQUENCY DEVIANT
FHdFreqs = [500,1000,2000];
FHdTone  = GenerateTone(sf, sDuration, FHdFreqs, sAmps);
FHdTone  = FHdTone .* CSramp'; 
FHdTone = [FHdTone'; FHdTone'];
%sound(FHdTone, sf); %plot(FHdTone')

%% Generate the stimulation sequence

stims = {};
stims{1,1} = sTone;     % 1. index == standard
stims{1,2} = FHdTone;   % 2. index == higher frequency deviant

Blocks = repmat(1:4,1,2);
Blocks = Blocks(randperm(numel(Blocks)));

SEQUENCE = {};

for bl = 1:numel(Blocks)
    
   SEQUENCE{bl} = [];
   
   % the standards and deviants for each block type; 
   % numbers will reference to stims struct later
   if Blocks(bl) == 1
       standard = [1 1 1 1 1];
       deviant = [1 1 1 1 2];
   elseif Blocks(bl) == 2
       standard = [2 2 2 2 2];
       deviant = [2 2 2 2 1];
   elseif Blocks(bl) == 3
       standard = [1 1 1 1 2];
       deviant = [1 1 1 1 1];
   elseif Blocks(bl) == 4
       standard = [2 2 2 2 1];
       deviant = [2 2 2 2 2];       
   end
   
   % introduce the standard before the first deviant
   SEQUENCE{bl} = [SEQUENCE{bl}, repmat(standard, 10,1)];
   
   % how many deviants and 4 times more standards == 20/80
   nrDeviants = randsample(22:30,1);
   nrStandard = 4*nrDeviants;
   
   % order of stimuli
   b2rand = [ones(nrStandard,1); ones(nrDeviants,1)*2];
   while ~isempty(findstr(b2rand',[2 2]))   
       b2rand = b2rand(randperm(numel(b2rand)));
   end
   
   % add to the SEQUENCE
   for i = 1:numel(b2rand)
       if b2rand(i) == 1
           SEQUENCE{bl} = [SEQUENCE{bl}; standard];
       elseif b2rand(i) == 2
           SEQUENCE{bl} = [SEQUENCE{bl}; deviant];
       end
   end
   
end
% sum(SEQUENCE{1,2}(:,5) ~= SEQUENCE{1,2}(1,1)) % sanity checks

% SAVE!!!!
save(['GlobLoc_stimulation_SEQUENCE_', timestamp], 'SEQUENCE')

%% Initiate triggers 

% %Open Device
% UseTriggerBox('Open','C:\Users\Neuro\Desktop\auditory_paradigms\TriggerBoxNET_Net3.dll','GeMSTR13007');
% disp('Open device...');
% 
% % Test trigger
% UseTriggerBox('Send',1,0,0,0,0,0,0,0);
% pause(0.5); % due to 10 ms output latency
% UseTriggerBox('Send',0,0,0,0,0,0,0,0);

%% Run the local-global sequence 

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

tz = zeros(1,6000);
count = 1;

tic;
for bl = 1:numel(Blocks)
    for trl = 1:size(SEQUENCE{bl},1)

        for i = 1:5
            PsychPortAudio('FillBuffer', pahandle, stims{1,SEQUENCE{bl}(trl,i)});
            tz(count) = PsychPortAudio('Start', pahandle, [], 0, 1);
            %UseTriggerBox('Send',1,0,0,0,0,0,0,0);
            %UseTriggerBox('Send',0,0,0,0,0,0,0,0);
            count = count+1;
            java.lang.Thread.sleep(140);  % in ms     % NB! Check output latency and adjust ISI if necessary
        end
        
        java.lang.Thread.sleep(randsample(700:50:1000,1));

    end
    
    if bl == 4
        [secs, keyCode, deltaSecs] = KbStrokeWait;
    else
        pause(5.0) % short delay before the next block begins
    end
end
toc; % Takes ca. 2000 seconds

PsychPortAudio('Close')
%UseTriggerBox('Close')
diary off

% SAVE!!!!
save(['GlobLoc_critical_events', timestamp], 'tz')

%% Save a backup of the entire Workspace 
%(can be removed after verifying that the procedure works properly)

save(['GlobLoc_backup_dump_', timestamp])
