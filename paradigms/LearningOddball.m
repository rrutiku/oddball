
%{
    Created July 2016
	@author: Renate Rutiku
	contact: renate.rutiku@gmail.com

    ---------------------------------------------------------------------
    ------------ Learning-oddball paradigm for P3b assessment -----------
    ---------------------------------------------------------------------

    This scrips runs the auditory learning-oddball paradigm (Jongsma et al., 
    2006, 2013) with minor adjustments to the stimuli.
    
    In order for the script to work:
    1. Psychtoolbox (http://psychtoolbox.org/) has to be installed.
    2. The 'GenerateTone.m' file from the PSYCHOACOUSTICS toolbox (Soranzo 
       & Grassi, 2014) has to be in the same folder as this script.   

    NB! Note that this is an active oddball sequence. Therefore, participants 
    should be instructed to pay attention to the rare deviant sounds and
    count them, for example.

    The most important output file has the prefix 'aP3b_stimulation_SEQUENCE_'.
    It contains a definition of each sound in the 682-sounds sequence. 

    In order to use this script for an EEG experiment, functioning trigger
    code has to be added. Set up whatever trigger system you are using by
    replacing the current code on lines 113-123. The triggers for the sounds 
    have to be sent on lines 149/150 and 161/162. 

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

diary(['aP3b_diary_', timestamp])           % a record of the Command Window output

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

% The FREQUENCY DEVIANT
FHdFreqs = [609, 609*2, 609*3];             % higher than the standard sound                             
FHdTone  = GenerateTone(sf, sDuration, FHdFreqs, sAmps);
FHdTone  = FHdTone .* CSramp'; 
FHdTone = [FHdTone'; FHdTone'];
%sound(FHdTone, sf); %plot(FHdTone')

%% Generate the stimulation sequence

stims = {};
stims{1,1} = sTone;         % 1. index == standard
stims{1,2} = FHdTone;       % 2. index == higher frequency deviant

SEQUENCE = [];
idis = repmat([2:5 7:10], 1, 6);
idis = idis(randperm(numel(idis))); 

% LEARNING ODDBALL SEQUENCE
count = 1;
for r = 1:8:48
    randfirst = idis(r:r+7);
    for i = 1:8
        SEQUENCE(count:count+randfirst(i)-1) = 1;
        count = count+randfirst(i);
        SEQUENCE(count) = 2;
        count = count+1;
    end
    syssec = repmat([ones(1,6),2],1,8);
    syssec(7:7:56) = 2;
    SEQUENCE(count:count+length(syssec)-1) = syssec;
    count = count+length(syssec);
end
% length(find(SEQUENCE==2))/length(SEQUENCE); % Sanity checks

pauses = zeros(1,length(SEQUENCE));
for p = 1:length(SEQUENCE)
    pauses(p) = randsample(800:50:1200,1);
end

% SAVE!!!!
save(['aP3b_stimulation_SEQUENCE_', timestamp], 'SEQUENCE', 'pauses')

%% Initiate triggers 

% %Open Device
% UseTriggerBox('Open','C:\Users\Neuro\Desktop\auditory_paradigms\TriggerBoxNET_Net3.dll','GeMSTR13007');
% disp('Open device...');
% 
% % Test trigger
% UseTriggerBox('Send',1,0,0,0,0,0,0,0);
% pause(0.5); % due to 10 ms output latency
% UseTriggerBox('Send',0,0,0,0,0,0,0,0);

%% Run the learning-oddball sequence

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
tz = zeros(1,length(SEQUENCE));
tic;
for first = 1:10
    
    PsychPortAudio('FillBuffer', pahandle, stims{1,1});
    % start it immediately (0) and wait for the playback to start, return onset timestamp.
    t0z(first) = PsychPortAudio('Start', pahandle, [], 0, 1);
    
    %UseTriggerBox('Send',1,0,0,0,0,0,0,0);
    %UseTriggerBox('Send',0,0,0,0,0,0,0,0);
    
    java.lang.Thread.sleep(randsample(800:50:1200,1))
    
end

for trl = 1:length(SEQUENCE)
    
    PsychPortAudio('FillBuffer', pahandle, stims{1, SEQUENCE(trl)});
    tz(trl) = PsychPortAudio('Start', pahandle, [], 0, 1);
    
    %UseTriggerBox('Send',1,0,0,0,0,0,0,0);
    %UseTriggerBox('Send',0,0,0,0,0,0,0,0);
    
    java.lang.Thread.sleep(pauses(trl))
    
end
toc; % Takes ca. 808 seconds
PsychPortAudio('Close')
%UseTriggerBox('Close')
diary off

% SAVE!!!!
save(['aP3b_critical_events', timestamp], 't0z', 'tz')

%% Save a backup of the entire Workspace 
%(can be removed after verifying that the procedure works properly)

save(['aP3b_backup_dump_', timestamp])
