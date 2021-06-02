function GASP_Formant_Reflex()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prior scripts utlizied
% PDGAP scripts from Defne Abur as at 8/11/2018 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GASP Audapter Formant Adaptation Subject Run Script
% Final edit : 11/03/2019
% Author : Hasini R Weerathunge
%
% 1) converted 24 control trials to shift down perturbations
% 2) changed uttenrences to 'bed','Ed','head' so shift up and shift down is possible
% 3) changed word randomization to have equal number of stimuli in each word perturbed
% 4) Added vocal onset time triggered reflex activity (online reflex perturbations)
% ALWAYS RUN GASP_initialization.m prior to running script and GASP_termination.m after subject run is complete. 
% Makes sure that this version of audapter mex file is isolated to GASP study only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart = tic;
close all;
clc;
cd('C:\GASP\speechres\audapter_matlab\mcode');

%% define experimental parameters
calib_check = input('Has this microphone position been calibrated? (0 = no,1 = yes):');

if calib_check ==0
    errordlg('Run the IntensityCalibration script');
    return
end

cue_check = input('Did you check the CueMix Fx configuration is correct? (0 = no,1 = yes):');
if cue_check == 0
    errordlg('Check CueMix setting!');
    return
end

cue_save = input('Did you save the CueMix Fx configuration to the subject folder? (0 = no,1 = yes):');
if cue_save == 0
    errordlg('Save CueMix setting!');
    return
end


prompt={'Subject ID:',...
        'Experimenter initials',...
        'Session ID:',...
        'group:',...
        'Gender ("male" or "female")','noise dB change'};
name='Subject Information';
numlines=1;
defaultanswer={'YCGASP','','Formant_Reflex','Control','female','-17'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
subjectID = answer{1};
experimenter = answer{2};
session   = answer{3};
group     = answer{4};
Gender    = answer{5};
dBchange  = str2num(answer{6});

dBratio = 10^(dBchange/20);

num_trials = questdlg('Practice or Full?','Length','Practice','Full','Full') ;
switch num_trials
    case 'Practice'
        numBlocks = 1;
        num_trials = 'Practice';
        load('pertInfo_practise.mat');
        vis =1 ;%HW 08/23/2018 visualization of f1 shifts
    case 'Full'
        numBlocks = 12;
        num_trials = 'Full';
        load('pertInfo.mat'); % make a randomization for this
        vis =0 ;%HW 08/23/2018 visualization of f1 shifts
end

% create the folder path to save the data files
baseFilename = ['data\' group, '\', subjectID, '\', session, '\',num_trials,'\'];

% check if the foler exists (to avoid overwrite)
if (exist(baseFilename,'dir'))
    h = errordlg('THE DIRECTORY EXISTS!');
    return;
else
    mkdir(baseFilename);
    while exist(baseFilename,'dir') ~= 7
            mkdir(baseFilename);
    end
end


%% initialization of the parameters
% Perturbed words
mainWords = {'BED'; 'ED'; 'HEAD' }; % HW 8/11/2018

trialDuration = 3;
trialInBlock = 9; % blocks 9 cannot change
PertInBlock = 4; % blocks 2 cannot change
maxGain = .3; % 30% upward
minGain = .3; % 30% downwards % HW 8/11/2018 only phase change required

Pert = zeros(1,numBlocks*trialInBlock);
phaseShift = zeros(1,numBlocks*trialInBlock);
numTrials = length(Pert);
rng('shuffle');

for i = 1:numTrials
    Words{i,1} = mainWords(pertInfo(i,2));
    Pert(1,i) = pertInfo(i,1);
    if Pert(1,i) == 2
        phaseShift(1,i) = 1;
        Pert(1,i) = minGain;
    elseif Pert(1,i) == 1
        Pert(1,i) = maxGain;
    end
end

Fs = 44100;

% calculate gain value for every trial
delay = randi([1, 10], 108); %OF NOTE: these yield actual delays of .5-1 due to equipment latency

%% Initialize Audapter
addpath c:\GASP\speechres\commonmcode;
cds('audapter_matlab');
audioInterfaceName = 'ASIO4ALL';%'MOTU MicroBook';%
sRate = 48000;  % Hardware sampling rate (before downsampling)
downFact = 3;
frameLen = 96;  % Before downsampling
noiseWavFN = 'SSN_ampChunk.wav';
Audapter('deviceName', audioInterfaceName);
Audapter('setParam', 'downFact', downFact, 0);
Audapter('setParam', 'sRate', sRate / downFact, 0);
Audapter('setParam', 'frameLen', frameLen / downFact, 0);
Audapter('ost', '', 0);
Audapter('pcf', '', 0);
params = getAudapterDefaultParams(lower(Gender));
params.f1Min = 0;
params.f1Max = 5000; %HW 08/13/2018
params.f2Min = 0;
params.f2Max = 5000;
params.pertF2 = linspace(0, 5000, 257);
params.pertAmp = 0.0 * ones(1, 257);
params.pertPhi = 0.0 * pi * ones(1, 257);
params.bTrack = 1;
params.bShift = 1;
params.bRatioShift = 1;
params.bMelShift = 0;
maxPBSize = Audapter('getMaxPBLen');
check_file(noiseWavFN);
[w, fs] = audioread(noiseWavFN);
if fs ~= params.sr * params.downFact
    w = resample(w, params.sr * params.downFact, fs);
end
if length(w) > maxPBSize
    w = w(1 : maxPBSize);
end

w = w*dBratio; %was set to 0.05

Audapter('setParam', 'datapb', w, 1);
%params.fb3Gain = 0.1;
params.fb = 1; % play speech only ; fb =3 play speech and noise - OST controls gain level (1)switch on noise (0) switch off noise
AudapterIO('init', params);
%%
%% Selecting booth and figure dimensions - HW - 08/13/2018
   figure1 = figure('NumberTitle','off','Color',[0 0 0],'Position',[1680 0 1920 1080],'MenuBar','none');

  
textBox1 = annotation(figure1,'textbox',...
    [0.3 0.45 0.01 0.1],... 
    'Color',[1 1 1],...
    'String','',...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontSize',130,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'BackgroundColor',[0 0 0],...
    'Visible','on');


drawnow;
pause(5)
set(textBox1,'String','READY');
drawnow;
pause(5)

%% Present words and collect the results
%vis =1 ;%HW 08/23/2018 visualization of f1 shifts

for tr = 1:numTrials+3
    
    if tr == 1 || tr == 2 || tr == 3
        AudapterIO('init',params);
        % start Audapter
        Audapter('reset');
        Audapter('start');
        % wait and receive data
        pause(trialDuration);
        fprintf(1, '\n');
        Audapter('stop');
        pause(1)
    else
    shift = pertInfo(tr-3,1);
    Shift_amount= Pert(tr-3)*((-1)^(phaseShift(tr-3)));  %HW 08/13/2018
    if pertInfo(tr-3,1) == 1% shift down
        params.dScale = 1.25; disp('set to 1.25');
    elseif pertInfo(tr-3,1) == 2
        params.dScale = 0.5;  disp('set to 0.45');
    elseif pertInfo(tr-3,1) == 0
        params.dScale = 0.9; disp('set to 0.9');
    end
        
    params.bgainadapt = 0 ; %set to 1 in C++ code
    
       
    
    fprintf('Trial %d, Shift = %1.3f , Word = %s , delay = %1.3f \n',tr-3,Shift_amount, char(Words{tr-3,1}), delay(tr-3)) %HW 08/13/2018
    %update onset time of perturbation from vocal onset %HW 01/22/2019
    if shift ~= 0
        switch delay(tr-3)
            case 1
              ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_m_1.ost';  
              baseline_delay = 0.3;
            case 2
              ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_m_2.ost';  
              baseline_delay = 0.3625;
            case 3
              ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_m_3.ost';  
              baseline_delay = 0.425;
            case 4
              ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_m_4.ost';
              baseline_delay = 0.4875;
            case 5
              ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_m_5.ost';  
              baseline_delay = 0.55;
            case 6
              ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_m_6.ost';  
              baseline_delay = 0.6125;
            case 7
              ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_m_7.ost';  
              baseline_delay = 0.675;
            case 8
              ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_m_8.ost';  
              baseline_delay = 0.7375;
            case 9
              ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_m_9.ost';  
              baseline_delay = 0.8;
            case 10
              ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_m_10.ost';  
              baseline_delay = 0.8375;
        end
    end
    
    %setting masking noise feedback
    params.fb = 3; % speech + noise -fb3Gain
    params.fb3Gain = 0.1; %default is 0.0
 
     
    switch shift
        case 0 
             pcfFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_n_m.pcf';
             ostFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_n_m.ost';  
             check_file(pcfFN); Audapter('pcf', pcfFN, 0);
        case 1
             pcfFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_su_m.pcf';
             check_file(pcfFN); Audapter('pcf', pcfFN, 0);
        case 2
             pcfFN = 'C:\GASP\speechres\audapter_matlab\mcode\example_data\formant_reflex_pert_sd_m.pcf';
             check_file(pcfFN); Audapter('pcf', pcfFN, 0);
    end
   
    check_file(ostFN);  
    Audapter('ost', ostFN, 0);
    
    set(textBox1,'String',Words{tr-3,1});
    AudapterIO('init', params);
    

    params.rmsThr =0;
    % start Audapter
    Audapter('reset');
    Audapter('start');
    % wait and receive data
    pause(trialDuration);
    fprintf(1, '\n');
    Audapter('stop');
    
    % save the data
    data = AudapterIO('getData');
    
%     if Shift_amount <0
%         Shift_amount = ['minus', num2str(abs(Shift_amount))];
%     else
%         Shift_amount = num2str(Shift_amount);
%     end
    
    
    fileNameData = [baseFilename 'Trial_' num2str(tr-3) '_' char(Words{tr-3,1}) '_' num2str(Shift_amount) '.mat']; %HW 08/13/2018
    save(fileNameData,'data');
    
    %debug visualiation
    if vis ==1
    visualize(data,sRate) 
    end
    
     ExpSave = [baseFilename 'ExperimenterInitials.mat']; %makes a new .mat file to save 
     ExperimenterInitials.experimenter = experimenter; %defines saving of user input mic trim
     save(ExpSave,'ExperimenterInitials');%saves data
%     CueMixSave = [baseFilename 'CueMixSettings.mat']; %makes a new .mat file to save 
%     CueMixSettings.cuemixMic = cuemixMic; %defines saving of user input mic trim
%     CueMixSettings.cuemixMainout = cuemixMainout; %defines saving of user input main out value
%     save(CueMixSave,'CueMixSettings');%saves data
    
    set(textBox1,'String',' ')
    
    % jitter the between trial pause
    jitterPause = [1,1.5,2];
    tempVarRand = randperm(3);
    pause(jitterPause(tempVarRand(2)));
    end

end
%% Saving and cleaning up
close all;
elapsed_time = toc(tstart)/60;
%clc;
fprintf('Total Time : %f\n',elapsed_time);
end


function visualize(data1,sRate)
%% visualize the perturbation if in visualize state- will take some time from the expected ISI
    OST_MULT = 1000;
    fs= data1.params.sr;
    % Visualization: input sound
     fig1=figure('Position', [100, 100, 1400, 600], 'Name', 'Input spectrogram');
     show_spectrogram(data1.signalIn, fs, 'noFig'); hold on;
     frameDur = data1.params.frameLen / data1.params.sr;
     tAxis = 0 : frameDur : frameDur * (size(data1.fmts, 1) - 1);
     frameDur = data1.params.frameLen / data1.params.sr;%check what is data
     plot(tAxis, data1.ost_stat * OST_MULT, 'k-');
     %legend({'F1 (original)', 'F2 (oringinal)', sprintf('OST status * %d', OST_MULT)});
     xlabel('Time (s)');
     ylabel('Frequency (Hz)');

    %% Visualization: output sound
    %fig2=figure('Position', [100, 100, 1400, 600], 'Name', 'Output spectrogram');
    show_spectrogram(data1.signalOut,  data1.params.sr, 'noFig');
    plot(tAxis, data1.fmts(:, 1 : 2), 'b');
    plot(tAxis, data1.sfmts(:, 1 : 2), 'k-');
    plot(tAxis, data1.ost_stat * OST_MULT, 'k*');
    legend({'F1 (original)', 'F2 (original)', 'F1 (shifted)', 'F2 (shifted)', sprintf('OST status * %d', OST_MULT)});
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');

    drawnow;
    pause;
    close(fig1);
    
end
