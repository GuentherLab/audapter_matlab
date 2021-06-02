function [ScanTrigger] = scanTrigger(timeLeft, y, fs, sinePlayer, curTrialTimer)

    tic
    while toc < timeLeft
    end
    
    
    
    y = y - mean(y);

vmax = max(abs(y));
y = y / vmax * 0.98;

%player = audioplayer(y, fs);
ScanTrigger = toc(curTrialTimer);
fprintf('\nScanner Triggered\n')
play(sinePlayer);

dur = length(y) / fs;
pause(dur);

end

% %written by Jenn Segawa
% %settings for testing in Bay 1
% %Audapter_test_2(50, 1, 1 / 50 / 2)
% 
%  Audapter(3,'wgfreq',1000,0); %trigger tone
%  Audapter(3,'wgamp',.2,0);
% % freq = 50;
% % amp = 1;
% % length = 1 / 50 / 2;
% % 
% % dur_0 = 0.05; 
% % dur_1 = length;
% % 
% % fs = 48000;
% % w = zeros(1, dur_0 * fs);
% % 
% % t = 0 : (1 / fs) : dur_1;
% % s = amp * sin(freq * 2 * pi * t);
% % w = [w, s];
% % % 
% % % 
% % % %generate tone------------------------------------
% % % %Sine wave (pure tone) generator. Plays a continuous pure tone of ...
% % % %frequency wgfreq (Hz), amplitude wgamp and initial time wgtime, 
% % % %that is, wgamp×sin(wgfreq×(t+wgtime)). No ramp is imposed.
% % % 
% % % 
% % % w = [w, 0 * ones(1, Audapter('getMaxPBLen') - size(w, 2))];
% % % 
% % % %%%check this
% % Audapter('setParam', 'datapb', w);
% % Audapter('reset');
% % 
% % Audapter('playWave');
% % pause(1);
% % Audapter('stop');
% % % fprintf('Scanner Triggered')

%Audapter('reset');