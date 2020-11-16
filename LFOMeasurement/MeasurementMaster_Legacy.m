%% This is a LFO measurment master file
% The file generates a test signal consisting of a series of chirps
% For carrying out LFO measurements on digital pedals, set 'Mode = 1'
% For carrying out LFO measurements on a real pedal, set 'Mode = 2'
% Mode = 2;
% % Pedal 1 = sweet tone phaser
% Pedal = 1;
% % Initial frequency estimation technique
% Freq_Est = 1;
% Ped_Name = 'VintagePhaser';
clear all
import LFOAnalyser

addpath('../DigitalEffectsPedals')

Device = 'PhaserSweetTone';
Digi = 1;
SNR = 40;
method = 1;
ch_ty = 1;
ch_len = 19;
ch_spc = 20;
T = 10;
fs = 44100;
rate = 0.3;
smo_f = 5;
hold off

if ~exist('pnum')
    LFOAnly = LFOAnalyser(Device, Digi);

    [LFOAnly, pnum] = LFOAnly.SigGen(method, ch_len, ch_spc, 20, 20000, T, fs);
    LFOAnly = LFOAnly.SigProc(pnum, rate, SNR);
    LFOAnly = LFOAnly.SpecExtract(pnum, 0.1);
    LFOAnly = LFOAnly.PrelimAnly(pnum);
end


fmin = round(0.9*LFOAnly.Measurements{1,'notch_min'});
fmax = round(1.05*LFOAnly.Measurements{1,'notch_max'});
[LFOAnly, num] = LFOAnly.SigGen(method, ch_len, ch_spc, fmin, fmax, T, fs);

LFOAnly = LFOAnly.SigProc(num, rate, SNR);
LFOAnly = LFOAnly.SpecExtract(num, 0.1);
if size(LFOAnly.Measurements{num,'Measured_LFO'}{1,1}) < 1
    LFOAnly = LFOAnly.LFOTrack(num);
end
LFOAnly = LFOAnly.LFOTrack2(num, smo_f);

LFOAnly = LFOAnly.RectSineFit(num);

params = LFOAnly.Measurements{num,'Sine_Params'}{1,1};
LFO_Measured = LFOAnly.Measurements{num,'Measured_LFO'}{1,1};

pred_LFO = LFOAnalyser.rectSineGen(params(1),params(2),params(3),params(4)*pi, LFO_Measured(:,2));

error1 = mean((LFO_Measured(:,1) - pred_LFO).^2);

plot(LFOAnly.Measurements{end,'Measured_LFO'}{1,1}(:,1));
hold on
plot(pred_LFO)

pred_LFO = LFOAnalyser.rectSineGenPow...
    (params(1),params(2),params(3),params(4)*pi, LFO_Measured(:,2), 2);


params = LFOAnalyser.LFOCorr(params, LFO_Measured(:,1));

pred_LFO = LFOAnalyser.rectSineGen(params(1),params(2),params(3),params(4)*pi, LFO_Measured(:,2));

power = LFOAnalyser.LFOReshapeFactor(params(1),params(2),params(3),params(4)*pi,...
                               LFO_Measured(:,1), LFO_Measured(:,2));

pred_LFO = LFOAnalyser.rectSineGenPow(params(1),params(2),params(3),params(4)*pi, LFO_Measured(:,2), power);
pred_LFO = LFOAnalyser.rectSineGenPow(params(1),params(2),params(3),params(4)*pi, LFO_Measured(:,2), power);

    

% LFOAnly = LFOAnly.PrelimAnly(num);



% 
% [LFO_Measured, LFO_True] = LFOMeasurer...
%                         (Device, Digi, ch_ty, ch_len, ch_spc, T, fs);


% pedal_name = functions(DigitalEffects{Pedal,1}).function;

% method = [1,1,1,2,2,2];
% ch_len = [9.5,19.5,29.5,9.5,19.5,29.5];
% ch_spc = [10,20,30,10,20,30];
% f_st   = [700,20,20,20,20,20];
% f_en   = [2000,2000,2000,2000,2000,2000];
% T      = [10,10,10,10,10,10];
% fs     = [44100,44100,44100,44100,44100,44100];


% switch Mode
%     case 1
% 
%         for m = 1
%             [test_sig, ch_i, n] = TestSignalGen(method(m), ch_len(m)...
%                                         ,ch_spc(m),f_st(m), f_en(m), T(m), fs(m));
% 
%             
% 
%             switch Pedal
%                 case 1
%                     speed = 80;
% 
%                     [proc_sig, LFO_real] =...
%                         PhaserSweetTone(test_sig, fs(m), speed, 0, 0);
% 
%                     [LFO_meas, tAx] = TestSignalAnl...
%                         (proc_sig, ch_i, method(m), f_st(m), f_en(m), fs(m), n);
% 
%                     plot(tAx, LFO_meas)
%                     disp('Hows that LFO looking? (press any key to continue)')
%                     pause
% 
%                     init_f = InitFreqEst(LFO_meas, tAx);
% 
%                     [A, B, C, f] = RectSineFit(LFO_Meas, tAx, init_f);
% 
%                     f_tru = 0.069*exp(0.040*speed);
%                     error = 0.069*exp(0.040*speed) - f
% 
% 
% 
%             end
% 
%         %     LFO_pred(:,1) = LFO_pred(:,1)./max(LFO_pred(:,1));
% 
%             tAx = 0:1/fs(m):10 - 1/fs(m);
%             plot(tAx, LFO_meas(1:length(tAx),4))
%             hold on
%             plot(LFO_pred(:,2), LFO_pred(:,1), 'x')
% 
%         end
%     case 2
%         
% %         stage = 1;
%         method = 1;
%         ch_len = 9.5;
%         ch_spc = 10;
%         f_st = 20;
%         f_en = 2000;
%         T = 10;
%         fs = 44100;
%         
%         switch stage
%             case 1
%             
%                 [test_sig, ch_i, n] = TestSignalGen(method, ch_len...
%                                         ,ch_spc,f_st, f_en, T, fs);
%                                             
%                 test_sig = test_sig/max(abs(test_sig));
% 
%                 audiowrite(strcat(Ped_Name,'/ChirpTrainInput.wav')...
%                     , [test_sig,test_sig], 44100)
%                 
%                 stage = 2;
%                 
%             case 2
%                 % Read the processed audio signal
%                 proc_sig = audioread(strcat(Ped_Name,'/ChirpTrainOutput.wav'));
%                 loopback_sig = proc_sig(:,2);
%                 proc_sig = proc_sig(:,1);
%                 [r, lags] = xcorr(test_sig, loopback_sig);
%                 [~, i] = max(r);
%                 proc_sig = proc_sig(-lags(i):end);
%                 
%                 if length(proc_sig) > length(test_sig)
%                     disp('The recorded signal is longer than the test signal...')
%                     pause
%                 end
%                 
%                 plot(test_sig(ch_i(1)- 100:ch_i(2) + 300))
%                 hold on
%                 plot(proc_sig(ch_i(1) - 100:ch_i(2) + 300))
%                 hold off
%                 disp('The two signals should be well alligned')
%                 pause
%                 
%                 stage  = 3;
%                 
%             case 3
%                 
%                 [LFO_meas, tAx] = TestSignalAnl...
%                     (proc_sig, ch_i, method, f_st, f_en, 44100, n);
%                 
%                 plot(tAx, LFO_meas)
%                 disp('Hows that LFO looking? (press any key to continue)')
%                 pause
%                 
%                 
% 
%         end
%         
%         
% end
% 
% function [init_f] = InitFreqEst(LFO, tAx, method)
% 
%     [~,locs] = findpeaks(LFO);
%     T = tAx(locs(2)) - tAx(locs(1));
%     init_f = 1/T;
%     
% %     [y, tAx] = resample(LFO_Measured(:,1),LFO_Measured(:,2));
% %     fs = 1/(tAx(2) - tAx(1));
% %     
% %     % Desired frequency resolution
% %     f_steps = 100000;
% %     fAx_full = 0:fs/f_steps:fs - 1/f_steps;
% %     f_res = fAx_full(2) - fAx_full(1);
% %     
% %     ax_st = floor((0.7*init_f/f_res) + 1);
% %     ax_en = ceil((1.3*init_f/f_res) + 1);
% %     fAx = fAx_full(ax_st:ax_en);
% %     
% %     spec = fft(y, f_steps);
% %     spec = spec(ax_st:ax_en);
% % 
% %     [~,k] = max(abs(spec));
% %     init_f = fAx(k);
%     
% end