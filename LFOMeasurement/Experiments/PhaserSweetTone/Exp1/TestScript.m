addpath('../../../MatlabClasses')
import SignalHolder
import SignalProcessor
import SignalAnalyser
import LFOFitter

%% Expriment 1 - Digital Model of Phaser Sweet Tone
% In this test, the accuracy of the LFO frequency measurement is evaluated
% for varying lengths of test signal, with varying levels of noise
% introduced into the signal to simulate real world test conditions

% Signal Processing Info
PedalName = 'PhaserSweetTone';
Digi = 1;
Rate = [0.1, 0.3, 0.6];
SNR = [10, 20, 30, 40];

% Test Signal Parameters
ch_type = 1;
ch_len = 19;
ch_spc = 20;
% Length, in s
T = [10, 20, 30, 40, 50, 60];
% These are fractions of the test signal to ignore, the % in between the
% two numbers is ignored, i.e [0.35, 0.65] discards from 0.35*T until
% 0.65*T
tfrac = [[0.5,0.5];[0.35, 0.6];[0.25, 0.75];[0.125,0.875]];

fs = 44100;
hold off

% If the object holding the sinals isn't loaded/doesn't exist...
if ~exist('Signals')
    try
        % Try to load it
        disp('Loading Measurments')
        load Measurements.mat
    catch ME
        % Or create it if it can't be loaded
        disp('No Measurements found, creating signals')
        
        % Iterate over the rates being tested
        for r = 1:length(Rate)
            rate = Rate(r);
            
            % Create the Signal Holder and Processor objects
            Signals(r) = SignalHolder();
            ProcSigs(r) = SignalProcessor(PedalName, Digi, rate);
            % Create and Process a preliminary signal
            Signals(r) = Signals(r).SigGen(1, 28, 30, 20, 22050, 2/rate, fs);
            ProcSigs(r) = ProcSigs(r).SigProc(Signals(r).SigGet(1), rate, 40, fs);

            % Analyse the preliminary signal to estiamte LFO freq and
            % parameters for test signal chirp
            AnlySig(r) = SignalAnalyser(Signals(r).Signals, ProcSigs(r).ProcessedSignals);
            AnlySig(r) = AnlySig(r).SpecExtract(0.5, 1);
            AnlySig(r) = AnlySig(r).PrelimAnly(1);
            st_f = 0.9*AnlySig(r).Min_f;
            en_f = 1.1*AnlySig(r).Max_f;

            disp('Preliminary Analysis Complete, Generating Test Signals')

            for t = 1:length(T)
                % Create test signals
                Signals(r) = Signals(r).SigGen(ch_type, ch_len, ch_spc, st_f, en_f, T(t), fs);
                i = size(Signals(r).Signals,1);
                
                % Process at varying SNRs
                for snr = 1:length(SNR)
                    ProcSigs(r) = ProcSigs(r).SigProc(Signals(r).SigGet(i), rate, SNR(snr), fs);
                end
            end

        % Delete some of the unneeded data and save the rest
        Signals(r).Signals = Signals(r).Signals(:, 2:end);
        ProcSigs(r).ProcessedSignals = ProcSigs(r).ProcessedSignals(:,1:4);
        AnlySig(r).Spectrograms = AnlySig(r).Spectrograms(2:end,:);
        AnlySig(r).Signals = 0;
        AnlySig(r).ProcSigs = 0;
        end
        disp('Saving Data')
        save ('Measurements.mat', 'Signals', 'ProcSigs', 'AnlySig') 
    end
end

% Check if an LFOFitter object has been created
if ~exist('LFOFits')
    for r = 1:length(Rate)
        rate = Rate(r);
        % Create LFO Fitter object
        LFOFits(r) = LFOFitter();
        AnlySig(r) = AnlySig(r).Update(Signals(r).Signals, ProcSigs(r).ProcessedSignals);
        
        % Extract Spectrograms and track the LFO
        for sigs = 2:size(ProcSigs(r).ProcessedSignals, 1)
            AnlySig(r) = AnlySig(r).SpecExtract(0.5, sigs);
            AnlySig(r) = AnlySig(r).LFOTrack(sigs - 1, 1);
            
%             lfo =  AnlySig(r).Measured_LFOs{end, 'Measured_LFO'}{1,1};
%             tAx =  AnlySig(r).Measured_LFOs{end, 'LFO_time_axis'}{1,1};
%             plot(tAx, lfo)
        end
        % Before saving, delete the spectrograms (keep table so each 
        % spec_i can be cross referenced to its signal/proc sig i
        AnlySig(r).Spectrograms = AnlySig(r).Spectrograms(:,3:end);
        
    end
    disp('Saving Data')
    save('Measurements.mat', 'Signals', 'ProcSigs', 'AnlySig', 'LFOFits') 
end

% If LFO fitting hasn't been carried out yet
if size(LFOFits(1).LFOs,1) < 1
    
    % Iterate over rates, load estimated LFO freq
    for r = 1:length(Rate)
        init_f = AnlySig(r).Initial_f;
        
        % For each measured LFO, load the lfo, time axis, length, strt time
        for lfo_i = 1:size(AnlySig(r).Measured_LFOs,1)
            lfo =  AnlySig(r).Measured_LFOs{lfo_i, 'Measured_LFO'}{1,1};
            tAx =  AnlySig(r).Measured_LFOs{lfo_i, 'LFO_time_axis'}{1,1};
            tot = tAx(end) - tAx(1);
            tst = tAx(1);
            
            for t = 1:size(tfrac)
                % Ignore parts of the measured LFO according to tfrac
                tf1 = tfrac(t,1);
                tf2 = tfrac(t,2);
                lfot = lfo(tAx<tot*tf1 + tst | tAx>=tot*tf2 + tst);
                tAxt = tAx(tAx<tot*tf1 + tst | tAx>=tot*tf2 + tst);
                
                % Fit an LFO to the measured LFO
                LFOFits(r) = LFOFits(r).RectSineGridSearch...
                    (lfot, tAxt, init_f, lfo_i, tf1 + 1 - tf2);
            end
        end 
    end
    disp('Saving Data')
    save('Measurements.mat', 'Signals', 'ProcSigs', 'AnlySig', 'LFOFits') 
end