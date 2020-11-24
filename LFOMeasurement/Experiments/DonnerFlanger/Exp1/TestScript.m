addpath('../../MatlabClasses')
addpath('../../../../Dataset')
import SignalHolder
import SignalProcessor
import SignalAnalyser
import LFOFitter

% Signal Processing Info
PedalName = 'DonnerFlanger';
Digi = 0;
rate = 9;
settings = {'Color', 12, 'Range', 12};
fs = 44100;

DataSetLocation = '../../../../Dataset';

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
        
        Signals = SignalHolder();
        ProcSigs = SignalProcessor(PedalName, Digi, rate);
        ProcSigs.Settings = settings;
        % Create Preliminary Measurement signal
        Signals = Signals.SigGen(1, 19, 20, 100, 1600, 10, fs);
        ProcSigs = ProcSigs.SigProc(Signals.SigGet(1), rate, 40, fs);
        while 1
            try
                ProcSigs = ProcSigs.SigLoad(1);
                break
            catch
                pause
                disp('-out file not found, process file and press any key')
            end
        end
        % Analyse the preliminary signal to estiamte LFO freq and
        % parameters for test signal chirp
        AnlySig = SignalAnalyser(Signals.Signals, ProcSigs.ProcessedSignals);
        AnlySig = AnlySig.SpecExtract(0.5, 1);
        AnlySig = AnlySig.PrelimAnly(1);
        


        Signals = Signals.SigGen(1, 15, 20, st_f, en_f, 10, fs);
        ProcSigs = ProcSigs.SigProc(Signals.SigGet(2), rate, 40, fs);
        while 1
            try
                ProcSigs = ProcSigs.SigLoad(2);
                break
            catch
                pause
                disp('-out file not found, process file and press any key')
            end
        end

        disp('Saving Data')
        save ('Measurements.mat', 'Signals', 'ProcSigs', 'AnlySig') 
    end

end

clear DataSignals

if ~exist('DataSignals')

    trdata = audioread(strcat(DataSetLocation, '/train/input.wav'));
    vadata = audioread(strcat(DataSetLocation, '/val/input.wav'));
    tedata = audioread(strcat(DataSetLocation, '/test/input.wav'));
    trdata = trdata(1:floor(length(trdata)/fs)*fs);
    vadata = vadata(1:floor(length(vadata)/fs)*fs);
    tedata = tedata(1:floor(length(tedata)/fs)*fs);
    
    data_full = [trdata; vadata; tedata];
    data_chirp_full = [];
    
    st_f = 0.8*AnlySig.Min_f;
    en_f = 1.1*AnlySig.Max_f;
    
    DataSignals = SignalHolder();
    
    data_ch = 40;
    T  = 10;
    
    for n = 1:ceil((length(data_full)/44100)/data_ch)
        
        aud_ch = data_full((n-1)*fs*data_ch + 1:n*fs*data_ch);
        
        DataSignals = DataSignals.SigGen(1, 15, 25, st_f, en_f, T, fs);
        chirps = DataSignals.SigGet(n);
        ch_sts = DataSignals.Signals{n, 'chirp_starts'}{1,1};
        DataSignals.Signals{n, 'chirp_starts'} = {ch_sts + (n-1)*(data_ch + T)};
     
        data_chirp_full = [data_chirp_full; chirps(1:end-1); aud_ch]; 
    end
    
end




% 
