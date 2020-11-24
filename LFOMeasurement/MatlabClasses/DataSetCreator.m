function [ProcSigs, AnlySigs, LFOFit] = ...
                    DataSetCreator(Signals, Rate, Settings)
   
    ProcSigs = SignalProcessor(PedalName, 0, Rate);
    ProcSigs.Settings = Settings;
    
    ProcSignal = SignalProcessor.SigLoad(Signals.SaveLoc);
    
    % Retrieve Chirp Trains from processed audio
    for n = 1:size(Signals.Signals,2)
        
        
        
        % Take chunk of audio
        aud_ch = data_full((n-1)*fs*data_ch + 1:n*fs*data_ch);
        
        % Generate a chirp train and get the signal
        ChirpTrains = ChirpTrains.SigGen(method, ch_len, ch_spc, f_st, f_en, T, fs);
        chirps = ChirpTrains.SigGet(n);
        
        % Adjust the chirp start times
        ch_sts = ChirpTrains.Signals{n, 'chirp_starts'}{1,1};
        ChirpTrains.Signals{n, 'chirp_starts'} = {ch_sts + (n-1)*(data_ch + T)};
        
        % Append the audio chunk and chirp train to the audio
        data_chirp_full = [data_chirp_full; chirps(1:end-1); aud_ch]; 
    end
    
end