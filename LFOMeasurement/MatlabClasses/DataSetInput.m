function [ChirpTrains, Boundaries] = DataSetInput...
    (DataDir, DataLocs, SaveLoc, ...
                        method, ch_len, ch_spc, f_st, f_en, T, seg_len, fs)

    Boundaries = [0];
    data_full = [];
    for subDir = 1:length(DataLocs)
        % Load each dataset
        data = audioread(strcat(DataDir, '/', DataLocs(subDir), '/input.wav'));
        % Round to nearest second and concatenate
        data = data(1:floor(length(data)/fs)*fs);
        data_full = [data_full; data];
        Boundaries = [Boundaries, length(data_full)];
    end

    % Create empty container for the dataset with chirps
    data_chirp_full = [];
    % Create holder for the chirp trains
    ChirpTrains = SignalHolder();
    
    % Split dataset into segments of seg_len length
    for n = 1:floor((length(data_full)/44100)/seg_len)
        % Take chunk of audio
        aud_ch = data_full((n-1)*fs*seg_len + 1:n*fs*seg_len);
        
        % Generate a chirp train and get the signal
        ChirpTrains = ChirpTrains.SigGen(method, ch_len, ch_spc, f_st, f_en, T, fs);
        chirps = ChirpTrains.SigGet(n);
        
        % Adjust the chirp start times
        ch_sts = ChirpTrains.Signals{n, 'chirp_starts'}{1,1};
        ChirpTrains.Signals{n, 'chirp_starts'} = {ch_sts + (n-1)*(seg_len + T)};
        
        % Append the audio chunk and chirp train to the audio
        data_chirp_full = [data_chirp_full; chirps(1:end-1); aud_ch]; 
    end
    
    % Generate a chirp train and get the signal
    ChirpTrains = ChirpTrains.SigGen(method, ch_len, ch_spc, f_st, f_en, T, fs);
    ch_sts = ChirpTrains.Signals{n+1, 'chirp_starts'}{1,1};
    ChirpTrains.Signals{n+1, 'chirp_starts'} = {ch_sts + (n)*(seg_len + T)};
    
    aud_ch = data_full(n*fs*seg_len + 1:end);
    ChirpTrains = ChirpTrains.SigGen(method, ch_len, ch_spc, f_st, f_en, T, fs);
    ch_sts = ChirpTrains.Signals{n+2, 'chirp_starts'}{1,1};
    ChirpTrains.Signals{n+2, 'chirp_starts'} = {ch_sts + (n)*(seg_len + T) + T + (length(aud_ch)/fs)};
    
    chirps = ChirpTrains.SigGet(n+1);
    chirpsfinal = ChirpTrains.SigGet(n+2);
    
    data_chirp_full = [data_chirp_full; chirps(1:end-1); aud_ch; chirpsfinal(1:end-1)]; 
    
    ChirpTrains.SaveLoc = SaveLoc;
    ChirpTrains.SegLen = seg_len;
    audiowrite(strcat(SaveLoc, '-input.wav'), data_chirp_full, fs);
    disp('process dataset with target device and save as ...-output.wav')   
end