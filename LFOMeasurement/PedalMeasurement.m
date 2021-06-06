addpath('MatlabClasses')
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
% By Alec Wright
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
% This script creates training data for LFO-modulated time-varying effects
% such as phaser and flanger pedals, which can then be used to create a
% neural network model of the effects pedal. The script has two stages

% First, the script takes a mono audio file, and periodically inserts an 
% LFO measurement signal into it, then saves the audio file and instructs 
% the user to play it through the effect and record the output, then place
% the resulting file in this directory. 

% Second, the script reads this file and analyses the LFO measurement
% signal, the user is prompted to click on some points on the spectrogram
% image to provide initial estimates for the algorithim. Then the script
% will estimate the LFO behaviour over the duration of the recording, and
% save the LFO-labelled data for use in training a neural network model
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 

% Stage = 1 to create input data for pedal
% Stage = 2 to read audio processed by pedal
% Stage = 3 to save the LFO-annotated data for neural net training
if ~exist('stage')
    stage = 1;
end

%% 0 - Dataset Info
% input data path, save path, length of data between measurement signals, 
% measurment signal lengths
InData = 'RawData/Input-Dataset.wav';
SaveAudioLoc = 'PedalData';
SaveDataLoc = 'MeasurementOutputs';
ChunkLen = 80;
ChirpSigLen = 20;
TargetPedal = 'BehPhaser';

% flag to determine whether to save multiple shorter datasets or 1 long one
splitOutputs = 1;
% frange, LFO will be normalised such that -1 and 1 corresponds to frange
frange = [50,1000];

%% 1 - Stage 1 - loads a data file, adds measurement signal then saves it,
% for the user to process with the target pedal
switch stage
    case 1
        % Read training data, normalise, calculate number of chunks
        [InputData, fs] = audioread(InData);
        InputData = 0.95*InputData./max(abs(InputData));
        AudioLen = length(InputData)/fs;
        TotalChunks = AudioLen/ChunkLen;
        
        % Create signal holder object, set save location
        Signals = SignalHolder();
        
        TestSig = [];
        %Generate chirp signals, insert into input data
        for n = 1:TotalChunks
            Signals = Signals.SigGen(1, 19, 20, [20, 22050], ChirpSigLen, fs, ChunkLen);
            Chirp_Train = Signals.SigGet(n);
            TestSig = [TestSig; Chirp_Train(1:end-1);...
                            InputData(fs*(n-1)*ChunkLen+1:fs*n*ChunkLen)];
        end
        Signals = Signals.SigGen(1, 19, 20, [20, 22050], ChirpSigLen, fs, 0);
        Chirp_Train = Signals.SigGet(n);
        TestSig = [TestSig; Chirp_Train(1:end-1)];
        
        save(strcat(SaveDataLoc, '/Signals'), 'Signals');
        audiowrite(strcat(SaveAudioLoc, '/', TargetPedal, '-Input.wav'), TestSig, fs)
        
        disp('Process audio with pedal and save as "Pedal-Output.wav"')
        stage = 2
        
%% 2 - Stage 2 - Loads the processed data file, asks user to initiate LFO 
% tracking
    case 2
        if ~exist('Signals','var')
            load(strcat(SaveDataLoc, '/Signals'));
        end
        
        ProcSigs = SignalProcessor(TargetPedal, '3');
        ProcSigs = ...
            ProcSigs.LoadFromFile(strcat(SaveAudioLoc,'/',TargetPedal), Signals);
        
        AnlySig = SignalAnalyser(Signals.Signals, ProcSigs.ProcessedSignals);
        % Extract the chirp spectra, with freq range [0 3000]
        for n = 1:size(Signals.Signals,1)
            AnlySig = AnlySig.SpecExtract(0.5, n, [0 3000]);
            if n == 1
                AnlySig = AnlySig.PrelimAnly(1);
            end
            AnlySig = AnlySig.LFOTrack(n, 1, 1, 2);
        end
        save(strcat(SaveDataLoc, '/AnlySig'), 'AnlySig', 'ProcSigs', 'Signals');
        
        plot(AnlySig.Measured_LFOs.LFO_time_axis{1,1},...
            AnlySig.Measured_LFOs.Measured_LFO{1,1})
        ylabel('Frequency (Hz)')
        xlabel('Time (s)')
        
        stage = 3
        
%% 2 - Stage 3 
    case 3
        if ~exist('AnlySig','var')
            load(strcat(SaveDataLoc, '/AnlySig'));
        end
        AnlySig = AnlySig.ResetFits();
        [InputData, fs] = audioread(InData);
        [TargetData, fs] = audioread(strcat(SaveAudioLoc,'/',TargetPedal, '-Output.wav'));
        InputData = 0.95*InputData./max(abs(InputData));
        
        input = [];
        target = [];
        st_ind = 1;
        LFO_en_t = 0;
        
        for n = 1:size(Signals.Signals,1)-1
            % Measurement signals to fit the LFO to
            fit_sigs = [n, n+1];
            % Fit the LFO
            AnlySig = AnlySig.SineFit(fit_sigs);
            
            % Find index for end of thsi chunk of audio
            en_ind = st_ind + AnlySig.Signals.chunk_len(n)*fs - 1;
            % Take chunk of audio
            data_chan1 = InputData(st_ind:en_ind);
            % Calculate start/end time for LFO
            LFO_st_t = LFO_en_t + AnlySig.Signals.T(n);
            LFO_en_t = LFO_st_t + AnlySig.Signals.chunk_len(n);
            
            LFO = AnlySig.Fitted_LFOs.LFO_Func{1}(LFO_st_t:1/fs:(LFO_en_t-1/fs));
            
            LFO = (LFO'-min(frange))/(max(frange)-min(frange));
            
            if splitOutputs == 1
                fileName = strcat(SaveDataLoc, '/',...
                    TargetPedal, '_set', num2str(n));
                audiowrite(strcat(fileName,'-input.wav'),...
                    [data_chan1, LFO],fs)
                audiowrite(strcat(fileName,'-target.wav'),...
                    ProcSigs.ProcessedAudio.processed_audio{n},fs)
            elseif splitOutput == 0
                output = [output; [data_chan1, LFO]];
            end
            
            st_ind = en_ind+1;
        end
        
        if splitOutputs == 0
            fileName = strcat(SaveDataLoc, '/', TargetPedal, '.wav');
            audiowrite(fileName, output, fs)
        end  
        
end
