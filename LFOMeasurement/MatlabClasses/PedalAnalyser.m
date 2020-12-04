classdef PedalAnalyser
    properties
        PedalName
        fs
        Settings
        
        ProcessLoc
        SaveDir
        SaveLocs
        LoadDir
        LoadLocs
        
        PmSigs
        PmProcSigs
        PmAnly
        Signals
        ProcSignals
        AnlySig
        Aud_Boundaries
        Chirp_Boundaries
        Data_Boundaries
        LFOs
        
        init_f
        f_st
        f_en
    end
    methods
        function obj = PedalAnalyser(PedalName, Settings, fs, ...
                ProcessLoc, SaveDir, SaveLocs, LoadDir, LoadLocs)
            obj.PedalName = PedalName;
            obj.fs = fs;
            obj.Settings = Settings;
            obj.ProcessLoc = ProcessLoc;
            obj.SaveDir = SaveDir;
            obj.SaveLocs = SaveLocs;
            obj.LoadDir = LoadDir;
            obj.LoadLocs = LoadLocs;
        end
        function obj = PrelimAnalysisGen(obj, c_ty, ch_l, ch_s, f_st, f_en, T)
            Sigs = SignalHolder();
            Sigs = Sigs.SigGen(c_ty, ch_l, ch_s, f_st, f_en, T, obj.fs);
            
            PSigs = SignalProcessor(obj.PedalName, 0, obj.Settings);
            PSigs.SaveLoc = strcat(obj.ProcessLoc, '/', obj.PedalName,...
                '_',[obj.Settings{:}], '-pre');
            PSigs = PSigs.SigProc(Sigs.SigGet(1),obj.Settings, 40, obj.fs);
            
            obj.PmSigs = Sigs;
            obj.PmProcSigs = PSigs;
        end
        function obj = PrelimSigAnalayse(obj)
            sig = SignalProcessor.SigLoad(obj.PmProcSigs.SaveLoc, 1);
            obj.PmProcSigs.ProcessedSignals{1,1} = {sig};
            
            obj.PmAnly = SignalAnalyser...
                    (obj.PmSigs.Signals, obj.PmProcSigs.ProcessedSignals);
            obj.PmAnly = obj.PmAnly.SpecExtract(0.1, 1);
            obj.PmAnly = obj.PmAnly.PrelimAnly(1);

            obj.init_f = obj.PmAnly.Initial_f;
            obj.f_st = obj.PmAnly.Min_f;
            obj.f_en = obj.PmAnly.Max_f;
        end
        function obj = ConstructInputData(obj, method, ch_len, ch_spc, T, seg_len)
            
            Boundaries = [0];
            data_full = [];
            for subDir = 1:length(obj.LoadLocs)
                % Load each dataset
                data = audioread(strcat(obj.LoadDir, '/', obj.LoadLocs(subDir)));
                % Round to nearest second and concatenate
                data = data(1:floor(length(data)/obj.fs)*obj.fs);
                data_full = [data_full; data];
                Boundaries = [Boundaries, length(data_full)];
            end
            
            % Create empty container for the dataset with chirps
            data_chirp_full = [];
            chirp_boundaries = [];
            audio_boundaries = [];
            % Create holder for the chirp trains
            ChirpTrains = SignalHolder();
            F_ST = 0.5*obj.f_st;
            F_EN = 1.25*obj.f_en;
            FS = obj.fs;

            % Split dataset into segments of seg_len length
            for n = 1:floor((length(data_full)/obj.fs)/seg_len)
                % Take chunk of audio
                aud_ch = data_full((n-1)*FS*seg_len + 1:n*FS*seg_len);

                % Generate a chirp train and get the signal
                ChirpTrains = ChirpTrains.SigGen(method, ch_len, ch_spc, F_ST, F_EN, T, FS);
                chirps = ChirpTrains.SigGet(n);

                % Adjust the chirp start times
                ch_sts = ChirpTrains.Signals{n, 'chirp_starts'}{1,1};
                ChirpTrains.Signals{n, 'chirp_starts'} = {ch_sts + (n-1)*(seg_len + T)};
                
                chirp_boundaries(end+1,1:2) = [length(data_chirp_full) + 1,length(data_chirp_full) + T*FS];
                audio_boundaries(end+1,1:2) = [length(data_chirp_full) + T*FS + 1, length(data_chirp_full) + T*FS + length(aud_ch)];
                % Append the audio chunk and chirp train to the audio
                data_chirp_full = [data_chirp_full; chirps(1:end-1); aud_ch]; 
            end
                
            % Generate a chirp train and get the signal
            ChirpTrains = ChirpTrains.SigGen(method, ch_len, ch_spc, F_ST, F_EN, T, FS);
            ch_sts = ChirpTrains.Signals{n+1, 'chirp_starts'}{1,1};
            ChirpTrains.Signals{n+1, 'chirp_starts'} = {ch_sts + (n)*(seg_len + T)};

            aud_ch = data_full(n*FS*seg_len + 1:end);
            ChirpTrains = ChirpTrains.SigGen(method, ch_len, ch_spc, F_ST, F_EN, T, FS);
            ch_sts = ChirpTrains.Signals{n+2, 'chirp_starts'}{1,1};
            ChirpTrains.Signals{n+2, 'chirp_starts'} = {ch_sts + (n)*(seg_len + T) + T + (length(aud_ch)/FS)};

            chirps = ChirpTrains.SigGet(n+1);
            chirpsfinal = ChirpTrains.SigGet(n+2);

            chirp_boundaries(end+1,1:2) = [length(data_chirp_full) + 1,length(data_chirp_full) + T*FS];
            audio_boundaries(end+1,1:2) = [length(data_chirp_full) + T*FS + 1, length(data_chirp_full) + T*FS + length(aud_ch)];
            data_chirp_full = [data_chirp_full; chirps(1:end-1); aud_ch]; 
            
            chirp_boundaries(end+1,1:2) = [length(data_chirp_full) + 1,length(data_chirp_full) + T*FS];
            data_chirp_full = [data_chirp_full; chirpsfinal(1:end-1)];

            audiowrite(strcat(obj.ProcessLoc,'/',obj.PedalName ,[obj.Settings{:}], '-input.wav'), data_chirp_full, FS);
            disp('process dataset with target device and save as ...-output.wav')   
            obj.Aud_Boundaries = audio_boundaries;
            obj.Chirp_Boundaries = chirp_boundaries;
            obj.Data_Boundaries = Boundaries;
            obj.Signals = ChirpTrains;
        end
        function obj = ReadTargetData(obj)
            ProcSigs = SignalProcessor(obj.PedalName, 0, 0);
%             ProcSigs.Settings = Settings;

            ProcSignal = SignalProcessor.SigLoad(strcat(obj.ProcessLoc,'/',obj.PedalName ,[obj.Settings{:}], '-output.wav'), 1);

            AnlySig = SignalAnalyser(obj.Signals.Signals, ProcSigs.ProcessedSignals);
            AnlySig = AnlySig.BatchSpecExtract(0.5, ProcSignal);

            AnlySig.Min_f = obj.f_st;
            AnlySig.Max_f = obj.f_en;
            AnlySig.Initial_f = obj.init_f;

            % Retrieve Chirp Trains from processed audio
            for n = 1:size(obj.Signals.Signals,1)
                AnlySig = AnlySig.LFOTrack(n, 5, 2);
            end
            
            obj.AnlySig = AnlySig;
        end
        function obj = LFOFitting(obj, inv)
            
%             full_lfo = [];
%             full_tAx = [];
%             for n = 1:size(obj.AnlySig.Measured_LFOs,1)
%                 full_lfo = [full_lfo; obj.AnlySig.Measured_LFOs{n,1}{1,1}];
%                 full_tAx = [full_tAx; obj.AnlySig.Measured_LFOs{n,2}{1,1}];
%             end
            % Invert the signal if inv is chosen
            
            param_mat = [];
            for n = 1:size(obj.AnlySig.Measured_LFOs,1) - 1
                
                full_lfo = [obj.AnlySig.Measured_LFOs{n,1}{1,1}; obj.AnlySig.Measured_LFOs{n+1,1}{1,1}];
                full_tAx = [obj.AnlySig.Measured_LFOs{n,2}{1,1}; obj.AnlySig.Measured_LFOs{n+1,2}{1,1}];
                
                if inv == 1
                    pre_lfo = full_lfo;
                    pre_C = min(full_lfo);
                    pre_Amp = max(full_lfo) - pre_C;
                    full_lfo = -full_lfo + max(full_lfo) + min(full_lfo);
                end

                C = min(full_lfo);
                Amp = max(full_lfo) - C;
            

           
                tst1 = obj.AnlySig.Measured_LFOs{n,2}{1,1}(1);
                ten1 = obj.AnlySig.Measured_LFOs{n,2}{1,1}(end);
                tst2 = obj.AnlySig.Measured_LFOs{n+1,2}{1,1}(1);
                ten2 = obj.AnlySig.Measured_LFOs{n+1,2}{1,1}(end);

                lfo1 = full_lfo(full_tAx >= tst1 & full_tAx <= ten1);
                lfo2 = full_lfo(full_tAx >= tst2 & full_tAx <= ten2);
                tAx1 = full_tAx(full_tAx >= tst1 & full_tAx <= ten1);
                tAx2 = full_tAx(full_tAx >= tst2 & full_tAx <= ten2);

                [lfoP, params] = LFOFitter.SimpleSineFit([lfo1;lfo2], [tAx1;tAx2], obj.init_f, Amp, C);
                lfoP1 = lfoP(1:length(tAx1));
                lfoP2 = lfoP(length(tAx2) + 1:end);
                error1 = mean((lfoP1 - lfo1).^2);
                error2 = mean((lfoP2 - lfo2).^2);
                maxe1 = max(abs(lfoP1 - lfo1));
                maxe2 = max(abs(lfoP2 - lfo2));
%                 error1 - error2
%                 maxe1 - maxe2

%                 plot([tAx1;tAx2], [lfo1;lfo2]);
%                 hold on
%                 plot([tAx1;tAx2], lfoP);
%                 hold off
                if inv == 1
                   Amp = pre_Amp;
                   C = pre_C + pre_Amp;
                end

                param_mat(end+1, :) = [params, Amp, C , error1 - error2, abs(maxe1 - maxe2)];

%                 param_mat(end+1, :) = params;
            end
            

            obj.LFOs = {param_mat};
        end
        function obj = DataMake(obj, inv)
            
            full_lfo = [];
            full_tAx = [];
            fs = obj.fs;
            
            Aud_Times = (obj.Aud_Boundaries - 1)/fs;
            Chirp_Times = (obj.Chirp_Boundaries - 1)/fs;
            
            if inv == 1
                si = -1;
            else
                si = 1;
            end

            
            for n = 1:size(obj.AnlySig.Measured_LFOs,1) - 1
                full_lfo = [full_lfo; obj.AnlySig.Measured_LFOs{n,1}{1,1}];
                full_tAx = [full_tAx; obj.AnlySig.Measured_LFOs{n,2}{1,1}];
                
                f = obj.LFOs{1,1}(n,1);
                phi = obj.LFOs{1,1}(n,2);
                Amp = obj.LFOs{1,1}(n,3);
                C = obj.LFOs{1,1}(n,4);
                aud_tAx = Aud_Times(n,1):1/obj.fs:Aud_Times(n,2);
                aud_lfo = si*Amp*abs(sin(f*pi*aud_tAx + phi)) + C;
                
                full_lfo = [full_lfo; aud_lfo'];
                full_tAx = [full_tAx; aud_tAx'];
                
                t1 = obj.AnlySig.Measured_LFOs{n,2}{1,1};
                t2 = obj.AnlySig.Measured_LFOs{n+1,2}{1,1};
%                 plot(t1, obj.AnlySig.Measured_LFOs{n,1}{1,1})
%                 hold on
%                 plot(t2, obj.AnlySig.Measured_LFOs{n+1,1}{1,1})
% %                 plot(aud_tAx, aud_lfo)
%                 plot([t1;t2], si*Amp*abs(sin(f*pi*[t1;t2] + phi)) + C)
%                 hold off
            end

            full_lfo = full_lfo - min(full_lfo);
            full_lfo = full_lfo/max(full_lfo);
            
            Input = [];
            Target = [];
            
            ProcSignal = SignalProcessor.SigLoad(obj.SaveLoc, 2);
            Input_Raw = ProcSignal(:, 2);
            Target_Raw = ProcSignal(:, 1);
            
            for n = 1:length(Aud_Times)
                st_t = Aud_Times(n,1);
                en_t = Aud_Times(n,2);
                st_s = obj.Aud_Boundaries(n,1);
                en_s = obj.Aud_Boundaries(n,2);
                
                In = Input_Raw(st_s:en_s);
                Lfo = full_lfo(full_tAx >= st_t & full_tAx <= en_t);
                Targ = Target_Raw(st_s:en_s);
                
                Input = [Input; [In, Lfo]];
                Target = [Target; Targ];
            end
            
            Input(:,1) = Input(:,1)/max(abs(Input(:,1)));
            Input(:,1) = Input(:,1)*0.95;
            Target = Target/max(abs(Target));
            Target = Target*0.95;
            
            for subDir = 1:length(obj.DataLocs)
                
                st = obj.Data_Boundaries(subDir) + 1;
                en = obj.Data_Boundaries(subDir + 1);
                
                loc = strcat(obj.DataDir, obj.DataLocs(subDir), '/', obj.PedalName);
                audiowrite(strcat(loc, '-input.wav'), Input(st:en, :), fs);
                audiowrite(strcat(loc, '-target.wav'), Target(st:en, :), fs);
            end
            
            

        end
        function obj = DataMakeCutFr(obj, inv, N)
            
            full_lfo = [];
            full_tAx = [];
            fs = obj.fs;
            
            Aud_Times = (obj.Aud_Boundaries - 1)/fs;
            Chirp_Times = (obj.Chirp_Boundaries - 1)/fs;
            
            if inv == 1
                si = -1;
            else
                si = 1;
            end
            


            frqs = obj.LFOs{1, 1}(:,1);
            [bin, bound] = discretize(frqs, N);
            win_bin = mode(bin);
            
%             things= zeros(20, 1);
%             for n = 1:20
%                 [bin, bound] = discretize(frqs, n);
%                 win_bin = mode(bin);
%                 things(n) = sum(bin == win_bin);
%             end
            
            for n = 1:size(obj.AnlySig.Measured_LFOs,1) - 1

                full_lfo = [full_lfo; obj.AnlySig.Measured_LFOs{n,1}{1,1}];
                full_tAx = [full_tAx; obj.AnlySig.Measured_LFOs{n,2}{1,1}];

                f = obj.LFOs{1,1}(n,1);
                phi = obj.LFOs{1,1}(n,2);
                Amp = obj.LFOs{1,1}(n,3);
                C = obj.LFOs{1,1}(n,4);

                aud_tAx = Aud_Times(n,1):1/obj.fs:Aud_Times(n,2);
                aud_lfo = si*Amp*abs(sin(f*pi*aud_tAx + phi)) + C;


                full_lfo = [full_lfo; aud_lfo'];
                full_tAx = [full_tAx; aud_tAx'];

                
%                 t1 = obj.AnlySig.Measured_LFOs{n,2}{1,1};
%                 t2 = obj.AnlySig.Measured_LFOs{n+1,2}{1,1};
%                 plot(t1, obj.AnlySig.Measured_LFOs{n,1}{1,1})
%                 hold on
%                 plot(t2, obj.AnlySig.Measured_LFOs{n+1,1}{1,1})
% %                 plot(aud_tAx, aud_lfo)
%                 plot([t1;t2], si*Amp*abs(sin(f*pi*[t1;t2] + phi)) + C)
%                 hold off
            end

            full_lfo = full_lfo - min(full_lfo);
            full_lfo = full_lfo/max(full_lfo);
            
            Input = [];
            Target = [];
            
            ProcSignal = SignalProcessor.SigLoad(obj.SaveLoc, 2);
            Input_Raw = ProcSignal(:, 2);
            Target_Raw = ProcSignal(:, 1);
            
            for n = 1:length(Aud_Times)
                st_t = Aud_Times(n,1);
                en_t = Aud_Times(n,2);
                st_s = obj.Aud_Boundaries(n,1);
                en_s = obj.Aud_Boundaries(n,2);
                
                In = Input_Raw(st_s:en_s);
                Lfo = full_lfo(full_tAx >= st_t & full_tAx <= en_t);
                Targ = Target_Raw(st_s:en_s);
                
                if bin(n) == win_bin
                    Input = [Input; [In, Lfo]];
                    Target = [Target; Targ];
                end
                

            end
            
            Input(:,1) = Input(:,1)/max(abs(Input(:,1)));
            Input(:,1) = Input(:,1)*0.95;
            Target = Target/max(abs(Target));
            Target = Target*0.95;
            
            t_l = length(Target);
            t_l_s = t_l/obj.fs;
            boundies = [0, 0.7, 0.85, 1];
            
            for subDir = 1:length(obj.DataLocs)
                
%                 st = obj.Data_Boundaries(subDir) + 1;
%                 en = obj.Data_Boundaries(subDir + 1);
                st_sec = floor(boundies(subDir)*t_l_s);
                st = (st_sec*obj.fs) + 1;
                en_sec = floor(boundies(subDir + 1)*t_l_s);
                en = en_sec*obj.fs;
                 
                loc = strcat(obj.DataDir, obj.DataLocs(subDir), '/', obj.PedalName);
                audiowrite(strcat(loc, 'frcut', num2str(N), '-input.wav'), Input(st:en, :), fs);
                audiowrite(strcat(loc, 'frcut', num2str(N), '-target.wav'), Target(st:en, :), fs);
            end
            winners = find(bin == win_bin);
            writematrix(winners,strcat('frcut', num2str(N), '.csv')) 
          
        end
        function obj = DataMakeCutMxEr(obj, inv, N)
            
            full_lfo = [];
            full_tAx = [];
            fs = obj.fs;
            
            Aud_Times = (obj.Aud_Boundaries - 1)/fs;
            Chirp_Times = (obj.Chirp_Boundaries - 1)/fs;
            
            if inv == 1
                si = -1;
            else
                si = 1;
            end
            
            MxErr = obj.LFOs{1, 1}(:,6);
            
            win_bin = MxErr < N;
            
            for n = 1:size(obj.AnlySig.Measured_LFOs,1) - 1

                full_lfo = [full_lfo; obj.AnlySig.Measured_LFOs{n,1}{1,1}];
                full_tAx = [full_tAx; obj.AnlySig.Measured_LFOs{n,2}{1,1}];

                f = obj.LFOs{1,1}(n,1);
                phi = obj.LFOs{1,1}(n,2);
                Amp = obj.LFOs{1,1}(n,3);
                C = obj.LFOs{1,1}(n,4);

                aud_tAx = Aud_Times(n,1):1/obj.fs:Aud_Times(n,2);
                aud_lfo = si*Amp*abs(sin(f*pi*aud_tAx + phi)) + C;


                full_lfo = [full_lfo; aud_lfo'];
                full_tAx = [full_tAx; aud_tAx'];

                
%                 t1 = obj.AnlySig.Measured_LFOs{n,2}{1,1};
%                 t2 = obj.AnlySig.Measured_LFOs{n+1,2}{1,1};
%                 plot(t1, obj.AnlySig.Measured_LFOs{n,1}{1,1})
%                 hold on
%                 plot(t2, obj.AnlySig.Measured_LFOs{n+1,1}{1,1})
% %                 plot(aud_tAx, aud_lfo)
%                 plot([t1,t2], si*Amp*abs(sin(f*pi*[t1,t2] + phi)) + C)
%                 hold off
            end

            full_lfo = full_lfo - min(full_lfo);
            full_lfo = full_lfo/max(full_lfo);
            
            Input = [];
            Target = [];
            
            ProcSignal = SignalProcessor.SigLoad(obj.SaveLoc, 2);
            Input_Raw = ProcSignal(:, 2);
            Target_Raw = ProcSignal(:, 1);
            
            for n = 1:length(Aud_Times)
                st_t = Aud_Times(n,1);
                en_t = Aud_Times(n,2);
                st_s = obj.Aud_Boundaries(n,1);
                en_s = obj.Aud_Boundaries(n,2);
                
                In = Input_Raw(st_s:en_s);
                Lfo = full_lfo(full_tAx >= st_t & full_tAx <= en_t);
                Targ = Target_Raw(st_s:en_s);
                
                if win_bin(n)
                    Input = [Input; [In, Lfo]];
                    Target = [Target; Targ];
                end
                

            end
            
            Input(:,1) = Input(:,1)/max(abs(Input(:,1)));
            Input(:,1) = Input(:,1)*0.95;
            Target = Target/max(abs(Target));
            Target = Target*0.95;
            
            t_l = length(Target);
            t_l_s = t_l/obj.fs;
            boundies = [0, 0.7, 0.85, 1];
            
            for subDir = 1:length(obj.DataLocs)
                
%                 st = obj.Data_Boundaries(subDir) + 1;
%                 en = obj.Data_Boundaries(subDir + 1);
                st_sec = floor(boundies(subDir)*t_l_s);
                st = (st_sec*obj.fs) + 1;
                en_sec = floor(boundies(subDir + 1)*t_l_s);
                en = en_sec*obj.fs;
                 
                loc = strcat(obj.DataDir, obj.DataLocs(subDir), '/', obj.PedalName);
                audiowrite(strcat(loc, 'MxErcut', num2str(N), '-input.wav'), Input(st:en, :), fs);
                audiowrite(strcat(loc, 'MxErcut', num2str(N), '-target.wav'), Target(st:en, :), fs);
            end
            
            winners = find(MxErr < N);
            writematrix(winners,strcat('MxErcut', num2str(N), '.csv')) 

        end
        function obj = DataMakeCutAvEr(obj, inv, N)
            
            full_lfo = [];
            full_tAx = [];
            fs = obj.fs;
            
            Aud_Times = (obj.Aud_Boundaries - 1)/fs;
            Chirp_Times = (obj.Chirp_Boundaries - 1)/fs;
            
            if inv == 1
                si = -1;
            else
                si = 1;
            end
            

            
            AvErr = abs(obj.LFOs{1, 1}(:,5));
            
            win_bin = AvErr < N;
            
            for n = 1:size(obj.AnlySig.Measured_LFOs,1) - 1

                full_lfo = [full_lfo; obj.AnlySig.Measured_LFOs{n,1}{1,1}];
                full_tAx = [full_tAx; obj.AnlySig.Measured_LFOs{n,2}{1,1}];

                f = obj.LFOs{1,1}(n,1);
                phi = obj.LFOs{1,1}(n,2);
                Amp = obj.LFOs{1,1}(n,3);
                C = obj.LFOs{1,1}(n,4);

                aud_tAx = Aud_Times(n,1):1/obj.fs:Aud_Times(n,2);
                aud_lfo = si*Amp*abs(sin(f*pi*aud_tAx + phi)) + C;


                full_lfo = [full_lfo; aud_lfo'];
                full_tAx = [full_tAx; aud_tAx'];

                
%                 t1 = obj.AnlySig.Measured_LFOs{n,2}{1,1};
%                 t2 = obj.AnlySig.Measured_LFOs{n+1,2}{1,1};
%                 plot(t1, obj.AnlySig.Measured_LFOs{n,1}{1,1})
%                 hold on
%                 plot(t2, obj.AnlySig.Measured_LFOs{n+1,1}{1,1})
% %                 plot(aud_tAx, aud_lfo)
%                 plot([t1,t2], si*Amp*abs(sin(f*pi*[t1,t2] + phi)) + C)
%                 hold off
            end

            full_lfo = full_lfo - min(full_lfo);
            full_lfo = full_lfo/max(full_lfo);
            
            Input = [];
            Target = [];
            
            ProcSignal = SignalProcessor.SigLoad(obj.SaveLoc, 2);
            Input_Raw = ProcSignal(:, 2);
            Target_Raw = ProcSignal(:, 1);
            
            for n = 1:length(Aud_Times)
                st_t = Aud_Times(n,1);
                en_t = Aud_Times(n,2);
                st_s = obj.Aud_Boundaries(n,1);
                en_s = obj.Aud_Boundaries(n,2);
                
                In = Input_Raw(st_s:en_s);
                Lfo = full_lfo(full_tAx >= st_t & full_tAx <= en_t);
                Targ = Target_Raw(st_s:en_s);
                
                if win_bin(n)
                    Input = [Input; [In, Lfo]];
                    Target = [Target; Targ];
                end
                

            end
            
%             Input(:,1) = Input(:,1)/max(abs(Input(:,1)));
%             Input(:,1) = Input(:,1)*0.95;
%             Target = Target/max(abs(Target));
%             Target = Target*0.95;
            
            t_l = length(Target);
            t_l_s = t_l/obj.fs;
            boundies = [0, 0.7, 0.85, 1];
            
            for subDir = 1:length(obj.DataLocs)
                
%                 st = obj.Data_Boundaries(subDir) + 1;
%                 en = obj.Data_Boundaries(subDir + 1);
                st_sec = floor(boundies(subDir)*t_l_s);
                st = (st_sec*obj.fs) + 1;
                en_sec = floor(boundies(subDir + 1)*t_l_s);
                en = en_sec*obj.fs;
                 
                loc = strcat(obj.DataDir, obj.DataLocs(subDir), '/', obj.PedalName);
                audiowrite(strcat(loc, 'AvErcut', num2str(N), '-input.wav'), Input(st:en, :), fs);
                audiowrite(strcat(loc, 'AvErcut', num2str(N), '-target.wav'), Target(st:en, :), fs);
            end
            
            winners = find(AvErr < N);
            writematrix(winners,strcat('AvErcut', num2str(N), '.csv')) 

        end
        function obj = DataMakeSingles(obj, inv)
            
            full_lfo = [];
            full_tAx = [];
            fs = obj.fs;
            
            Aud_Times = (obj.Aud_Boundaries - 1)/fs;
            Chirp_Times = (obj.Chirp_Boundaries - 1)/fs;
            
            if inv == 1
                si = -1;
            else
                si = 1;
            end
            
            for n = 1:size(obj.AnlySig.Measured_LFOs,1) - 1

                full_lfo = [full_lfo; obj.AnlySig.Measured_LFOs{n,1}{1,1}];
                full_tAx = [full_tAx; obj.AnlySig.Measured_LFOs{n,2}{1,1}];

                f = obj.LFOs{1,1}(n,1);
                phi = obj.LFOs{1,1}(n,2);
                Amp = obj.LFOs{1,1}(n,3);
                C = obj.LFOs{1,1}(n,4);

                aud_tAx = Aud_Times(n,1):1/obj.fs:Aud_Times(n,2);
                aud_lfo = si*Amp*abs(sin(f*pi*aud_tAx + phi)) + C;

                full_lfo = [full_lfo; aud_lfo'];
                full_tAx = [full_tAx; aud_tAx'];

                
%                 t1 = obj.AnlySig.Measured_LFOs{n,2}{1,1};
%                 t2 = obj.AnlySig.Measured_LFOs{n+1,2}{1,1};
%                 plot(t1, obj.AnlySig.Measured_LFOs{n,1}{1,1})
%                 hold on
%                 plot(t2, obj.AnlySig.Measured_LFOs{n+1,1}{1,1})
% %                 plot(aud_tAx, aud_lfo)
%                 plot([t1,t2], si*Amp*abs(sin(f*pi*[t1,t2] + phi)) + C)
%                 hold off
            end

            full_lfo = full_lfo - min(full_lfo);
            full_lfo = full_lfo/max(full_lfo);
            
            
            ProcSignal = SignalProcessor.SigLoad(obj.SaveLoc, 2);
            Input_Raw = ProcSignal(:, 2);
            Target_Raw = ProcSignal(:, 1);
            
            for n = 1:length(Aud_Times)
                st_t = Aud_Times(n,1);
                en_t = Aud_Times(n,2);
                st_s = obj.Aud_Boundaries(n,1);
                en_s = obj.Aud_Boundaries(n,2);
                
                In = Input_Raw(st_s:en_s);
                Lfo = full_lfo(full_tAx >= st_t & full_tAx <= en_t);
                Targ = Target_Raw(st_s:en_s);
                
                Input = [In, Lfo];
                Target = Targ;
            
                t_l = length(Target);
                t_l_s = t_l/obj.fs;
                boundies = [0, 0.85, 0.975, 1];

                for subDir = 1:length(obj.DataLocs)

    %                 st = obj.Data_Boundaries(subDir) + 1;
    %                 en = obj.Data_Boundaries(subDir + 1);
                    st_sec = floor(boundies(subDir)*t_l_s);
                    st = (st_sec*obj.fs) + 1;
                    en_sec = floor(boundies(subDir + 1)*t_l_s);
                    en = en_sec*obj.fs;

                    loc = strcat(obj.DataDir, obj.DataLocs(subDir), '/', obj.PedalName);
                    audiowrite(strcat(loc, 'Singles', num2str(n), '-input.wav'), Input(st:en, :), fs);
                    audiowrite(strcat(loc, 'Singles', num2str(n), '-target.wav'), Target(st:en, :), fs);
                end
            end


        end
    end
end