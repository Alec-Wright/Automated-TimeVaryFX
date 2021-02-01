classdef SignalAnalyser
    properties
        Spectrograms
        Measured_LFOs
        Signals
        ProcSigs
        Initial_f
        Min_f
        Max_f
        Ntch_Jump_Fac
    end
    methods
        function obj = SignalAnalyser(Signals, ProcSigs)
            
            obj.Signals = Signals;
            obj.ProcSigs = ProcSigs;
            obj.Ntch_Jump_Fac = 2;
            
            table_header = [["spectrogram", "cell"]; ...
                        ["fAx", "cell"]; ...
                        ["f_to_ind", "cell"]; ...
                        ["f_round", "cell"]; ...
                        ["processed_sig", "int16"];...
                        ["sig_num", "int16"]
                        ["notch_start_fr", "double"];...
                        ["notch_start_i", "int16"]];
            
            % Make table using fieldnames & value types from above
            obj.Spectrograms = table('Size',[0,size(table_header,1)],... 
                'VariableNames', table_header(:,1),...
                'VariableTypes', table_header(:,2));
            
            table_header = [["Measured_LFO", "cell"];...
            ["LFO_time_axis", "cell"];...
            ["smooth_factor", "double"]
            ["spec_num", "int16"]];

            % Make table using fieldnames & value types from above
            obj.Measured_LFOs = table('Size',[0,size(table_header,1)],... 
                'VariableNames', table_header(:,1),...
                'VariableTypes', table_header(:,2));
            
        end
        function obj = SplitSpecExtract(obj, f_res, proc_i, split_num)
            
            proc_sig = obj.ProcSigs{proc_i, 'processed_signal'}{1,1};
            sig_i = obj.ProcSigs{proc_i,'signal_number'};
            
            fs = obj.Signals{sig_i,'fs'};
            ch_i = obj.Signals{sig_i,'chirp_starts'}{1,1};

            %Desired frequency resolution
            fAx_full = 0:f_res:fs - f_res;
            % find the frequency axis start/end bins on the full frequency axis
            ax_st = round(obj.Signals{sig_i,'f_st'}/f_res) + 1;
            ax_en = round(obj.Signals{sig_i,'f_en'}/f_res) + 1;
            % Truncate the frequency axis to just cover frequencies in the chirp
            fAx = fAx_full(ax_st:ax_en);
            % Calcultate the required FFT length to achieve f_res freq resolution
            N = fs/f_res;
            
            
          % Find frequency resolution and starting frequency from fAx
            f_st = fAx(1);
            f_to_ind = @(freq) round(1 + (freq - f_st)/f_res);
            f_round = @(freq) round(freq/f_res)*f_res;
            
            ch_len = round(fs*(ch_i(2) - ch_i(1)));
            
            for S = 1:split_num
            
                chirps = zeros(round(fs*(ch_i(2) - ch_i(1))), (length(ch_i)/2 - 1));
                spec = zeros(ax_en - ax_st + 1, length(ch_i)/2 - 1);

                % Iterate over each chirp
            
                for m = (S-1)*(length(ch_i)/2) + 1:S*(length(ch_i)/2 - 1)
                    % Isolate chirps
                    chirps(:,m - (S-1)*(length(ch_i)/2)) = proc_sig(round(fs*ch_i(m)):round(fs*ch_i(m)) + ch_len - 1);
    %                 plot(chirps(:,m))
                    spec_temp = fft(chirps(:,m - (S-1)*(length(ch_i)/2)), N);
                    spec(:,m - (S-1)*(length(ch_i)/2)) = spec_temp(ax_st:ax_en,:);
                end
                % fft and truncate spectrogram
                spec = abs(spec);

                new_row = ...
                    {{spec}, {fAx}, f_to_ind, f_round, proc_i, sig_i, 0, 0};
                obj.Spectrograms = [obj.Spectrograms; new_row];
            end
            
        end
        function obj = SpecExtract(obj, f_res, proc_i)
            
            proc_sig = obj.ProcSigs{proc_i, 'processed_signal'}{1,1};
            sig_i = obj.ProcSigs{proc_i,'signal_number'};
            
            fs = obj.Signals{sig_i,'fs'};
            ch_i = obj.Signals{sig_i,'chirp_starts'}{1,1};

            %Desired frequency resolution
            fAx_full = 0:f_res:fs - f_res;
            % find the frequency axis start/end bins on the full frequency axis
            ax_st = round(obj.Signals{sig_i,'f_st'}/f_res) + 1;
            ax_en = round(obj.Signals{sig_i,'f_en'}/f_res) + 1;
            % Truncate the frequency axis to just cover frequencies in the chirp
            fAx = fAx_full(ax_st:ax_en);
            % Calcultate the required FFT length to achieve f_res freq resolution
            N = fs/f_res;
            chirps = zeros(round(fs*(ch_i(2) - ch_i(1))), length(ch_i) - 1);
            
          % Find frequency resolution and starting frequency from fAx
            f_st = fAx(1);
            f_to_ind = @(freq) round(1 + (freq - f_st)/f_res);
            f_round = @(freq) round(freq/f_res)*f_res;
            
            spec = zeros(ax_en - ax_st + 1, length(ch_i) - 1);
            
            ch_len = round(fs*(ch_i(2) - ch_i(1)));
            
            % Iterate over each chirp
            for m = 1:length(ch_i) - 1
                % Isolate chirps
                chirps(:,m) = proc_sig(round(fs*ch_i(m)):round(fs*ch_i(m)) + ch_len - 1);
%                 plot(chirps(:,m))
                spec_temp = fft(chirps(:,m), N);
                spec(:,m) = spec_temp(ax_st:ax_en,:);
            end
            % fft and truncate spectrogram
            spec = abs(spec);
            
            new_row = ...
                {{spec}, {fAx}, f_to_ind, f_round, proc_i, sig_i, 0, 0};
            obj.Spectrograms = [obj.Spectrograms; new_row];
        end
        function obj = Update(obj, Signals, ProcSigs)
            obj.Signals = Signals;
            obj.ProcSigs = ProcSigs;
        end
        function obj = PrelimAnly(obj, proc_i)
            % Plot the spectrum so the user can provide an initial freq est
            hold off
            obj.PlotSpec(proc_i)
            disp('click on two conectutive peaks/troughs in the LFO')
            [x,~] = ginput(2); 
            
            % Calculate the initial frequency estimate
            dx = abs(x(2) - x(1));
            
            init_T = dx*obj.Signals{...
                obj.ProcSigs{proc_i, 'signal_number'}, 'ch_spc'}/1000;
            obj.Initial_f = 1/init_T;
            
            % Get user to select min/max frequencies of the notch to track
            disp('click on max and min frequencies of the notch to track')        
            [~,y] = ginput(2);
            obj.Max_f = max(y);
            obj.Min_f = min(y);     
        end
        function PlotSpec(obj, i)
%             i = find(obj.Spectrograms{:,'processed_sig'} == proc_i);
            imagesc(obj.Spectrograms{i,'spectrogram'}{1,1},...
                'YData', (obj.Spectrograms{i,'fAx'}{1,1}));
            set(gca,'YDir', 'normal')
            hold off
        end
        % LFOTrack function tracks the notch/peak in the spectrogram
        % type == 1 for notch, 2 for peak
        function obj = LFOTrack(obj, spec_num, smooth_f, type)

            if ~obj.Spectrograms{spec_num, "notch_start_fr"}
                [st_fr, st_i] = obj.LFOTrackInit(spec_num, type);
                obj.Spectrograms{spec_num, "notch_start_fr"} = st_fr;
                obj.Spectrograms{spec_num, "notch_start_i"} = st_i;
            end
            
            n_ind = obj.Spectrograms{spec_num, "notch_start_i"};
            spec = obj.Spectrograms{spec_num, "spectrogram"}{1,1};
            % Track the LFO through before the selected notch
            [backwards_f] = obj.NotchTracker...
                (spec_num, fliplr(spec(:,1:n_ind-1)), smooth_f, type);
            [forwards_f] = obj.NotchTracker...
                (spec_num, spec(:,n_ind+1:end), smooth_f, type);
            
            sig_i = obj.Spectrograms{spec_num, 'sig_num'};
            grp_d = obj.Signals{sig_i, 'chirp_offset'}{1,1};
            ch_sts = obj.Signals{sig_i, 'chirp_starts'}{1,1};
            
            LFO_freqs = [flip(backwards_f);...
                obj.Spectrograms{spec_num, "notch_start_fr"}; forwards_f];
            LFO_t = ch_sts(1:end-1)' + grp_d(2*pi*LFO_freqs);

            new_row = {{LFO_freqs}, {LFO_t}, smooth_f, spec_num};
            obj.Measured_LFOs = [obj.Measured_LFOs; new_row];
        end
        function obj = LFOTrackSplit(obj, spec_num, smooth_f, type, splits, split_num)

            if ~obj.Spectrograms{spec_num, "notch_start_fr"}
                [st_fr, st_i] = obj.LFOTrackInit(spec_num, type);
                obj.Spectrograms{spec_num, "notch_start_fr"} = st_fr;
                obj.Spectrograms{spec_num, "notch_start_i"} = st_i;
            end
            
            n_ind = obj.Spectrograms{spec_num, "notch_start_i"};
            spec = obj.Spectrograms{spec_num, "spectrogram"}{1,1};
            % Track the LFO through before the selected notch
            [backwards_f] = obj.NotchTracker...
                (spec_num, fliplr(spec(:,1:n_ind-1)), smooth_f, type);
            [forwards_f] = obj.NotchTracker...
                (spec_num, spec(:,n_ind+1:end), smooth_f, type);
            
            sig_i = obj.Spectrograms{spec_num, 'sig_num'};
            grp_d = obj.Signals{sig_i, 'chirp_offset'}{1,1};
            
            ch_sts = obj.Signals{sig_i, 'chirp_starts'}{1,1};
            ch_sts = ch_sts((split_num-1)*(length(ch_sts)/splits) + 1:(split_num)*(length(ch_sts)/splits));
            
            LFO_freqs = [flip(backwards_f);...
                obj.Spectrograms{spec_num, "notch_start_fr"}; forwards_f];
            LFO_t = ch_sts(1:end-1)' + grp_d(2*pi*LFO_freqs);

            new_row = {{LFO_freqs}, {LFO_t}, smooth_f, spec_num};
            obj.Measured_LFOs = [obj.Measured_LFOs; new_row];
        end
        function [no_fr, time_ind] = LFOTrackInit(obj, spec_num, type)
            % Plot the spectrum and ask user to select where to start track
            spectrogram = obj.Spectrograms{spec_num,'spectrogram'}{1,1};
            fAx = obj.Spectrograms{spec_num,'fAx'}{1,1};
            f_round = obj.Spectrograms{spec_num,'f_round'}{1,1};
            f_to_ind = obj.Spectrograms{spec_num,'f_to_ind'}{1,1};
            
            obj.PlotSpec(spec_num);
            disp('click on a clearly defined point in the notch')
            
            % Round to nearest chirp
            [time_ind,~] = ginput(1); 
            time_ind = round(time_ind);
            spectrogram = spectrogram(:,time_ind);
            
            % Plot the FR of that chirp and ask user to click on the notch
            plot(fAx, spectrogram)
            disp('click on the notch')
            [x,~] = ginput(1);
            x = f_round(x);
            % Set percentage range around user selection for notch search
            perc = 5;
     
            % Set range in hz
            x_range = [x*(1-perc/100), x*(1+perc/100)];
            % Convert range to indices on the fAx
            x_ind = f_to_ind(x_range);
            x_ind = [max(x_ind(1), 1), min(x_ind(2), length(fAx))];
            
            % Find the notch (spectrum min/max over the chosen range)
            if type == 1
                [~,no_i] = min(spectrogram(x_ind(1):x_ind(2)));
            elseif type == 2
                [~,no_i] = max(spectrogram(x_ind(1):x_ind(2)));
            end
            no_fr = fAx(x_ind(1) + no_i - 1);
        end
        function obj = LFOTrackInitAuto(obj, spec_num, type, st_f, en_f)
            % Plot the spectrum and ask user to select where to start track
            spectrogram = obj.Spectrograms{spec_num,'spectrogram'}{1,1};
            fAx = obj.Spectrograms{spec_num,'fAx'}{1,1};
            f_round = obj.Spectrograms{spec_num,'f_round'}{1,1};
            f_to_ind = obj.Spectrograms{spec_num,'f_to_ind'}{1,1};
            
            % Round to nearest chirp
            spectrogram = spectrogram(:,1);
            

            % Convert range to indices on the fAx
            x_ind = f_to_ind([st_f, en_f]);
            x_ind = [max(x_ind(1), 1), min(x_ind(2), length(fAx))];
            
            % Find the notch (spectrum min/max over the chosen range)
            if type == 1
                [~,no_i] = min(spectrogram(x_ind(1):x_ind(2)));
            elseif type == 2
                [~,no_i] = max(spectrogram(x_ind(1):x_ind(2)));
            end
            no_fr = fAx(x_ind(1) + no_i - 1);
            obj.Spectrograms{spec_num, "notch_start_fr"} = no_fr;
            obj.Spectrograms{spec_num, "notch_start_i"} = 1;
            
        end
        function [ntch_fr] = NotchTracker(obj, spec_num, spec, smo_f, type)
            
            A = (obj.Max_f - obj.Min_f);
            C = obj.Min_f;
            lfo_f = obj.Initial_f;
            
            J_fac = obj.Ntch_Jump_Fac;
            
            prev_f = obj.Spectrograms{spec_num, "notch_start_fr"};
            
            sig_i = obj.Spectrograms{spec_num, "sig_num"};
            dt = obj.Signals{sig_i, 'ch_spc'}/1000;
            grp_d = obj.Signals{sig_i, 'chirp_offset'}{1,1};
            fAx = obj.Spectrograms{spec_num, "fAx"}{1,1};

            f_st = fAx(1);
            f_res = fAx(2) - fAx(1);
            f_to_ind = @(f) round(1 + (f - f_st)/f_res);
            
            if type == 1
                func = @min;
            elseif type == 2
                func = @max;
            end
            
            min_f = obj.Min_f;
            max_f = obj.Max_f;

            ntch_fr = zeros(size(spec,2), 1);
            
%             lfo_w = lfo_w*1.2;

%             df = 0;
%             if dt > 0.015
%                 J_fac = 1.5;
%             end
%             if dt > 0.04
%                 J_fac = 1;
%             end
            for chirp = 1:length(ntch_fr)
                ch_spec = smooth(spec(:,chirp), smo_f);
                
            % Based on the inital LFO estimate, find the sine argument
%                 arg = asin((prev_w - C)/A);

                f_case1 = prev_f + 2*dt*A*lfo_f*pi*J_fac;
                f_case2 = prev_f - 2*dt*A*lfo_f*pi*J_fac;
%                 f_case1 = prev_f + prev_f*0.1 + df;
%                 f_case2 = prev_f - prev_f*0.1 + df;
% %                 
%                 if f_case1 > max_f
%                     f_case_diff = max_f - f_case1; 
%                     f_case1 = max_f + f_case_diff;
%                     f_case2 = max_f;
%                 end
                
%                 if f_case2 < min_f
%                     f_case_diff = min_f - f_case2; 
%                     f_case2 = min_f + f_case_diff;
%                     f_case1 = min_f;
%                 end
%                 
%                 plot(fAx, ch_spec)
%                 xline(f_case1)
%                 xline(f_case2)
%                 xline(max_f, 'color', 'r')
%                 xline(min_f, 'color', 'r')



                x_range = sort([f_case1, f_case2]);

                x_range = [max(x_range(1), min_f),...
                    min(x_range(2), max_f)];

                x_ind = f_to_ind(x_range);
                [~,i] = func(ch_spec(x_ind(1):x_ind(2)));
                
                df = fAx(x_ind(1) + i - 1) - prev_f; 
                prev_f = fAx(x_ind(1) + i - 1);

                ntch_fr(chirp) = prev_f;
            end
%             ntch_t = ch_off(ntch_fr);
%             ntch_fr = ntch_fr/(2*pi);
        end
        
        function obj = BatchSpecExtract(obj, f_res, proc_sig)
            
            fs = obj.Signals{1,'fs'};

            %Desired frequency resolution
            fAx_full = 0:f_res:fs - f_res;
            % find the frequency axis start/end bins on the full frequency axis
            ax_st = round(obj.Signals{1,'f_st'}/f_res) + 1;
            ax_en = round(obj.Signals{1,'f_en'}/f_res) + 1;
            % Truncate the frequency axis to just cover frequencies in the chirp
            fAx = fAx_full(ax_st:ax_en);
            % Find frequency resolution and starting frequency from fAx
            f_st = fAx(1);
            f_to_ind = @(freq) round(1 + (freq - f_st)/f_res);
            f_round = @(freq) round(freq/f_res)*f_res;
            
                        % Calcultate the required FFT length to achieve f_res freq resolution
            N = fs/f_res;
            
            
            for n = 1:size(obj.Signals,1)
                ch_i = obj.Signals{n,'chirp_starts'}{1,1};
                chirps = zeros(round(fs*(ch_i(2) - ch_i(1))), length(ch_i)-1);
                
                spec = zeros(ax_en - ax_st + 1, length(ch_i) - 1);
                % Iterate over each chirp
                for m = 1:length(ch_i) - 1
                    % Isolate chirps
                    chirps(:,m) = proc_sig(round(fs*ch_i(m)):round(fs*ch_i(m+1)) - 1);
                    spec_temp = fft(chirps(:,m), N);
                    spec(:,m) = spec_temp(ax_st:ax_en,:);
%                     plot(chirps(:,m))
                end
                % fft and truncate spectrogram
                spec = abs(spec);
                
                new_row = ...
                {{spec}, {fAx}, f_to_ind, f_round, 1, n, 0, 0};
                obj.Spectrograms = [obj.Spectrograms; new_row];
                
            end
           
        end
    end
    
    methods (Access = 'protected', Static = true)

    end

        
    methods (Access = 'public', Static = true)

    end
end