classdef LFOAnalyser
    properties
        PedalName
        Measurements
        Digi
    end
    methods
        function obj = LFOAnalyser(PedalName, Digi)
            obj.PedalName = PedalName;
            
            table_header = [["ch_type", "int16"]; ...
                        ["ch_len", "double"]; ...
                        ["ch_spc", "double"]; ...
                        ["f_st", "double"]; ...
                        ["f_en", "double"]; ...
                        ["T", "double"]; ...
                        ["fs", "double"]; ...
                        ["chirp_sts", "cell"]; ...
                        ["proc_sig", "cell"]; ...
                        ["rate", "double"]
                        ["spectro", "cell"]
                        ["init_f", "double"]
                        ["notch_min", "double"]
                        ["notch_max", "double"]
                        ["notch_start", "double"]
                        ["n", "double"]
                        ["grp_delay", "cell"]
                        ["Measured_LFO", "cell"]
                        ["Analysis_Results", "table"]];
            
            if Digi
                obj.Digi = 1;
                table_header = [table_header; ...
                               ["SNR", "double"]; ...
                               ["LFO_real", "cell"]; ...
                               ["nois_sig", "cell"]];
            end
            % Make table using fieldnames & value types from above
            obj.Measurements = table('Size',[0,size(table_header,1)],... 
                'VariableNames', table_header(:,1),...
                'VariableTypes', table_header(:,2));
            
        end
        function [obj, num] = SigGen(obj, method, ch_len, ch_spc, f_st, f_en, T, fs)
            
            
           switch method
                % Create anonymous functions for the group and phase delay
                case 1
                    grp_d = @(n,w) n*w;
                case 2
                    grp_d = @(n,w) n*log(w);
           end 

           analysis_header = [["smo_f", "int16"]; ...
                        ["params", "cell"]];
           analysis = table('Size',[0,size(analysis_header,1)],... 
                'VariableNames', analysis_header(:,1),...
                'VariableTypes', analysis_header(:,2));
            
           new_row = {method, ch_len, ch_spc, f_st, f_en, T, fs,...
                {}, {}, {}, 0, 0, 0, 0, 0, 0, grp_d, {}, analysis};
           
           if obj.Digi
               new_row = [new_row, {0, {}, {}}];
           end
           obj.Measurements = [obj.Measurements; new_row];
           num = size(obj.Measurements,1);
        end
        function obj = SigProc(obj, sig_num, rate, SNR)
            
            % Generate the test signal
            output = rowfun(@LFOAnalyser.TestSignalGen,...
                              obj.Measurements(sig_num,1:7),... 
                              'OutputVariableNames',...
                            {'test_sig' 'chp_sts' 'n'});
            tst_sig = cell2mat(output{1,1});
            obj.Measurements{sig_num,'chirp_sts'} = output{1,2};
            obj.Measurements{sig_num,'rate'} = {rate};
            obj.Measurements{sig_num,'SNR'} = SNR;
            obj.Measurements{sig_num,'n'} = output{1,3};
 
            if obj.Digi
                [proc_sig, LFO_real] = feval(obj.PedalName,...
                                    tst_sig, rate,...
                                    obj.Measurements{sig_num, 'fs'});               
                obj.Measurements{sig_num,'LFO_real'} = {LFO_real};
            % If an SNR was given, add noise equivalent to SNR
                if SNR
                    S_P = sum(proc_sig.^2)./length(proc_sig);
                    N_P = S_P*10^(-SNR/10);
                    noise = randn(length(proc_sig),1)*sqrt(N_P);
                    proc_sig = proc_sig + noise;
                    obj.Measurements{sig_num,'nois_sig'} = {noise};
                end
            else
            % Add the bit where the test signal is saved to file so it can
            % be processed by the pedal
            end
            obj.Measurements{sig_num,'proc_sig'} = {proc_sig};
        end
        function obj = SpecExtract(obj, sig_num, f_res)
            
            fs = obj.Measurements{sig_num, 'fs'};
            f_st = obj.Measurements{sig_num, 'f_st'};
            f_en = obj.Measurements{sig_num, 'f_en'};
            ch_i = cell2mat(obj.Measurements{sig_num, 'chirp_sts'});
            proc_sig = cell2mat(obj.Measurements{sig_num, 'proc_sig'});
          
            %Desired frequency resolution
            fAx_full = 0:f_res:fs - f_res;
            % find the frequency axis start/end bins on the full frequency axis
            ax_st = (f_st/f_res) + 1;
            ax_en = (f_en/f_res) + 1;
            % Truncate the frequency axis to just cover frequencies in the chirp
            fAx = fAx_full(ax_st:ax_en);
            % Calcultate the required FFT length to achieve f_res freq resolution
            N = fs/f_res;
            spectro = zeros(length(fAx), length(ch_i));
            
            % Iterate over each chirp
            for m = 1:length(ch_i) - 1
                % Isolate chirp and then fft it
                chrp = proc_sig(ch_i(m):ch_i(m+1) - 1);
                spec = fft(chrp, N);
                % Truncate the spectrum to only cover the chirp frequencies
                spec_shot = spec(ax_st:ax_en);
                spectro(:, m) = abs(spec_shot);
            end
            obj.Measurements{sig_num,'spectro'} = {[spectro, fAx']};
        end
        function obj = PrelimAnly(obj, sig_num)
            % Plot the spectrum so the user can provide an initial freq est
            LFOAnalyser.plot_spec(obj.Measurements{sig_num,'spectro'});
            disp('click on two conectutive peaks/troughs in the LFO')
            [x,~] = ginput(2); 
            
            % Calculate the initial frequency estimate
            dx = abs(x(2) - x(1));
            init_T = dx*obj.Measurements{sig_num,'ch_spc'}/1000;
            obj.Measurements{sig_num, 'init_f'} = 1/init_T;
            
            % Get user to select min/max frequencies of the notch to track
            disp('click on max and min frequencies of the notch to track')        
            [~,y] = ginput(2);
            obj.Measurements{sig_num, 'notch_max'} = max(y);
            obj.Measurements{sig_num, 'notch_min'} = min(y);     
        end
        function obj = LFOTrack(obj, sig_num)
            
            [spectrogram, fAx, f_to_ind, f_round] =...
                LFOAnalyser.specLoad(obj.Measurements{sig_num,'spectro'}{1,1});
            
            % Initialise matrix to hold the notch frequency and time
            LFO_traj = zeros(size(spectrogram,2),2);
            
            % Get test signal parameters
            fs = obj.Measurements{sig_num,'fs'};
            n = obj.Measurements{sig_num,'n'};
            chp_sts = obj.Measurements{sig_num,'chirp_sts'}{1,1};
            grp_delay = obj.Measurements{sig_num,'grp_delay'}{1,1};

            % Plot the spectrum and ask user to select where to start track
            LFOAnalyser.plot_spec(obj.Measurements{sig_num,'spectro'});
            disp('click on a clearly defined point in the notch')
            
            % Round to nearest chirp
            [time_ind,~] = ginput(1); 
            time_ind = round(time_ind);
            
            % Plot the FR of that chirp and ask user to click on the notch
            plot(fAx, spectrogram(:,time_ind))
            disp('click on the notch')
            [x,~] = ginput(1);
            x = f_round(x);
            % Set percentage range around user selection for notch search
            perc = 3;
            while 1
                % Set range in hz
                x_range = [x*(1-perc/100), x*(1+perc/100)];
                % Convert range to indices on the fAx
                x_ind = f_to_ind(x_range);
                x_ind = [max(x_ind(1), 1), min(x_ind(2), length(fAx))];
                % Ask user if notch search range covers the notch
                xline(fAx(x_ind(1)))
                xline(fAx(x_ind(2)))
                pr = input('does this include the notch? (any key or n)', 's');
                % If search range is too narrow, expand it
                if strcmp(pr, 'n')
                    perc = perc + 1;
                    plot(fAx, spectrogram(:,time_ind))
                else
                    break
                end
            end
            
            % Find the group delay of the initial chirp frequency
            st_off = grp_delay(n, 2*pi*obj.Measurements{sig_num,'f_st'});
            % Find the notch (spectrum min over the chosen range)
            [~,i] = min(spectrogram(x_ind(1):x_ind(2),time_ind));
            start_notch_freq = fAx(x_ind(1) + i - 1);
            % Find the group delay at the notch frequency
            gd = grp_delay(n, 2*pi*start_notch_freq);
            % convert notch start to seconds and adjust for group delay
            start_notch_time = (chp_sts(time_ind) - 1)/fs + gd - st_off;
            LFO_traj(time_ind, :) = [start_notch_freq, start_notch_time];
            obj.Measurements{sig_num,'Measured_LFO'} = {LFO_traj};
        end
        function obj = LFOTrack2(obj, sig_num, smo_f)
            % Get test signal parameters
            fs = obj.Measurements{sig_num,'fs'};
            n = obj.Measurements{sig_num,'n'};
            chp_sts = obj.Measurements{sig_num,'chirp_sts'}{1,1};
            grp_delay = obj.Measurements{sig_num,'grp_delay'}{1,1};
            st_off = grp_delay(n, 2*pi*obj.Measurements{sig_num,'f_st'});
            LFO_traj = obj.Measurements{sig_num,'Measured_LFO'}{1,1};
            time_ind = find(LFO_traj,1);
            start_notch_freq = LFO_traj(time_ind,1);
            
            [spectrogram, fAx, ~, ~] =...
                LFOAnalyser.specLoad(obj.Measurements{sig_num,'spectro'}{1,1});
            % Create initial estimate for LFO signal parameters
            A = 2*pi*(obj.Measurements{1,'notch_max'} ...
                - obj.Measurements{1,'notch_min'});
            C = 2*pi*obj.Measurements{1,'notch_min'};
            lfo_w = pi*obj.Measurements{1,'init_f'};
            % Track the LFO through before the selected notch
            [backwards_f, backwards_t] = LFOAnalyser.NotchTracker(A, C, lfo_w,...
                    fliplr(spectrogram(:,1:time_ind-1)), fAx, flip(chp_sts(:,1:time_ind-1))...
                    , start_notch_freq*2*pi, grp_delay, n, fs, st_off, smo_f);
            % Track the LFO through after the selected notch
            [bront_f, bront_t] = LFOAnalyser.NotchTracker(A, C, lfo_w,...
                    spectrogram(:,time_ind+1:end), fAx, chp_sts(:,time_ind+1:end)...
                    , start_notch_freq*2*pi, grp_delay, n, fs, st_off, smo_f);
            % Save the measured LFO Trajectory
            LFO_traj(1:time_ind-1, 1) = flip(backwards_f);
            LFO_traj(1:time_ind-1, 2) = flip(backwards_t);
            LFO_traj(time_ind+1:end, 1) = bront_f;
            LFO_traj(time_ind+1:end, 2) = bront_t;
            obj.Measurements{sig_num,'Measured_LFO'} = {LFO_traj};
            obj.Measurements{sig_num,'Analysis_Results'}{1, 'smo_f'}
        end
        function obj = RectSineFit(obj, sig_num)
            
            y = obj.Measurements{sig_num, 'Measured_LFO'}{1,1};
            tAx = y(:,2);
            y = y(:,1);
            init_f = obj.Measurements{1, 'init_f'};
            
            [A, B, C, f] = LFOAnalyser.rectSineGridSearch(y, tAx, init_f);
            
            obj.Measurements{sig_num, 'Sine_Params'} = {[A, B, C, f]};
            
        end
    end
    
    methods (Access = 'protected', Static = true)
        
        function [tst_sig, chp_sts, n]...
            = TestSignalGen(method, ch_len, ch_spc, f_st, f_en, T, fs)

        %   Amount of silence at the start and end of the test signal (s)
            st_sil = 0.5;
            en_sil = 0.5;

        %   Convert chirp length and spacing in samples
            ch_spc_s = ceil(ch_spc*fs/1000);

        %   Synthesise Chirp
            [chirp_signal, n] = ...
                LFOAnalyser.freq_design_chirp(method, ch_len, f_st, f_en, 10, fs);

        %   Create the chirp train (choo choo!)
            num_chirps = floor((T - st_sil - en_sil)*1000/ch_spc);
            if ch_spc_s > length(chirp_signal)
                chirp_signal =[chirp_signal,...
                    zeros(1, ch_spc_s - length(chirp_signal))];
            end
            % Concatenate chirps and normalise
            tst_sig = repmat(chirp_signal, [1,num_chirps]);
            tst_sig = [zeros(1,st_sil*fs), tst_sig, zeros(1, en_sil*fs)]';
            chp_sts = {st_sil*fs:ch_spc*fs/1000:T*fs - en_sil*fs - 1};
            tst_sig = {tst_sig/max(abs(tst_sig))};
        end
        
        function [chirp, n] = ...
                freq_design_chirp(method, ch_len, f_st, f_en, T, fs)

            % Generate frequency axis
            w_st = f_st*2*pi;
            w_en = f_en*2*pi;
            wAx = 2*pi*(0:1/T:(fs/2 - 1/T));

            switch method
                % Create anonymous functions for the group and phase delay
                case 1
                    grp_d = @(n,w) n*w;
                    phs_d = @(n,w) -0.5*n*w.^2;
                case 2
                    grp_d = @(n,w) n*log(w);
                    phs_d = @(n,w) -n*w.*(log(w) - 1);
            end

            % Calculate required n to make chirp ch_len ms long
            n = (ch_len/1000)/(grp_d(1,w_en) - grp_d(1,w_st));
            % Calculate phase of each frequency component
            phases = phs_d(n,wAx);
            % Calculate chirp start/end time (time when it reaches w_st/en)
            st_t = grp_d(n,w_st);
            en_t = grp_d(n,w_en);

            % Calculate frequency domain signal
            fr_dom_sig = exp(1i*phases);
            % Zero the DC component
            fr_dom_sig(1) = 0;

            % Synthesise signal in the frequency domain
            fr_dom_sig = [fr_dom_sig, 0, flip(conj(fr_dom_sig(2:end)))];
            % ifft
            chirp = ifft(fr_dom_sig);

            % truncate the signal so it starts at w_st and ends and w_en
            chirp = chirp(round(st_t*fs):round(en_t*fs));
            % Apply a fade in and fade out to the signal
            f_steps = 40;
            chirp(1:f_steps) = chirp(1:f_steps).*linspace(0, 1, f_steps);
            chirp(end-f_steps + 1:end) ...
                = chirp(end-f_steps + 1:end).*linspace(1, 0, f_steps);
        end
        
        function [notch_freq, notch_time] =  NotchTracker(A, C, lfo_w,...
                    spectrogram, fAx, chp_sts, prev_w, grp_d, n, fs, st_of, smo_f)

                dt = (chp_sts(2) - chp_sts(1))/fs;
                wAx = fAx*2*pi;
                w_st = wAx(1);
                w_res = wAx(2) - wAx(1);
                w_to_ind = @(w) round(1 + (w - w_st)/w_res);
                
                notch_freq = zeros(length(chp_sts), 1);
                
                for chirp = 1:length(chp_sts)
%                     plot(fAx, spectrogram(:,chirp))
                % Based on the inital LFO estimate, find the sine argument
                    arg = asin((prev_w - C)/A);
                % In case one, arg is < pi/2, new freq is approx one time
                % step further along the lfo
                    w_case1 = A*sin(arg + 2*dt*lfo_w) + C;
%                     xline(w_case1/(2*pi))
                % The new notch freq affects where the actual time
                % difference
                    w_case1 = A*sin(arg + 2*dt*lfo_w +...
                        grp_d(n, w_case1) - grp_d(n, prev_w)) + C;
%                     xline(w_case1/(2*pi))
                    
                    w_case2 = A*sin(pi - arg + 2*dt*lfo_w) + C;
%                     xline(w_case2/(2*pi))
                    w_case2 = A*sin(pi - arg + 2*dt*lfo_w +...
                        grp_d(n, w_case2) - grp_d(n, prev_w)) + C;
%                     xline(w_case2/(2*pi))
                    
                    x_range = sort([w_case1, w_case2]);
                    
                    x_range = [max(x_range(1), wAx(1)),...
                        min(x_range(2), wAx(end))];
                    
                    x_ind = w_to_ind(x_range);
                    [~,i] = min(smooth(spectrogram(x_ind(1):x_ind(2),chirp), smo_f));
                    prev_w = wAx(x_ind(1) + i - 1);
%                     xline(prev_w/(2*pi))
                    notch_freq(chirp) = prev_w;
                end
                notch_time = ((chp_sts - 1)/fs)' + grp_d(n, notch_freq) - st_of;
                notch_freq = notch_freq/(2*pi);
        end
    end
    methods (Access = 'public', Static = true)
        function specplot = plot_spec(spectrogram)
            spectrogram = spectrogram{1,1};
            fAx = spectrogram(:,end);
            spectrogram = spectrogram(:,1:end-1);
            specplot = imagesc(spectrogram, 'YData', fAx);
            set(gca,'YDir', 'normal')
        end
        function [spectrogram, fAx, f_to_ind, f_round] = ...
                specLoad(spectrogram)
            % Extract Spectrogram and fAx from processed measurement signal
            fAx = spectrogram(:,end);
            spectrogram = spectrogram(:,1:end-1);
            % Find frequency resolution and starting frequency from fAx
            f_res = fAx(2) - fAx(1);
            f_st = fAx(1);
            
            f_to_ind = @(freq) round(1 + (freq - f_st)/f_res);
            f_round = @(freq) round(freq/f_res)*f_res;
        end
        function [A, B, C, f] = rectSineGridSearch(y, tAx, init_f)
            % Set the grid search starting bounds, 
            % init_f*(1-freq_bound) to init_f*(1+freq_bound)
            freq_bound = 0.3;
            upper_f = init_f*(1+freq_bound);
            lower_f = init_f*(1-freq_bound);
    % Set number of steps in the grid search, and generate test frequencies
            steps = 50;
            test_freqs = linspace(lower_f, upper_f, steps);
            % Initialise the global best, local best, and error matrices
            global_best_error = 10e12;
            local_best_error = [10e12];
            error_mat = [];
            
    % Set loop to iterate over, decreasing the grid boundaries each time
            while 1
                round_best = 1e12;
            % Find the grid frequency step size
                f_step = test_freqs(2) - test_freqs(1);
            % Iterate over each frequency in test_freqs
                for n = 1:length(test_freqs)
                % Run least squares to predict Rectified Sine Wave parameters
                    [x, y_pred] = LFOAnalyser.rectifiedSineFit(y, tAx, test_freqs(n));
                    error = mean((y - y_pred).^2);
                % add frequency and corresponding error to the error_matrix
                    error_mat(end+1, 1) =  test_freqs(n);
                    error_mat(end, 2) =  error;
                    % if this is a new global best error, save the freq and params
                    if global_best_error > error
                        global_best_error = error;
                        best_freq = test_freqs(n);
                        best_x = x;
                    end
                    % If this is best in this grid search, save freq and error
                    if round_best > error
                       round_best = error; 
                       round_freq = test_freqs(n);
                    end
                end
            % Save the best error achieved on this round
                local_best_error(end+1) = round_best;

                % Optional, sort by freqs and plot the error v freqs
                hold off
                error_mat = sortrows(error_mat);
                plot(error_mat(:,1), error_mat(:,2), '-o')

                % If this grid yielded no improvement in performance, end the loop
                if local_best_error(end) >= local_best_error(end-1)
                    break
                end

                % Set the test frequencies for the next grid
                upper_f = round_freq + f_step;
                lower_f = round_freq - f_step;
                steps = steps*2;
                test_freqs = linspace(lower_f, upper_f, steps);
            end
    
            % Return the parameters of the winning rectified sine wave
            A = best_x(1);
            B = best_x(2);
            C = best_x(3);
            f = best_freq;
            
        end
        function [x, y_pred] = rectifiedSineFit(y, tAx, f_prev)

            w_prev = pi*f_prev;
            alph_prev = (max(y) - min(y))*2;
            phi_prev = 0;
            A_prev = alph_prev*sin(phi_prev);
            B_prev = alph_prev*cos(phi_prev);

            D = ones(size(y, 1), 3);

            phases = mod(w_prev*tAx + phi_prev, 2*pi);

            D(phases<pi,1) = sin(w_prev*tAx(phases<pi));
            D(phases<pi,2) = cos(w_prev*tAx(phases<pi));
            D(phases>=pi,1) = -sin(w_prev*tAx(phases>pi));
            D(phases>=pi,2) = -cos(w_prev*tAx(phases>pi));

            x = D\y;

            A_prev = x(1);
            B_prev = x(2);

            y_pred = abs(A_prev*sin(w_prev*tAx) + B_prev*cos(w_prev*tAx)) + x(3);
        end
        function [LFO_pred] = rectSineGen(A,B,C,w,tAx)
            LFO_pred = abs(A*sin(w*tAx) + B*cos(w*tAx)) + C;
        end
        function [LFO_pred] = rectSineGenPow(A,B,C,w,tAx,n)
            LFO_pred = abs(A*sin(w*tAx) + B*cos(w*tAx)) + C;
            predmax = sqrt(A^2 + B^2);
            LFO_pred = predmax.*((LFO_pred-C)./predmax).^n + C;
        end
        function [params] = LFOCorr(params, LFO_measured, smo_f)
            params(3) = min(LFO_measured);
            
            alph = max(smooth(LFO_measured, smo_f)) - min(LFO_measured);
            prev_alph = sqrt(params(1).^2 + params(2).^2);
            params(1) = params(1)*alph/prev_alph;
            params(2) = params(2)*alph/prev_alph;
        end
        function [power] = LFOReshapeFactor(A, B, C, w, LFO_measured, tAx)
            n = [0.1, 3];
            step_size = 0.1;
            besterror = 1e20;
            prevbest = besterror;
            bestpower = 1;
            while 1
                steps = n(1):step_size:n(2);
                for power = 1:length(steps)
                    LFO_pred = LFOAnalyser.rectSineGenPow(A,B,C,w,tAx,steps(power));
%                     plot(LFO_pred);
%                     hold on
%                     plot(LFO_measured);
%                     hold off
                    error = mean((LFO_pred - LFO_measured).^2);
                    if error < besterror
                        besterror = error;
                        bestpower = steps(power);
                    end
                end
                
                if besterror == prevbest
                    power = bestpower;
                    break
                end
                prevbest = besterror;
                n = [bestpower - 2*step_size, bestpower + 2*step_size];
                step_size = step_size/5;
            end
        end

        
        
        
    end
end