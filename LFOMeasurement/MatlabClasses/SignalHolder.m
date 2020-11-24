classdef SignalHolder
    properties
        Signals
        SaveLoc
        SegLen
    end
    methods
        function obj = SignalHolder()
            
            table_header = [["signal", "cell"];
                        ["ch_type", "int16"]; ...
                        ["ch_len", "double"]; ...
                        ["ch_spc", "double"]; ...
                        ["f_st", "double"]; ...
                        ["f_en", "double"]; ...
                        ["T", "double"]; ...
                        ["fs", "double"]; ...
                        ["chirp_starts", "cell"];... 
                        ["chirp_offset", "cell"]];
            
            % Make table using fieldnames & value types from above
            obj.Signals = table('Size',[0,size(table_header,1)],... 
                'VariableNames', table_header(:,1),...
                'VariableTypes', table_header(:,2));
            
        end
        function obj = SigGen(obj, method, ch_len, ch_spc, f_st, f_en, T, fs)
            
            [tst_sig, chp_sts, n] = SignalHolder.TestSignalGen...
                (method, ch_len, ch_spc, f_st, f_en, T, fs);
            
           switch method
                % Create anonymous function for the time as a funciton of f
                case 1
                    chp_off = @(w) n*w - n*f_st*2*pi;
                case 2
                    chp_off = @(w) n*log(w) - n*log(f_st*2*pi);
           end

           new_row = {tst_sig, method, ch_len, ch_spc, f_st, f_en, T,...
               fs, chp_sts, chp_off};
           
           obj.Signals = [obj.Signals; new_row];
        end
        function signal = SigGet(obj, sig_num)
            signal = [obj.Signals{sig_num,'signal'}{1,1}; sig_num];
        end
    end
        
    methods (Access = 'public', Static = true)
        function [tst_sig, chp_sts, n] = ...
                TestSignalGen(method, ch_len, ch_spc, f_st, f_en, T, fs)
        % Amount of silence at the start and end of the test signal (s)
            s_sil = 0.5;
            e_sil = 0.5;

        %   Convert chirp length and spacing into samples
            ch_spc_s = ceil(ch_spc*fs/1000);
            ch_spc = 1000*ch_spc_s/44100;

        %   Synthesise Chirp
            [chirp_signal, n] = SignalHolder.freq_design_chirp...
                                (method, ch_len, f_st, f_en, 10, fs);

        %   Create the chirp train (choo choo!)
            num_chirps = ceil((T - s_sil - e_sil)*1000/ch_spc);
            if ch_spc_s > length(chirp_signal)
                chirp_signal =[chirp_signal,...
                    zeros(1, ch_spc_s - length(chirp_signal))];
            end
            % Concatenate chirps and normalise
            tst_sig = repmat(chirp_signal, [1,num_chirps]);
            tst_sig = [zeros(1, s_sil*fs), tst_sig, zeros(1, e_sil*fs)]';
            tst_sig = tst_sig(1:T*fs);
            chp_sts = {s_sil:ch_spc_s/fs:T - e_sil};
            tst_sig = {tst_sig/max(abs(tst_sig))};
        end
        function [chirp, n] = ...
            freq_design_chirp(method, ch_len, f_st, f_en, T, fs)

            % Generate frequency axis
            w_st = f_st*2*pi;
            w_en = f_en*2*pi;
            wAx = 2*pi*(0:1/T:(fs/2 - 1/T));

            switch method
                % Create anon functions for the group and phase delay
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
            % Calculate chirp start/end time (when it reaches w_st/en)
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

            % truncate signal so it starts at w_st and ends and w_en
            chirp = chirp(round(st_t*fs):round(en_t*fs));
            % Apply a fade in and fade out to the signal
            f_steps = 20;
            chirp(1:f_steps)...
                = chirp(1:f_steps).*linspace(0, 1, f_steps);
            chirp(end-f_steps + 1:end) ...
                = chirp(end-f_steps + 1:end).*linspace(1, 0, f_steps);
        end
    end
end