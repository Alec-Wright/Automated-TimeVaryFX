classdef SignalHolder
    properties
        Signals
    end
    methods
        function obj = SignalHolder()
            
            table_header = [["signal", "cell"];
                        ["ch_type", "int16"]; ...
                        ["ch_len", "double"]; ...
                        ["ch_spc", "double"]; ...
                        ["f_st", "double"]; ...
                        ["f_en", "double"]; ...
                        ["a", "double"]; ...
                        ["T", "double"]; ...
                        ["fs", "double"]; ...
                        ["chirp_starts", "cell"];... 
                        ["chirp_offset", "cell"];...
                        ["chunk_len", "double"]];
            
            % Make table using fieldnames & value types from above
            obj.Signals = table('Size',[0,size(table_header,1)],... 
                'VariableNames', table_header(:,1),...
                'VariableTypes', table_header(:,2));
            
        end
        
        %% Signal Generator function, adds a chirp signal to Signal list
        % Method = 1/2, generates the chirp in the frequency domain using a
        % linear/exponential phase delay function. Method = 3 uses first
        % order all pass filter freq response.
        %
        % ch_len/ch_spc is the length and space between chirp starts in ms
        % frng is chirp start/end frequency if Method = 1/2, otherwise it
        % is the all-pass filter parameter a.
        
        % T is the test singal total length (includes 0.5s silence at
        % the start/end). chunk_len is the length in seconds of the non
        % chirp signal audio that will follow the chirp signal.
        function obj = SigGen(obj, method, ch_len, ch_spc, frng, T, fs,...
                chunk_len)
            
            % subfunction that generates the chirp train, returns the chirp
            % train, and chirp start positions, and n (see below)
            [tst_sig, chp_sts, n] = SignalHolder.TestSignalGen...
                (method, ch_len, ch_spc, frng, T, fs);
            
            switch method
            % Create anonymous function which uses the group delay to
            % calculate where a frequency appears in the chirp relative to
            % the chrip start (in ms)
            % Then save parameters used to generate the chirp train
                case 1
                % n is a parameter controlling total chirp length
                    chp_off = @(w) (n*w - n*frng(1)*2*pi)';
                    new_row = {tst_sig, method, ch_len, ch_spc, frng(1),...
                    frng(2), 0, T, fs, chp_sts, chp_off, chunk_len};
                case 2
                    chp_off = @(w) (n*log(w) - n*log(frng(1)*2*pi))';
                    new_row = {tst_sig, method, ch_len, ch_spc, frng(1),...
                    frng(2), 0, T, fs, chp_sts, chp_off, chunk_len};
                case 3
                % n is the number of cascaded allpass filters
                    if frng > 0
                        chp_off = @(w) -n*(((1 - frng^2)/(1 + 2*frng*cos(w/(fs)) + frng^2)) - ...
                            ((1 - frng^2)/(1 + 2*frng*cos(0) + frng^2)))/fs;
                    else
                        chp_off = @(w) -n*(((1 - frng^2)/(1 + 2*frng*cos(w/(fs)) + frng^2)) - ...
                            ((1 - frng^2)/(1 + 2*frng*cos(pi) + frng^2)))/fs;
                    end
                    new_row = {tst_sig, method, ch_len, ch_spc, 0,...
                        22050, frng, T, fs, chp_sts, chp_off, chunk_len};
            end
           
            

           offset_ch_sts = new_row{1,10}{1,1};
           for n = 1:size(obj.Signals,1)
               offset_ch_sts = offset_ch_sts + obj.Signals{n,'T'} + obj.Signals{n,'chunk_len'};
           end
           new_row{1,10}{1,1} = offset_ch_sts;
            
           obj.Signals = [obj.Signals; new_row];
        end
        
        function signal = SigGet(obj, sig_num)
            signal = [obj.Signals{sig_num,'signal'}{1,1}; sig_num];
        end
    end
        
    methods (Access = 'public', Static = true)
        function [tst_sig, chp_sts, n] = ...
                TestSignalGen(method, ch_len, ch_spc, frng, T, fs)
        % Amount of silence at the start and end of the test signal (s)
            s_sil = 0.5;
            e_sil = 0.5;

        %   Convert chirp length and spacing into samples
            ch_spc_s = ceil(ch_spc*fs/1000);
            ch_spc = 1000*ch_spc_s/44100;

        %   Synthesise Chirp
            if method == 1 || method == 2
                assert(length(frng) == 2)
                [chirp_signal, n] = SignalHolder.freq_design_chirp...
                                    (method, ch_len, frng(1), frng(2), 10, fs);
            elseif method == 3
                assert(length(frng) == 1)
                [chirp_signal, n] = SignalHolder.fo_ap_chirp...
                    (ch_len, frng, 10, fs);
            elseif method == 4
                assert(length(frng) == 2)
                chirp_signal = [1,zeros(1,10)];
                n = 0;
            end
            
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
                    if w_st < 100
                        w_st = 100;
                    end
                    grp_d = @(n,w) n*log(w);
                    phs_d = @(n,w) -n*w.*(log(w) - 1);
            end

            % Calculate required n to make chirp ch_len ms long
            n = (ch_len/1000)/(grp_d(1,w_en) - grp_d(1,w_st));
            if method == 2
                n = 0.9*n;
            end

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
            fr_dom_sig = [fr_dom_sig,0 ,flip(conj(fr_dom_sig(2:end)))];
            % ifft
            chirp = ifft(fr_dom_sig);

            % truncate signal so it starts at w_st and ends and w_en
            chirp = chirp(max([round(st_t*fs),1]):round(en_t*fs));
            %windy = tukeywin(length(chirp),0.05);
            %windy = hann(length(chirp));
            %chirp = chirp.*windy';
        end
        
        % Generate chirp using cascaded first-order all pass signals
        function [chirp, n] = ...
            fo_ap_chirp(ch_len, a, T, fs)
            
            delhf = (1 - a^2)/(1 + 2*a*cos(pi) + a^2);
            dellf = (1 - a^2)/(1 + 2*a*cos(0) + a^2);
            
            grpdiff = max(delhf,dellf) - min(delhf,dellf);
            
            ch_len_samps = fs*ch_len/1000;
            
            M = floor(0.8*ch_len_samps/grpdiff);
            n = M;

            x = [1,zeros(1,T*44100)];

            for n = 1:M
                x = filter([a, 1], [1, a], x);
            end
           
            chirp = x(find(x>1e-3,1):find(x>1e-3,1) + floor(ch_len_samps) - 1);
        end
        
    end
end