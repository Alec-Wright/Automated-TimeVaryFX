classdef SignalProcessor
    properties
        ProcessedSignals
        ProcessedAudio
        PedalName
        Rate
        Settings
        SaveLoc
        
    end
    methods
        function obj = SignalProcessor(PedalName, Rate)
            obj.PedalName = PedalName;
            obj.Rate = Rate;
           
            table_header = [["processed_signal", "cell"]; ...
                        ["signal_number", "int16"]; ...
                        ["rate", "double"]];
                       
            
            % Make table using fieldnames & value types from above
            obj.ProcessedSignals = ...
                table('Size',[0,size(table_header,1)],... 
                'VariableNames', table_header(:,1),...
                'VariableTypes', table_header(:,2));
            
            table_header = [["processed_audio", "cell"]; ...
            ["preceeding_sig", "cell"]];
                       
            % Make table using fieldnames & value types from above
            obj.ProcessedAudio = ...
                table('Size',[0,size(table_header,1)],... 
                'VariableNames', table_header(:,1),...
                'VariableTypes', table_header(:,2));
            
        end
        function obj = SigProc(obj, tst_sig, rate, SNR, fs)
            
            new_row = {{}, tst_sig(end), rate, SNR, {}, {}};
            obj.ProcessedSignals = [obj.ProcessedSignals; new_row];
            
            if obj.Digi == 1
                [proc_sig, LFO_real] = feval(obj.PedalName,...
                                    tst_sig(1:end-1), rate,fs);               
                obj.ProcessedSignals{end,'LFO_real'} = {LFO_real};
            % If an SNR was given, add noise equivalent to SNR
                if SNR
                    S_P = sum(proc_sig.^2)./length(proc_sig);
                    N_P = S_P*10^(-SNR/10);
                    noise = randn(length(proc_sig),1)*sqrt(N_P);
                    proc_sig = proc_sig + noise;
                    obj.ProcessedSignals{end,'nois_sig'} = {noise};
                end
                obj.ProcessedSignals{end,'processed_signal'} = {proc_sig};
            elseif obj.Digi == 2
                
            else
            % Add the bit where the test signal is saved to file so it can
            % be processed by the pedal
                audiowrite(strcat(obj.SaveLoc, '-input.wav'), tst_sig(1:end-1), fs)
                disp('process file and add -output to end of file name')
                pause
            end
        end
        function obj = LoadFromFile(obj, loc, Signals)
            
            in = audioread(strcat(loc ,'-input.wav'));
            [out, fs] = audioread(strcat(loc ,'-output.wav'));
            out = 0.95*out./max(abs(out));
            
            b = find(abs(out(1:44100*5)) > 1e-2, 1);
            
            if b < fs/2
                out = [zeros((fs/2)-b, 1); out];
            else
                out = out(b - (fs/2):end);
            end
            
            if length(out) > length(in)
                out = out(1:length(in));
            elseif length(out) < length(in)
                out = [out; zeros(length(in) - length(out), 1)];
            end
            
            st = 0;
            for n = 1:size(Signals.Signals,1)
                T = Signals.Signals{n,'T'};
                fs = Signals.Signals{n,'fs'};
                en = st + T;
                
                new_row = {{out((st*fs) + 1:en*fs)}, n, obj.Rate};
                obj.ProcessedSignals = [obj.ProcessedSignals; new_row];
                
                if Signals.Signals{n,'chunk_len'} > 0 
                    aud_st = en;
                    aud_en = en + Signals.Signals{n,'chunk_len'};
                    new_row = {{out((aud_st*fs) + 1:aud_en*fs)}, n};
                    obj.ProcessedAudio = [obj.ProcessedAudio; new_row];
                end
                
                st = en + Signals.Signals{n,'chunk_len'};
            end
            
        end
            
        function signal = SigGet(obj, sig_num)
            signal = obj.ProcessedSignals{sig_num,'processed_signal'}{1,1};
        end
    end
    methods (Access = 'public', Static = true)
        function [out] = SigLoad(loc, channels)
            in = audioread(strcat(loc ,'-input.wav'));
            out = audioread(strcat(loc ,'-output.wav'));
            b =  find( out(1:44100*2,2) > 1e-3, 1 );

            out = out(b - 22050 + 1:b - 22050 + length(in),1);
        end
    end
end