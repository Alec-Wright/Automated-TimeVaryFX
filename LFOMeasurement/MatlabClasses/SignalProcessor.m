classdef SignalProcessor
    properties
        ProcessedSignals
        PedalName
        Digi
        Rate
        Settings
        SaveLoc
        
    end
    methods
        function obj = SignalProcessor(PedalName, Digi, Rate)
            obj.Digi = Digi;
            obj.PedalName = PedalName;
            obj.Rate = Rate;
           
            table_header = [["processed_signal", "cell"]; ...
                        ["signal_number", "int16"]; ...
                        ["rate", "double"]];
                    
            table_header = [table_header; ...
                           ["SNR", "double"]; ...
                           ["LFO_real", "cell"]; ...
                           ["nois_sig", "cell"]];
                       
            
            % Make table using fieldnames & value types from above
            obj.ProcessedSignals = ...
                table('Size',[0,size(table_header,1)],... 
                'VariableNames', table_header(:,1),...
                'VariableTypes', table_header(:,2));
            
        end
        function obj = SigProc(obj, tst_sig, rate, SNR, fs)
            
            new_row = {{}, tst_sig(end), rate, SNR, {}, {}};
            obj.ProcessedSignals = [obj.ProcessedSignals; new_row];
            
            if obj.Digi
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
            else
            % Add the bit where the test signal is saved to file so it can
            % be processed by the pedal
                audiowrite(strcat(obj.SaveLoc, '-input.wav'), tst_sig(1:end-1), fs)
                disp('process file and add -output to end of file name')
                pause
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
            [r, lags] = xcorr(in(1:44100*10), out(1:44100*10,2));
            [~, i] = max(r);
            out = out(-lags(i) + 1:-lags(i) + length(in),1:channels);
        end
    end
end