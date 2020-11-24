classdef LFOFitter
    properties
        LFOs
    end
    
    methods
        function obj = LFOFitter()
            
            table_header = [["true_freq", "double"]; ...
            ["pred_freq", "double"]; ...
            ["pred_amp", "double"]; ...
            ["pred_phase", "double"]; ...
            ["pred_offset", "double"]; ...
            ["LFO_i", "int16"]; ...
            ["method", "double"]];
            
            % Make table using fieldnames & value types from above
            obj.LFOs = table('Size',[0,size(table_header,1)],... 
                'VariableNames', table_header(:,1),...
                'VariableTypes', table_header(:,2));
            
        end
        function obj = RectSineGridSearch...
                (obj, LFO, tAx, init_f, LFO_i, method)
            
            % Set the grid search starting bounds, 
            % init_f*(1-freq_bound) to init_f*(1+freq_bound)
            freq_bound = 0.1;
            upper_f = init_f*(1+freq_bound);
            lower_f = init_f*(1-freq_bound);
    % Set number of steps in the grid search, and generate test frequencies
            steps = 50;
            test_freqs = linspace(lower_f, upper_f, steps);
            % Initialise the global best, local best, and error matrices
            global_best_error = 10e12;
            local_best_error = [10e12];
            error_mat = [];
            
            MaxTrue = max(smooth(LFO, 5));
            MinTrue = min(LFO);
            AmpTrue = (MaxTrue - MinTrue);
            
    % Set loop to iterate over, decreasing the grid boundaries each time
            while 1
                round_best = 1e12;
            % Find the grid frequency step size
                f_step = test_freqs(2) - test_freqs(1);
            % Iterate over each frequency in test_freqs
                for n = 1:length(test_freqs)
                % Run least squares to predict Rectified Sine Wave parameters
                    [A,B,C] = LFOFitter.rectifiedSineFit(LFO, tAx, pi*test_freqs(n));
                    
                    
                    if method > 1
                        
                        prev_alph = sqrt(A.^2 + B.^2);
                        A = A*AmpTrue/prev_alph;
                        B = B*AmpTrue/prev_alph;
                        
%                         [~, phase] = LFOFitter.AmpsToSine(A,B);
%                         [A, B] = LFOFitter.SineToAmps(AmpTrue, phase);
                        C = MinTrue;
                    end
                    
                    LFO_p = LFOFitter.LFOGenAmps(A, B, C, pi*test_freqs(n), tAx);
            
                    
                    
                    
                    error = mean((LFO - LFO_p).^2);
                % add frequency and corresponding error to the error_matrix
                    error_mat(end+1, 1) =  test_freqs(n);
                    error_mat(end, 2) =  error;
                    % if this is a new global best error, save the freq and params
                    if global_best_error > error
                        global_best_error = error;
                        best_freq = test_freqs(n);
                        best_x = [A,B,C];
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
%                 plot(error_mat(:,1), error_mat(:,2), '-o')

                % If this grid yielded no improvement in performance, end the loop
                if local_best_error(end) >= local_best_error(end-1)
                    break
                end

                % Set the test frequencies for the next grid
                upper_f = round_freq + 2*f_step;
                lower_f = round_freq - 2*f_step;
                steps = steps*2;
                test_freqs = linspace(lower_f, upper_f, steps);
            end
    
            % Return the parameters of the winning rectified sine wave
            A = best_x(1);
            B = best_x(2);
            C = best_x(3);
            f = best_freq;
            
            [amp, phase] = obj.AmpsToSine(A,B);
            
            new_row = {0, f, amp, phase, C, LFO_i, method};
            obj.LFOs = [obj.LFOs; new_row];
            
        end
        function LFO = LFOLookup(obj, LFO_i, tAx)
            amp = obj.LFOs{LFO_i, 'pred_amp'};
            phase = obj.LFOs{LFO_i, 'pred_phase'};
            C = obj.LFOs{LFO_i, 'pred_offset'};
            w = obj.LFOs{LFO_i, 'pred_freq'};
            LFO = LFOFitter.LFOGenSine(amp, phase, C, pi*w, tAx);
        end
    end
    methods (Access = 'public', Static = true)
        function [amp, phase] = AmpsToSine(A,B)
            amp = sqrt(A^2 + B^2);
            phase = asin(A/amp);
        end
        function [A, B] = SineToAmps(amp,phase)
            A = amp*cos(phase);
            B = amp*sin(phase);
        end
        function [LFO] = LFOGenSine(amp, phase, C, w, tAx)
            LFO = abs(amp*sin(w*tAx + phase)) + C;
        end
        function [LFO] = LFOGenAmps(A, B, C, w, tAx)
            LFO = abs(A*cos(w*tAx) + B*sin(w*tAx)) + C;
        end
        function [A,B,C] = rectifiedSineFit(y, tAx, w_init)

            C_init = min(y);
            alph_init = (max(y) - C_init)*2;
            phi_init = asin((y(1) - C_init)/alph_init) - w_init*tAx(1);
            
            if (y(10) - y(1))  < 0
                phi_init = pi - phi_init;
            end
            
%             phi_init = 0;
            
            D = ones(size(y, 1), 3);

            phases = mod(w_init*tAx + phi_init, 2*pi);

            D(phases<pi,1) = cos(w_init*tAx(phases<pi));
            D(phases<pi,2) = sin(w_init*tAx(phases<pi));
            D(phases>=pi,1) = -cos(w_init*tAx(phases>=pi));
            D(phases>=pi,2) = -sin(w_init*tAx(phases>=pi));

            x = D\y;

            A = x(1);
            B = x(2);
            C = x(3);

%             y_pred = abs(A*cos(w_init*tAx) + B*sin(w_init*tAx)) + x(3);
%             plot(tAx, y_pred)
%             hold on
%             plot(tAx, y)
%             hold off
            
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
                    LFO_pred = LFOFitter.rectSineGenPow(A,B,C,w,tAx,steps(power));
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
                
                if besterror >= prevbest
                    power = bestpower;
                    break
                end
                prevbest = besterror;
                n = [bestpower - 2*step_size, bestpower + 2*step_size];
                step_size = step_size/5;
            end
        end
        function [LFO_pred] = rectSineGenPow(A,B,C,w,tAx,n)
            LFO_pred = abs(A*sin(w*tAx) + B*cos(w*tAx)) + C;
            predmax = sqrt(A^2 + B^2);
            LFO_pred = predmax.*((LFO_pred-C)./predmax).^n + C;
        end
    end
end