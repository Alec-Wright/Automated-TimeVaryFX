classdef LFOFitter
    properties
        LFOs
    end
    
    methods
        function obj = LFOFitter()
            
            table_header = [["true_freq", "double"]; ...
            ["pred_freq", "double"]; ...
            ["pred_amp", "double"]; ...
            ["pred_phase", "double"]];
            
            % Make table using fieldnames & value types from above
            obj.LFOs = table('Size',[0,size(table_header,1)],... 
                'VariableNames', table_header(:,1),...
                'VariableTypes', table_header(:,2));
            
        end
        function obj = RectSineFit(obj, LFO, tAx, sig_i, proc_i, anly_i, init_f)
            
            [A, B, C, f] = LFOAnalyser.rectSineGridSearch(y, tAx, init_f);
            
            obj.Measurements{sig_num, 'Sine_Params'} = {[A, B, C, f]};
            
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
            
%             A_init = alph_init*sin(phi_init);
%             B_init = alph_init*cos(phi_init);

            D = ones(size(y, 1), 3);

            phases = mod(w_init*tAx + phi_init, 2*pi);

            D(phases<pi,1) = sin(w_init*tAx(phases<pi));
            D(phases<pi,2) = cos(w_init*tAx(phases<pi));
            D(phases>=pi,1) = -sin(w_init*tAx(phases>pi));
            D(phases>=pi,2) = -cos(w_init*tAx(phases>pi));

            x = D\y;

            A = x(1);
            B = x(2);
            C = x(3);

            y_pred = abs(A*sin(w_init*tAx) + B*cos(w_init*tAx)) + x(3);
            
            plot(tAx, y_pred)
            hold on
            plot(tAx, y)
            hold off
            
        end
    end
end