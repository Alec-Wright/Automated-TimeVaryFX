function [ LFO_single ] = generate_LFO( speed, fs, LFO_switch )
%GENERATE_LFO Generate a single period of LFO-signal according to the parameters

%Target time vector for a single period of LFO:
frequency = 0.069*exp(0.040*speed);
len = ceil((1/frequency)*fs);

c1 = zeros(len,1)-0.89;               % First allpass filter coefficient is constant
%Choose correct coefficients depending on the LFO_switch:
if LFO_switch == 1
    lfo = sawtooth(2*pi*frequency*(1:len)/fs,0.5);  % Triangular LFO
    c2max = -0.39;                              % Largest c1 value
    c2min = -0.84;                             % Smallest c1 value
    lfo = c2min + abs(c2max-c2min)*(lfo+1)/2;  
    c2 = lfo;                                  % c2 is modulated with LFO
else % Default case ('off')
    lfo = abs(sin(2*pi*frequency*(1:len)/fs/2));  % Rectified sine wave
    c2max = 0.77;                             % Largest c1 value
    c2min = -0.49;                            % Smallest c1 value
    lfo = c2min + abs(c2max-c2min)*lfo;  
    c2 = lfo;                                % c2 is modulated with LFO
end

LFO_single = zeros(len,10);
LFO_single(:,1) = c1;
LFO_single(:,2) = c1;
LFO_single(:,3) = c2;
LFO_single(:,4) = c2;
LFO_single(:,5) = c2;
LFO_single(:,6) = c2;
LFO_single(:,7) = c2;
LFO_single(:,8) = c2;
LFO_single(:,9) = c1;
LFO_single(:,10) = c1;


end

