function [signal_phasered, LFO] = PhaserSweetTone(signal, speed, fs)
    
    [signal_phasered, LFO] = PhaserSweetToneSub(signal, fs, speed, 0, 0);
    
%     LFO = LFO(:,4);
    
end

function [signal_phasered, LFO] = PhaserSweetToneSub(signal, fs, speed, LFO_switch, feedback)
%Phaser model for DAFX 16 paper "Time-variant Gray-box modeling of
%a phaser pedal".
%   This function models Fame Sweet Tone Phaser pedal according to the
%   measurements, as described in the aforementioned paper.

%Inputs:
%   signal      input signal to be processed
%   fs          sampling frequency
%   speed       the value of the "Speed"-potentiometer, between 0 and 100
%   LFO_switch  if the LFO-switch is toggled on or off. Value 0 or 1.
%   feedback    if the feedback-switch is on or off. Value 0 or 1.
%
%Output:
%   signal_phasered     the processed signal

%Default values:
g = 0.5;    %This corresponds to the case g = w = 0.5
fb_n = 2;   %By changing this value, the allpass filter to which feedback
            %returns is changed. This changes the tone of the feedback etc.
            %Can accept values between 1 and 10.

%If feedback-switch is on:
if feedback == 1
    fb = 0.1;   %Change this value if different amounts of feedback is wanted
else
    fb = 0;     %If no feedback is wanted, then it is set to 0.
end

%Generate a single period of LFO
LFO_single = generate_LFO(speed, fs, LFO_switch);

%Now we place multiple of the single period LFO-signals after one another
%to generate long LFO-signal that can be used to model the phaser.
LFO = LFO_single;

while length(LFO(:,1)) < length(signal)*1.2  %We want to make sure it's long enough
   LFO = [LFO; LFO_single]; 
end

% LFO = LFO(randi(round(fs/speed)):end,:);
% LFO = LFO(1:round(length(signal)*1.05),:);
% plot(LFO(:,5))

%Then pass the input signal through a phaser:
signal_phasered = phasing_algorithm(signal, fs, g, fb, LFO, fb_n);

end

function [ LFO_single ] = generate_LFO( speed, fs, LFO_switch )
%GENERATE_LFO Generate a single period of LFO-signal according to the parameters

%Target time vector for a single period of LFO:
% frequency = 0.069*exp(0.040*speed);
frequency = speed;
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

function [signal_phasered, coeff] = phasing_algorithm(signal, fs, g, feedback, LFO, fb_n)
%PHASER_WITH_FEEDBACK This is an implementation of a simple digital phaser for
%seeing how feedback changes the amplitude spectrum.
%   This phaser consists of all-pass filters in cascade and a feedback
%   loop from the end of the all-pass filters to somewhere in the all-pass 
%   chain. 
%
%Parameters:
%   signal      the signal to be phasered (must ve a column vector)
%   fs          sampling frequency
%   g           the amount of original and filtered signal at the output
%   feedback    the amount of feedback in the system
%   LFO         the coefficients for each time
%   fb_n        to which all-pass filter the feedback-loop returns, such
%               that 1 returns to the first, 2 to second etc.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

coeff = LFO;
w = 1 - g;

%Generate a time vector (that is overly long)
time = (1:length(signal))/fs;

%Calculate a rough tau value to generate long enough vectors
tau = 4*abs(coeff(1,1))/(1-coeff(1,1)^2);
   
signal = [signal; zeros(ceil(tau*1.3-length(signal)),1)];
signal_filtered = zeros(ceil(length(signal)+tau*1.3), 1);

%Processing:
delay_line = zeros(size(LFO,2),1);

for i = 1:length(signal)
    %Input from the signal:
    in = signal(i);
    %But if the feedback comes to the first filter:
    if fb_n == 1 && i > 1
        in = signal(i) + signal_filtered(i-1)*feedback;
    end
    %All-pass filter chain:
    for n = 1:length(LFO(1,:))
        %The "output" of the filter
        temp = delay_line(n) + in*coeff(i,n);
        %New value for the delay_line
        delay_line(n) = -temp*coeff(i,n) + in;
        %The input for the next filter:
        in = temp;
        if n == fb_n-1 && i > 1
            in = temp + signal_filtered(i-1)*feedback;
        end
    end
    %Now that all the filters are processed, the output of the chain:
    signal_filtered(i) = in;
end

%And the "phasering" itself, ie combining the filtered and original signal:
signal = [signal; zeros(length(signal_filtered)-length(signal),1)];
signal_phasered = signal_filtered*w + signal*g;

end

