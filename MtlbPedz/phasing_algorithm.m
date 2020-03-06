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

