function [ signal_phasered ] = Sweet_Tone_model( signal, fs, speed, LFO_switch, feedback )
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
fb_n = 1;   %By changing this value, the allpass filter to which feedback
            %returns is changed. This changes the tone of the feedback etc.
            %Can accept values between 1 and 10.

%If feedback-switch is on:
if feedback == 1
    fb = 0.4;   %Change this value if different amounts of feedback is wanted
else
    fb = 0;     %If no feedback is wanted, then it is set to 0.
end

%Generate a single period of LFO
LFO_single = generate_LFO(speed, fs, LFO_switch);

%Now we place multiple of the single period LFO-signals after one another
%to generate long LFO-signal that can be used to model the phaser.
LFO = LFO_single;

while length(LFO(:,1)) < length(signal)*1.4  %We want to make sure it's long enough
   LFO = [LFO; LFO_single]; 
end

%Then pass the input signal through a phaser:
signal_phasered = phasing_algorithm(signal, fs, g, fb, LFO, fb_n);

end

