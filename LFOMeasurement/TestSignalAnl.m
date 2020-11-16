%% Test Signal Analyser
% Function that analyses a chirp train signal thats been processed by a 
% periodically modulated time-varying audio effects, the purpose is convert
% the chirp train into a sequence of points describing the trajectory of
% the notch frequency over the measurment signal
% by Alec Wright

%% Arguments
% proc_sig is the test signal after being processed by the audio effect
% method selects the method, 1 = linear swept sine, 2 = exp swept sine

% ch_len is chirp length (ms) and ch_spc is chirp spacing (ms),
% f_st and f_en are the chirp start and end frequencies in hz
% T (s) is total signal time (s) and fs is sample rate

function [LFO_Pos, tAx]...
    = TestSignalAnl(proc_sig, ch_i, method, f_st, f_en, fs, n)

    switch method
    % Create anonymous functions for the group and phase delay
    case 1
        grp_d = @(n,w) n*w;
    case 2
        grp_d = @(n,w) n*log(w);
    end
   
    %Desired frequency resolution
    f_res = 0.1;
    fAx_full = 0:f_res:fs - f_res;
    % find the frequency axis start/end bins on the full frequency axis
    ax_st = (f_st/f_res) + 1;
    ax_en = (f_en/f_res) + 1;
    % Truncate the frequency axis to just cover frequencies in the chirp
    fAx = fAx_full(ax_st:ax_en);
    % Calcultate the required FFT length to achieve f_res freq resolution
    N = fs/f_res;

    % Find the number of samples between the start of each chirp
    ch_g = ch_i(2) - ch_i(1);
    spectro = zeros(length(fAx), length(ch_i));

    % Initialise matrix of notch frequency time when notch occurs
    LFO_Pos = zeros(length(ch_i), 1);
    tAx = zeros(size(LFO_Pos));
    
    % Calculate the group delay at the start of the chirp signal
    st_off = grp_d(n, 2*pi*f_st);

    % Iterate over each chirp
    for m = 1:length(ch_i)
        % Isolate chirp and then fft it
        chrp = proc_sig(ch_i(m):ch_i(m) + ch_g - 1);
        spec = fft(chrp, N);
        
        % Truncate the spectrum to only cover the chirp frequencies
        spec_shot = spec(ax_st:ax_en);
        
        spectro(:, m) = spec_shot;

        % Find where the notch occurs in the chirp
        [~,k] = min(abs(spec_shot));
        notch_freq = fAx(k);
        % Save notch location frequency
        LFO_Pos(m) = notch_freq;
        % Find and save the time when that frequency occured in the test
        % signal
        tAx(m) = (grp_d(n, 2*pi*notch_freq) - st_off)...
                            + (ch_i(m)/fs);
        
    end

    imagesc(flipud(abs(spectro)))
    ticklocs = round(1:length(fAx)/6:length(fAx));
    yticks(ticklocs + length(fAx) - ticklocs(end))
    yticklabels(num2cell(flip(round(fAx(ticklocs),-1))));

end
