%% Test Signal Generator
% Function that generates a chirp signal for measurement of time-varying
% systems over time
% by Alec Wright

%% Arguments
% method selects the method, 1 = linear swept sine, 2 = exp swept sine

% ch_len is chirp length (ms) and ch_spc is chirp spacing (ms),
% f_st and f_en are the chirp start and end frequencies in hz
% T (s) is total signal time (s) and fs is sample rate
function [tst_sig, chp_sts, n]...
    = TestSignalGen(method, ch_len, ch_spc, f_st, f_en, T, fs)

%   Set the amount of silence at the start and end of the test signal (s)
    st_sil = 0.5;
    en_sil = 0.5;
    
%   Convert chirp length and spacing in samples
%     ch_len_s = ceil(ch_len*fs/1000);
    ch_spc_s = ceil(ch_spc*fs/1000);

%   Synthesise Chirp
    [chirp_signal, n] = freq_design_chirp(method, ch_len, f_st, f_en, 10, fs);

%   Create the chirp train (choo choo!)
    num_chirps = floor((T - st_sil - en_sil)*1000/ch_spc);
    if ch_spc_s > length(chirp_signal)
        chirp_signal =...
            [chirp_signal, zeros(1, ch_spc_s - length(chirp_signal))];
    end
    tst_sig = repmat(chirp_signal, [1,num_chirps]);
    tst_sig = [zeros(1,st_sil*fs), tst_sig, zeros(1, en_sil*fs)]';
    chp_sts = st_sil*fs:ch_spc*fs/1000:T*fs - en_sil*fs - 1;
    
    tst_sig = tst_sig/max(abs(tst_sig));
end

function [chirp, n] = freq_design_chirp(method, ch_len, f_st, f_en, T, fs)

    % Generate frequency axis
    w_st = f_st*2*pi;
    w_en = f_en*2*pi;
    wAx = 2*pi*(0:1/T:(fs/2 - 1/T));
    
    switch method
        % Create anonymous functions for the group and phase delay
        case 1
            grp_d = @(n,w) n*w;
            phs_d = @(n,w) -0.5*n*w.^2;
        case 2
            grp_d = @(n,w) n*log(w);
            phs_d = @(n,w) -n*w.*(log(w) - 1);
    end
    
    % Calculate required n to make chirp ch_len ms long
    n = (ch_len/1000)/(grp_d(1,w_en) - grp_d(1,w_st));
    % Calculate phase of each frequency component
    phases = phs_d(n,wAx);
    % Calculate the start time of the chirp (time when it reaches w_st)
    st_t = grp_d(n,w_st);
    en_t = grp_d(n,w_en);
    
    % Calculate frequency domain signal
    fr_dom_sig = exp(1i*phases);
    % Zero the DC component
    fr_dom_sig(1) = 0;
    
    % Synthesise signal in the frequency domain
    fr_dom_sig = [fr_dom_sig, 0, flip(conj(fr_dom_sig(2:end)))];
    % ifft
    chirp = ifft(fr_dom_sig);
    
    % truncate the signal so it starts at w_st and ends and w_en
    chirp = chirp(round(st_t*fs):round(en_t*fs));
    % Apply a fade in and fade out to the signal
    f_steps = 40;
    chirp(1:f_steps) = chirp(1:f_steps).*linspace(0, 1, f_steps);
    chirp(end-f_steps + 1:end) ...
        = chirp(end-f_steps + 1:end).*linspace(1, 0, f_steps);
end
