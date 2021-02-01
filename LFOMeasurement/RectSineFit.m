%% Rectified Sine Wave Fit
% Function that takes an irregularly sampled rectified sine wave and fits
% a rectified sine wave to it, estimating the phase offset, frequency,
% amplitude and a DC offset
% by Alec Wright

%% Arguments
% proc_sig is the test signal after being processed by the audio effect
% method selects the method, 1 = linear swept sine, 2 = exp swept sine

% ch_len is chirp length (ms) and ch_spc is chirp spacing (ms),
% f_st and f_en are the chirp start and end frequencies in hz
% T (s) is total signal time (s) and fs is sample rate

function [A, B, C, f] = RectSineFit(y, tAx, init_f)
    
% Set the grid search starting bounds, 
% init_f*(1-freq_bound) to init_f*(1+freq_bound)
    freq_bound = 0.5;
    upper_f = init_f*(1+freq_bound);
    lower_f = init_f*(1-freq_bound);
% Set number of steps in the grid search, and generate test frequencies
    steps = 40;
    test_freqs = linspace(lower_f, upper_f, steps);
% Initialise the global best, local best, and error matrices
    global_best_error = 10e12;
    local_best_error = [10e12];
    error_mat = [];
    
    % Set loop to iterate over, decreasing the grid boundaries each time
    while 1
        round_best = 1e12;
    % Find the grid frequency step size
        f_step = test_freqs(2) - test_freqs(1);
    % Iterate over each frequency in test_freqs
        for n = 1:length(test_freqs)
        % Run least squares to predict Rectified Sine Wave parameters
            [x, y_pred] = RectifiedSineFit(y, tAx, test_freqs(n));
            error = mean((y - y_pred).^2);
        % add frequency and corresponding error to the error_matrix
            error_mat(end+1, 1) =  test_freqs(n);
            error_mat(end, 2) =  error;
            % if this is a new global best error, save the freq and params
            if global_best_error > error
                global_best_error = error;
                best_freq = test_freqs(n);
                best_x = x;
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
%         hold off
%         error_mat = sortrows(error_mat);
%         plot(error_mat(:,1), error_mat(:,2), '-o')

        % If this grid yielded no improvement in performance, end the loop
        if local_best_error(end) >= local_best_error(end-1)
            break
        end
        
        % Set the test frequencies for the next grid
        upper_f = round_freq + f_step;
        lower_f = round_freq - f_step;
        test_freqs = linspace(lower_f, upper_f, steps);
    end
    
    % Return the parameters of the winning rectified sine wave
    A = best_x(1);
    B = best_x(2);
    C = best_x(3);
    f = best_freq;
end





function [x, y_pred] = RectifiedSineFit(y, tAx, f_prev)

    w_prev = pi*f_prev;
    alph_prev = (max(LFO_Measured(:,1)) - min(LFO_Measured(:,1)))*2;
    phi_prev = 0.5;
    A_prev = alph_prev*sin(phi_prev);
    B_prev = alph_prev*cos(phi_prev);
    
    
    D = ones(size(LFO_Measured, 1), 3);

    phases = mod(w_prev*tAx + phi_prev, 2*pi);

    D(phases<pi,1) = sin(w_prev*tAx(phases<pi));
    D(phases<pi,2) = cos(w_prev*tAx(phases<pi));
    D(phases>=pi,1) = -sin(w_prev*tAx(phases>pi));
    D(phases>=pi,2) = -cos(w_prev*tAx(phases>pi));

    x = D\y;

    A_prev = x(1);
    B_prev = x(2);

    y_pred = abs(A_prev*sin(w_prev*tAx) + B_prev*cos(w_prev*tAx)) + x(3);
end


% This was a function that just fits a sinewave to the data

% function [x] = SineFit(y, f_prev)
% 
%     y = LFO_Measured(:,1);
%     tAx = LFO_Measured(:,2);
% 
%     w_prev = 2*pi*f_prev;
%     alph_prev = (max(LFO_Measured(:,1)) - min(LFO_Measured(:,1)));
%     A_prev = 0;
%     B_prev = alph_prev;
%     
% %     x0 = [A_prev, B_prev, C]'
%     
%     D = ones(size(LFO_Measured, 1), 4);
%     
%     while 1
%     
% 
%         D(:,1) = cos(w_prev*tAx);
%         D(:,2) = sin(w_prev*tAx);
%         D(:,4) = -A_prev*tAx.*sin(w_prev*tAx) + B_prev*tAx.*cos(w_prev*tAx);
% 
%         x = D\y;
%         
%         A_prev = x(1);
%         B_prev = x(2);
%         w_prev = w_prev + x(4);
%         
%         y_pred = A_prev*cos(w_prev*tAx) + B_prev*sin(w_prev*tAx) + x(3);
%         
%         hold off
%         plot(y)
%         hold on
%         plot(y_pred)
%         
%     end
%     
% end
