%% This is an LFO measurment function
% Device, is a string of the effect name
% Enter Digi = 0 if the device being modelled is a physical device,
% otherwise Digi is [SNR, LFO_f], where SNR determines how much noise to 
% add to the signal (to test the robustness of the LFO measurement) and
% LFO_f is the LFO frequency in hz

% ch_len, ch_spc are the chirp length and spacing in ms
% ch_ty is the chirp type (1 for lin and 2 for exp increasing freq)

% T is total length of the chirp train (including 0.5 s silence at the
% beginning and end) in s, fs is the sample rate



function [LFO_Measured, LFO_True] = LFOMeasurer...
                        (Device, Digi, ch_ty, ch_len, ch_spc, T, fs)
                    
%% Preliminary testing - first, the device is measured with a chirp train
% covering the whole frequency range, so the user can pick the desired
% notch to follow the LFO by, and provide an initial LFO frequency estimate

% Generate the chirp train, test_sig, ch_i is the chirp start locations
% indices, and n is the factor used in the group delay formula to achieve
% the desired chirp length whilst covering the frequency range
%     LFOAnalyser = LFOAnalyser(Device)
    [test_sig, ch_i, n] = TestSignalGen...
                          (ch_ty, ch_len, ch_spc, 20, 10000, T, fs);
                      
  if Digi
     [proc_sig, LFO_real] = feval(Device, test_sig, Digi(2), fs);
     
     % If an SNR was given, add noise to measurement equivalent to SNR
     if Digi(1)
         S_P = sum(proc_sig.^2)./length(proc_sig);
         N_P = S_P*10^(-Digi(1)/10);
         noise = randn(length(proc_sig),1)*sqrt(N_P);
         proc_sig = proc_sig + noise;
     end
     
  else
      
  end
                    
end