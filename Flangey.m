addpath('LFOMeasurement/MatlabClasses')

% Sample rate
fs = 44100;

PedalName = 'fakeFlanger';

% Set the 1st notch minimum and maximum frequencies
Ntc1_st_f = 200;
Ntc1_en_f = 800;
% Set the 1st peak minimum and maximum frequencies
Pk1_st_f = 400;
Pk1_en_f = 1600;
% Find the maximum and minimum delay line length for the required freq lims
MaxD = DelayLenNotch(Ntc1_st_f, 1);
MinD = DelayLenNotch(Ntc1_en_f, 1);
MaxD2 = DelayLenPeak(Pk1_st_f, 1);
MinD2 = DelayLenPeak(Pk1_en_f, 1);
assert(MaxD == MaxD2 & MinD == MinD2)
% Find the 
notchMax = [NotchFreq(MaxD, 1), NotchFreq(MaxD, 2), NotchFreq(MaxD, 3)];
notchMin = [NotchFreq(MinD, 1), NotchFreq(MinD, 2), NotchFreq(MinD, 3)];
peakMax = [PeakFreq(MaxD2, 1), PeakFreq(MaxD2, 2), PeakFreq(MaxD2, 3)];
peakMin = [PeakFreq(MinD2, 1), PeakFreq(MinD2, 2), PeakFreq(MinD2, 3)];

b0s = [0.95, 0.5, 0.05];
aMs = [0.05, 0.5, 0.95];

rate = 0.25;

T = 20;

ch_len = [14, 24, 34;
          14, 24, 34;
          14, 24, 34;];
ch_spc = [15, 25, 35];

JF = ones(8,3,3)*2;

% Siggy = SignalHolder();
% Proccy = SignalProcessor(PedalName, 2, rate);
% Siggy = Siggy.SigGen(1, 0, 19, 20, 20, 4000, 6, fs);
% Proccy = Proccy.SigProc(Siggy.SigGet(1), 0.25, 40, fs);
% [y,LFO] = Flange(Siggy.SigGet(1), rate, 2.5, 0.5, 0.9, 0.1, fs);
% Proccy.ProcessedSignals{1,1} = {y};
% Proccy.ProcessedSignals{1,5} = {LFO};
% Anly = SignalAnalyser(Siggy.Signals, Proccy.ProcessedSignals);
% Anly = Anly.SpecExtract(1, 1);
% Anly.PlotSpec(1)

% Measured_LFO = DelayLenPeak(AnlySig(r).Measured_LFOs.Measured_LFO{n*2,1}, o);
% AnlySig(r).Measured_LFOs.Measured_LFO{n*2,1} = Measured_LFO;
% tAx = AnlySig(r).Measured_LFOs.LFO_time_axis{n*2,1};
% tAx_samps = round(tAx*fs) + 1;
% plot(tAx, Measured_LFO, 'LineWidth', 2)




table_header = [["chp_len", "int16"]; ...
            ["chp_spc", "int16"];  ...
            ["Pk/Tgh", "string"]; ...
            ["Pk/Tgh_No", "int16"]; ...
            ["b0", "double"]; ...
            ["aM", "double"]; ...
            ["Min_Fr", "double"];...
            ["Max_Fr", "double"];...
            ["Max_Err", "double"];...
            ["Mean_Err", "double"];...
            ["MAE", "double"];...
            ["MPE", "double"];...
            ["LFO_Freq_Err", "double"];...
            ["LFO_Freq_Perc", "double"];...
            ["Measured_LFO", "cell"];...
            ["Time_Axis", "cell"];...
            ["Actual_LFO", "cell"]];

% Make table using fieldnames & value types from above
Results = table('Size',[0,size(table_header,1)],... 
    'VariableNames', table_header(:,1),...
    'VariableTypes', table_header(:,2));

mults = [0.5, 1.25; 0.75, 1.1; 0.8, 1.1];

for o = 1:length(notchMax)-1
    
    st_f_ntch = notchMax(o)*mults(o,1);
    en_f_ntch = notchMin(o)*mults(o,2);

    st_f_peak = peakMax(o)*mults(o,1);
    en_f_peak = peakMin(o)*mults(o,2);

    for r = 1:length(b0s)
        b0 = b0s(r);
        aM = aMs(r);

        % Create the Signal Holder and Processor objects
        Signals(r) = SignalHolder();
        ProcSigs(r) = SignalProcessor(PedalName, 2, rate);
        AnlySig(r) = SignalAnalyser(Signals(r).Signals, ProcSigs(r).ProcessedSignals);
        AnlySig(r).Initial_f = rate;


        for n = 1:length(ch_spc)
            
            
            Signals(r) = Signals(r).SigGen(1, 0, ch_len(r,n), ch_spc(n), st_f_ntch, en_f_ntch, T, fs);
            Signals(r) = Signals(r).SigGen(1, 0, ch_len(r,n), ch_spc(n), st_f_peak, en_f_peak, T, fs);

            ProcSigs(r) = ProcSigs(r).SigProc(Signals(r).SigGet(n*2 - 1), rate, 40, fs);
            [y,LFO] = Flange(Signals(r).SigGet(n*2 - 1), rate, MaxD, MinD, b0, aM, fs);
            ProcSigs(r).ProcessedSignals{n*2 - 1,1} = {y};
            ProcSigs(r).ProcessedSignals{n*2 - 1,5} = {LFO};
            LFO_st = LFO(0.5*fs + 1); 
            LFO_en = LFO(round((0.5 + (ch_len(r,n)/1000))*fs) + 1);
            Ntch = [NotchFreq(LFO_st, o), NotchFreq(LFO_en, o)];


            ProcSigs(r) = ProcSigs(r).SigProc(Signals(r).SigGet(n*2), rate, 40, fs);
            [y,LFO] = Flange(Signals(r).SigGet(n*2), rate, MaxD, MinD, b0, aM, fs);
            ProcSigs(r).ProcessedSignals{n*2,1} = {y};
            ProcSigs(r).ProcessedSignals{n*2,5} = {LFO};
            Peak = [PeakFreq(LFO_st, o), PeakFreq(LFO_en, o)];

            AnlySig(r) = AnlySig(r).Update(Signals(r).Signals, ProcSigs(r).ProcessedSignals);
            
    %          AnlySig(r).PlotSpec(n*2 - 1);
            
    %         AnlySig(r).PlotSpec(n*2);

            AnlySig(r).Min_f = notchMax(o) - 25;
            AnlySig(r).Max_f = notchMin(o) + 25;
            AnlySig(r) = AnlySig(r).SpecExtract(0.5, n*2 - 1);
            AnlySig(r) = AnlySig(r).LFOTrackInitAuto(n*2 - 1, 1, min(Ntch)*0.92, max(Ntch)*1.08);
            AnlySig(r) = AnlySig(r).LFOTrack(n*2 - 1, 1, 1);
            Measured_LFO = DelayLenNotch(AnlySig(r).Measured_LFOs.Measured_LFO{n*2 - 1,1}, o);
            AnlySig(r).Measured_LFOs.Measured_LFO{n*2 - 1,1} = Measured_LFO;
            tAx = AnlySig(r).Measured_LFOs.LFO_time_axis{n*2 - 1,1};
            tAx_samps = round(tAx*fs) + 1;
            subplot(2,1,1)
            plot(tAx, Measured_LFO)
            hold on
            plot(tAx, LFO(tAx_samps))
            hold off
            
            [lfo_p, params] = LFOFitter.SimpleSineFit(Measured_LFO, tAx, 0.49, max(Measured_LFO) - min(Measured_LFO), min(Measured_LFO));

            error = Measured_LFO - LFO(tAx_samps)';
            perc_error = abs((Measured_LFO - LFO(tAx_samps)')./LFO(tAx_samps)')*100;
            
            Freq = 100*abs((params(1)/2) - rate)/rate;

            new_row = ...
                {ch_len(r,n), ch_spc(n), "Notch", o, b0, aM,...
                    st_f_ntch, en_f_ntch,...
                    max(abs(error)), mean(error), mean(abs(error)), mean(perc_error),...
                    params(1) - (2*rate), Freq, {Measured_LFO}, {tAx}, {LFO(tAx_samps)'}};
            Results = [Results; new_row];

            AnlySig(r).Min_f =  peakMax(o) - 25;
            AnlySig(r).Max_f =  peakMin(o) + 25;
            AnlySig(r) = AnlySig(r).SpecExtract(0.5, n*2);
            AnlySig(r) = AnlySig(r).LFOTrackInitAuto(n*2, 2, min(Peak)*0.92, max(Peak)*1.08);
            AnlySig(r) = AnlySig(r).LFOTrack(n*2, 1, 2);
            Measured_LFO = DelayLenPeak(AnlySig(r).Measured_LFOs.Measured_LFO{n*2,1}, o);
            AnlySig(r).Measured_LFOs.Measured_LFO{n*2,1} = Measured_LFO;
            tAx = AnlySig(r).Measured_LFOs.LFO_time_axis{n*2,1};
            tAx_samps = round(tAx*fs) + 1;
            subplot(2,1,2)
            plot(tAx, Measured_LFO)
            hold on
            plot(tAx, LFO(tAx_samps))
            hold off
            
            [lfo_p, params] = LFOFitter.SimpleSineFit(Measured_LFO, tAx, 0.49, max(Measured_LFO) - min(Measured_LFO), min(Measured_LFO));

            error = Measured_LFO - LFO(tAx_samps)';
            perc_error = abs((Measured_LFO - LFO(tAx_samps)')./LFO(tAx_samps)')*100;
            
            Freq = 100*abs((params(1)/2) - rate)/rate;

            new_row = ...
                {ch_len(r,n), ch_spc(n), "Peak", o, b0, aM,...
                    st_f_peak , en_f_peak ,...
                    max(abs(error)), mean(error), mean(abs(error)), mean(perc_error),...
                    params(1) - (2*rate), Freq, {Measured_LFO}, {tAx}, {LFO(tAx_samps)'}};
            Results = [Results; new_row];

        end
        

        
    end
    
    SignalsComplete(o,1:r) = Signals(1:r);
    ProcSigsComplete(o,1:r) = ProcSigs(1:r);
    AnlySigComplete(o,1:r) = AnlySig(1:r);
    
end







function PlotFilt(b0, aM, Delayms, fs)
    Delay_Samps = fs*Delayms/1000;
    Whole = floor(Delay_Samps);
    Frac = mod(Delay_Samps,1);
    freqz([b0, zeros(1,Whole - 1), 1 - Frac, Frac],...
        [1, zeros(1,Whole - 1), -aM*(1 - Frac), -aM*Frac],...
        5000, fs)
end

function [notch] = NotchFreq(Delayms, Notch_Num) 
    notch = 500/Delayms;
    notch = notch*(Notch_Num*2 - 1);
end

function [peak] = PeakFreq(Delayms, Peak_Num) 
    peak = 1000/Delayms;
    peak = peak*Peak_Num;
end

function [Del] = DelayLenNotch(Notch_Freq, Notch_Num) 
    Notch_Freq = Notch_Freq./(Notch_Num*2 - 1);
    Del = 500./Notch_Freq;
end

function [Del] = DelayLenPeak(Peak_Freq, Peak_Num) 
    Peak_Freq = Peak_Freq./Peak_Num;
    Del = 1000./(Peak_Freq);
end

function [y, LFO] = Flange(x, Rate, Max_ms, Min_ms, b0, aM, fs)
    
    b0 = abs(b0);
    aM = abs(aM);

    t = 0:1/fs:(length(x)-1)/fs;
    LFO = (Max_ms-Min_ms)*abs(sin(2*pi*Rate*t)) + Min_ms;
    LFO_Samps = fs*LFO/1000;
    
    y = zeros(length(x),1);
    
    Max_D = ceil(fs*Max_ms/1000);

    y(1:Max_D) = b0.*x(1:Max_D);

    Delays = [Max_D+1:length(x)] - LFO_Samps(Max_D+1:length(x));

    vq = interp1(1:length(x),x,Delays,'cubic');
    y(Max_D+1:end) = b0*x(Max_D+1:end) + vq';

    for n = Max_D + 1:length(x)
        TargetSamp = Delays(n - Max_D);
        Whole = floor(TargetSamp);

        vq = interp1(Whole - 5: Whole + 4,y(Whole-5:Whole+4),TargetSamp,'cubic');
        y(n) = y(n) + aM*vq;
    end
    
end
