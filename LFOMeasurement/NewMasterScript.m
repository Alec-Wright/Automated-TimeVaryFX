

import SignalHolder
import SignalProcessor
import SignalAnalyser
import LFOFitter
addpath('../DigitalEffectsPedals')

% Signal Processing Info
PedalName = 'PhaserSweetTone';
Digi = 1;
Rate = [0.1, 0.3, 0.6];
SNR = [10, 20, 30, 40];

% Test Signal Parameters
ch_typ = [1,2];
ch_len = [8, 18, 28, 38];
ch_spc = [10, 20, 30, 40];
T = 30;

tfrac = [1, 0.75, 0.5, 0.25];

smoo = [1, 10, 20];

fs = 44100;
hold off



for r = 1:length(Rate)
    i = 2;
    Signals(r) = SignalHolder();
    ProcSigs(r) = SignalProcessor(PedalName, Digi, ra);
    AnlySig(r) = SignalAnalyser(Signals(r).Signals, ProcSigs(r).ProcessedSignals);
    LFOFits(r) = LFOFitter();
    
    for smo_f = 1:length(smoo)
    
        ra = Rate(r);

        Signals = Signals.SigGen(1, 28, 30, 20, 22050, 10, fs);
        ProcSigs(r) = ProcSigs(r).SigProc(Signals(r).SigGet(1), ra, 60, fs);

        AnlySig(r) = AnlySig(r).Update(Signals(r).Signals, ProcSigs(r).ProcessedSignals);
        AnlySig(r) = AnlySig(r).SpecExtract(0.5, 1);
        AnlySig(r) = AnlySig(r).PrelimAnly(1);
        st_f = 0.9*AnlySig.Min_f;
        en_f = 1.1*AnlySig.Max_f;

        Signals(r) = Signals(r).SigGen(type, 18, 20, st_f, en_f, T, fs);
        ProcSigs(r) = ProcSigs(r).SigProc(Signals(r).SigGet(i), ra, 40, fs);

        AnlySig(r) = AnlySig(r).Update(Signals(r).Signals, ProcSigs(r).ProcessedSignals);
        AnlySig(r) = AnlySig(r).SpecExtract(0.5, i);
        AnlySig(r) = AnlySig(r).LFOTrack(i, smoo(smo_f));
        
        
        i = i + 1;
    end
    
    
    
%     for type = 1:length(ch_typ)
%         for ch_le = 1:length(ch_len)
%             
%             type = ch_typ(ty);
%             
%         end
%     end
end

for snr = 1:length(SNR)






    
    
    

    

    
    LFOFits(r) = LFOFitter();
    

    
    init_f = AnlySig(r).Initial_f;
    
    i = 2;
    

                
                
                sn = SNR(snr);
                
                ch_l = ch_len(ch_le);
                ch_s = ch_spc(ch_le);
                
                Signals(r) = Signals(r).SigGen(type, ch_l, ch_s, st_f, en_f, T, fs);
                
                i = size(Signals(r).Signals,1);

                ProcSigs(r) = ProcSigs(r).SigProc(Signals.SigGet(i), ra, sn, fs);
                
                AnlySig(r) = AnlySig(r).Update(Signals(r).Signals, ProcSigs(r).ProcessedSignals);
                AnlySig(r) = AnlySig(r).SpecExtract(0.5, i);
                
                for smo_f = 1:length(smoo) 
                    AnlySig(r) = AnlySig(r).LFOTrack(i, smoo(smo_f));

                    
                end

    Signals_lite(r).Signals = 
end

Measured_LFOs = AnlySig(1).Measured_LFOs;



save('PhaserSweetTone/Measurements.mat', 'Signals', 'ProcSigs', 'Measured_LFOs')
    
LFO_i = size(AnlySig(r).Measured_LFOs, 1);
lfo =  AnlySig(r).Measured_LFOs{LFO_i, 'Measured_LFO'}{1,1};
tAx =  AnlySig(r).Measured_LFOs{LFO_i, 'LFO_time_axis'}{1,1};
%                 plot(tAx, lfo);

%                 spec_i = AnlySig(r).Measured_LFOs{LFO_i, 'spec_num'};
%                 proc_i = AnlySig(r).Spectrograms(spec_i, 'processed_sig');
%                 sig_i = AnlySig(r).Spectrograms(spec_i, 'sig_num');

for time = 1:length(LFOFitMethods)



%                         LFOFits(r) = LFOFits(r).RectSineGridSearch(lfo(tAx < tt), tAx(tAx < tt), init_f, LFO_i, time);
    LFOFits(r) = LFOFits(r).RectSineGridSearch(lfo, tAx, init_f, LFO_i, time);
%                         LFO_p = LFOFits.LFOLookup(size(LFOFits.LFOs,1), tAx(tAx < tt));
%                         plot(tAx(tAx < tt), LFO_p)
%                         hold on
%                         plot(tAx(tAx < tt), lfo(tAx < tt))
%                         hold off
end
                














