

import SignalHolder
import SignalProcessor
import SignalAnalyser
addpath('../DigitalEffectsPedals')

PedalName = 'PhaserSweetTone';
Digi = 1;
Rate = [0.3, 0.6];
SNR = [20,40];

ch_ty = [1,2];

ch_len = [8,18,28,38];
ch_spc = [10,20,30,40];

fs = 44100;

for r = 1:length(Rate)
    ra = Rate(r);
    Signals(r) = SignalHolder();
    Signals(r) = Signals(r).SigGen(1, 28, 30, 20, 22050, 10, fs);
    ProcSigs(r) = SignalProcessor(PedalName, Digi, ra);
    ProcSigs(r) = ProcSigs(r).SigProc(Signals(r).SigGet(1), ra, 20, fs);
     
    AnlySig(r) = SignalAnalyser(Signals.Signals, ProcSigs.ProcessedSignals);
    AnlySig(r) = AnlySig(r).SpecExtract(0.5, 1);
    AnlySig(r) = AnlySig(r).PrelimAnly(1);
    
    st_f = 0.9*AnlySig.Min_f;
    en_f = 1.2*AnlySig.Max_f;
    
    i = 2;
    
    for snr = 1:length(SNR)
        for ty = 1:length(ch_ty)
            for ch_le = 1:length(ch_len)
                
                type = ch_ty(ty);
                sn = SNR(snr);
                
                ch_l = ch_len(ch_le);
                ch_s = ch_spc(ch_le);
                
                Signals(r) = Signals(r).SigGen(type, ch_l, ch_s, st_f, en_f, 10, fs);
                
                ProcSigs(r) = ProcSigs(r).SigProc(Signals.SigGet(i), ra, sn, fs);
                AnlySig(r) = AnlySig(r).Update(Signals(r).Signals, ProcSigs(r).ProcessedSignals);
                AnlySig(r) = AnlySig(r).SpecExtract(0.5, i);
                AnlySig(r) = AnlySig(r).LFOTrack(i, 5);
                
                lfo =  AnlySig(r).Measured_LFOs{i-1, 'Measured_LFO'}{1,1};
                tAx =  AnlySig(r).Measured_LFOs{i-1, 'LFO_time_axis'}{1,1};
                plot(tAx{1,1}, lfo{1,1});
                
                LFOFits(r) = LFOFitter();
                    
                i = i + 1;
            end
        end
    end
end














