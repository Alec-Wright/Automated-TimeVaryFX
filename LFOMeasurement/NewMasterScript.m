

import SignalHolder
import SignalProcessor
import SignalAnalyser
addpath('../DigitalEffectsPedals')

PedalName = 'PhaserSweetTone';
Digi = 1;
Rate = 0.3;
SNR = 40;

fs = 44100;


Signals = SignalHolder();
Signals = Signals.SigGen(1, 28, 30, 20, 22050, 10, fs);


ProcSigs = SignalProcessor(PedalName, Digi, Rate);
ProcSigs = ProcSigs.SigProc(Signals.SigGet(1), Rate, SNR, fs);


AnlySig = SignalAnalyser(Signals.Signals, ProcSigs.ProcessedSignals);
AnlySig = AnlySig.SpecExtract(0.5, 1);


AnlySig = AnlySig.PrelimAnly(1);

st_f = 0.9*AnlySig.Min_f;
en_f = 1.2*AnlySig.Max_f;
Signals = Signals.SigGen(1, 28, 30, st_f, en_f, 10, fs);
ProcSigs = ProcSigs.SigProc(Signals.SigGet(2), Rate, SNR, fs);
AnlySig = AnlySig.Update(Signals.Signals, ProcSigs.ProcessedSignals);
AnlySig = AnlySig.SpecExtract(0.5, 2);

AnlySig = AnlySig.LFOTrack(2, 5);


