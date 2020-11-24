addpath('../../MatlabClasses')
import SignalHolder
import SignalProcessor
import SignalAnalyser
import LFOFitter

% Signal Processing Info
PedalName = 'DonnerFlanger';
Digi = 0;
rate = 9;
fs = 44100;

Signals = SignalHolder();
ProcSigs = SignalProcessor(PedalName, Digi, rate);

Signals = Signals.SigGen(1, 28, 30, 20, 22050, 10, fs);
ProcSigs = ProcSigs.SigProc(Signals.SigGet(1), rate, 40, fs);
ProcSigs = ProcSigs.SigLoad(1);

AnlySig = SignalAnalyser(Signals.Signals, ProcSigs.ProcessedSignals);
AnlySig = AnlySig.SpecExtract(0.5, 1);
AnlySig = AnlySig.PrelimAnly(1);
st_f = 0.9*AnlySig.Min_f;
en_f = 1.1*AnlySig.Max_f;

Signals = Signals.SigGen(1, 19, 20, st_f, en_f, 10, fs);
ProcSigs = ProcSigs.SigProc(Signals.SigGet(2), rate, 40, fs);
ProcSigs = ProcSigs.SigLoad(2);
AnlySig = SignalAnalyser(Signals.Signals, ProcSigs.ProcessedSignals);
AnlySig = AnlySig.SpecExtract(0.5, 2);
AnlySig = AnlySig.LFOTrackPeak(2, 10);