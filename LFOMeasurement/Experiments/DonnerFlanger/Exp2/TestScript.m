

DataDir = '../../../../Dataset/';
DataLocs = ["train", "val", "test"];
SaveLoc = 'Data/DonnerFlanger';
method = 1;
ch_len = 10;
ch_spc = 20;
f_st = 20;
f_en = 3000;
T = 10;
seg_len = 40;
fs = 44100;
Settings = {'color', '12', 'range', '12'};

% Signals = SignalHolder();
% ProcSigs = SignalProcessor(PedalName, Digi, rate);
% ProcSigs.Settings = settings;
% ProcSigs.SaveLoc = SaveLoc;
% % Create Preliminary Measurement signal
% Signals = Signals.SigGen(1, 15, 20, 20, 3000, 10, fs);
% ProcSigs = ProcSigs.SigProc(Signals.SigGet(1), rate, 40, fs);
% 
% 
% ProcSignal = ProcSigs.SigLoad(strcat(SaveLoc, 'prelim'));
% 
% ProcSigs.ProcessedSignals{1,1} = {ProcSignal};
% AnlySig = SignalAnalyser(Signals.Signals, ProcSigs.ProcessedSignals);
% AnlySig = AnlySig.SpecExtract(0.1, 1);
% AnlySig = AnlySig.PrelimAnly(1);

init_f = AnlySig.Initial_f;
f_st = AnlySig.Min_f;
f_en = AnlySig.Max_f;

% [ChirpTrains, Boundaries] = DataSetInput...
%     (DataDir, DataLocs, SaveLoc, ...
%                         method, ch_len, ch_spc, f_st*0.5, f_en*1.25, T, seg_len, fs);
                    
% [ProcSigs, AnlySigs] = ...
%             DataSetCreator(ChirpTrains, 'Flanger', 12, Settings, f_st, f_en, init_f);
% full_lfo = [];
% full_tAx = [];
% for n = 1:size(AnlySigs.Measured_LFOs,1)
%     full_lfo = [full_lfo; AnlySigs.Measured_LFOs{n,1}{1,1}];
%     full_tAx = [full_tAx; AnlySigs.Measured_LFOs{n,2}{1,1}];
% end
% full_lfo_flipped = -full_lfo + max(full_lfo) + min(full_lfo);
% 
% LFOFit1 = LFOFitter();
% 
% lfoP = LFOFitter.SimpleSineFit(full_lfo_flipped, full_tAx, init_f);

% LFOFit1 = LFOFit1.RectSineGridSearch...
%                 (full_lfo_flipped, full_tAx, init_f, 1 ,1);
%             
% elf = LFOFit1.LFOLookup(1, full_tAx);
% 
% LFOFit2 = LFOFitter();

for n = 1:size(AnlySigs.Measured_LFOs,1) - 1
    max(AnlySigs.Measured_LFOs{n,1}{1,1})
    lfo1 = AnlySigs.Measured_LFOs{n,1}{1,1};
    lfo2 = AnlySigs.Measured_LFOs{n+1,1}{1,1};
    lfo_sig = [lfo1; lfo2];
    
    tAx1 = AnlySigs.Measured_LFOs{n,2}{1,1};
    tAx2 = AnlySigs.Measured_LFOs{n+1,2}{1,1};
    tAx_sig = [tAx1; tAx2];
    lfo_sig = -lfo_sig + max(lfo_sig) + min(lfo_sig);
    lfo1 = lfo_sig(tAx_sig <= tAx1(end));
    lfo2 = lfo_sig(tAx_sig >= tAx2(1));
%     LFOFit2 = LFOFit2.RectSineGridSearch...
%                 (lfo_sig, tAx_sig, init_f, 1 ,2);
%     elf = LFOFit2.LFOLookup(n, tAx_sig);
    lfoP = LFOFitter.SimpleSineFit(lfo_sig, tAx_sig, init_f);
    lfoP1 = lfoP(tAx_sig <= tAx1(end));
    lfoP2 = lfoP(tAx_sig >= tAx2(1));
    error1 = mean(lfoP1 - lfo1).^2;
    error2 = mean(lfoP2 - lfo2).^2;
    maxe1 = max(abs(lfoP1 - lfo1));
    maxe2 = max(abs(lfoP2 - lfo2));
    error1 - error2
    maxe1 - maxe2
    
    plot(tAx_sig, lfo_sig);
    hold on
    plot(tAx_sig, lfoP);
    hold off
end
     