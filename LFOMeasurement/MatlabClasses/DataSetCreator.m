function [ProcSigs, AnlySigs] = ...
                    DataSetCreator(Signals, PedalName, Rate, Settings, f_st, f_en, init_f)
   
    ProcSigs = SignalProcessor(PedalName, 0, Rate);
    ProcSigs.Settings = Settings;
    
    ProcSignal = SignalProcessor.SigLoad(Signals.SaveLoc);
    
    AnlySig = SignalAnalyser(Signals.Signals, ProcSigs.ProcessedSignals);
    AnlySig = AnlySig.BatchSpecExtract(0.5, ProcSignal);
    
    AnlySig.Min_f = f_st;
    AnlySig.Max_f = f_en;
    AnlySig.Initial_f = init_f;
    
    % Retrieve Chirp Trains from processed audio
    for n = 1:size(Signals.Signals,1)
      
        AnlySig = AnlySig.LFOTrack(n, 5, 2);
        
    end
    
end