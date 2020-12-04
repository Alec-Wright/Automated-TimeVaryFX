
PrelimSaveLoc = 'Data/EDonnerFlangerPrelim';
SaveLoc = 'Data/EDonnerFlanger';

DataDir = '../../../../Dataset/Data/';
DataLocs = ["train", "val", "test"];


PedalName = 'EDonnerFlanger';
method = 1;
ch_len = 9;
ch_spc = 15;
f_st = 20;
f_en = 3000;
T = 10;
seg_len = 40;
fs = 44100;
Settings = {'rate', '12' ; 'color', '12'; 'range', '12'};

Pedal(1) = PedalAnalyser(PedalName, Settings, fs, SaveLoc, ...
    PrelimSaveLoc, DataDir, DataLocs);
Pedal(1) = Pedal(1).PrelimAnalysisGen(1, ch_len, ch_spc, f_st, f_en, T);
Pedal(1) = Pedal(1).PrelimSigAnalayse();

Pedal(1) = Pedal(1).ConstructInputData(method, ch_len, ch_spc, T, seg_len);
Pedal(1) = Pedal(1).ReadTargetData();

Pedal(1) = Pedal(1).LFOFitting(1);

Pedal(1) = Pedal(1).DataMake(1);

Pedal(1) = Pedal(1).DataMakeSingles(1);

frcuts = [2,3,4,6];
for n = 1:length(frcuts)
    Pedal(1) = Pedal(1).DataMakeCutFr(1, frcuts(n));
end

maxerrcuts = [2, 5, 10, 15];
for n = 1:length(maxerrcuts)
    Pedal(1) = Pedal(1).DataMakeCutMxEr(1, maxerrcuts(n));
end

averrcuts = [200, 400, 800, 1600, 2000];
for n = 1:length(averrcuts)
    Pedal(1) = Pedal(1).DataMakeCutAvEr(1, averrcuts(n));
end

save('DatasetInfo.mat', 'Pedal')