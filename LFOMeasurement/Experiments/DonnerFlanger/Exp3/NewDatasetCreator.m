addpath(genpath('../../../..'))

FileLen = 30;
split_num = 6;
spl_len = FileLen/split_num;

dataDir = '../../../../DatasetCreation';
subDirs = {'train', 'val', 'test'};

dataNames = {{'Rock1', 'Rock2', 'Rock4', 'Metal1', 'Metal3', 'Metal4',...
    'Country1', 'Country2', 'Country3', 'Classical2', 'Classical3',...
    'Classical4'},...
    {'Rock3', 'Classical1'},...
    {'Metal2', 'Clean1'}};

saveDir = '../../../../Dataset/Data/';
saveSubDirs = {'train', 'val', 'test'};
saveName = 'newDataInput';

for m = 1:length(subDirs)
    for n = 1:length(dataNames{m})
        [file, fs] = audioread(strcat...
            (dataDir,'/', subDirs{m},'/', dataNames{m}{1,n}, '.wav'));
        
        dataNames{m}{2,n} = file(1:fs*FileLen);
    end
end

dataset = [];
for o = 1:split_num
    for m = 1:length(subDirs)
        for n = 1:length(dataNames{m})
            data = dataNames{m}{2,n};
            dataset = [dataset; data((o-1)*fs*spl_len + 1:o*fs*spl_len)];
        end
    end
end

audiowrite('InputDataset2.wav',dataset, fs);
