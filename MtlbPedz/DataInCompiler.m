
clear all

chords = dir('../../FraunhoferDatasets/IDMT-SMT-GUITAR_V2/dataset1/*Chords');

for n = 1:length(chords)
    files = dir(strcat('../../FraunhoferDatasets/IDMT-SMT-GUITAR_V2/dataset1/',...
    chords(n).name,'/audio/*.wav'));

    for each = 1:length(files)
        if mod(each,6) == 0 
            [x,fs] = audioread(strcat(files(each).folder,'/',files(each).name));
            
            if exist('val', 'var') == 1
                val = [val;x];
            else
                val = x; 
            end
            
        elseif mod(each,7) == 0 
            [x,fs] = audioread(strcat(files(each).folder,'/',files(each).name));
            
            if exist('test', 'var') == 1
                test = [test;x];
            else
                test = x;
            end
        else
            [x,fs] = audioread(strcat(files(each).folder,'/',files(each).name));
            
            if exist('train', 'var') == 1
                train = [train;x];
            else
                train = x;
            end
            
        end
    end
 
end

for data = 1:3
    switch data
        case 1
            dset = 'train';
            data = train;
        case 2
            dset = 'val';
            data = val;
        case 3
            dset = 'test';
            data = test;
    end
    
    files = dir(strcat('../../FraunhoferDatasets/CutGuitarForTimeVaryFX/', dset,'/*.wav'));
    for each = 1:length(files)      
        [x,fs] = audioread(strcat(files(each).folder,'/',files(each).name));
        data = [data;x];
    end
   audiowrite(strcat('../../TimeVarFx/Dataset/Data/',dset,'/input.wav'),data,fs)
    
end

