addpath(genpath('Dataset'))
addpath(genpath('NNetz'))

Pedal = 'DonnerFlanger';
Code = 'K';

SaveDir = 'Dataset/Data/';
SaveLocs = ["train", "val", "test"];

Files = dir(strcat('NNetz/Results/',Code,'*'));

table_header = [["Job_Name", "string"];
                ["Data_Name", "string"];
            ["Hidden_Size", "int16"]; ...
            ["Range", "int16"]; ...
            ["Single_Num", "int16"]; ...
            ["Best_v_loss", "double"]; ...
            ["Best_t_loss", "double"]];

% Make table using fieldnames & value types from above
Results = table('Size',[0,size(table_header,1)],... 
    'VariableNames', table_header(:,1),...
    'VariableTypes', table_header(:,2));

for n = 1:length(Files)
    
    job_name = Files(n).name;
    data_name = Files(n).name(1:end-3);
    
    val = jsondecode(fileread(strcat('NNetz/Results/',job_name,'/config.json')));
    tloss = fileread(strcat('NNetz/Results/',job_name,'/tloss.txt'));
    tloss = str2num(tloss(2:end-1));
    
    Key = 'rg';
    Index = strfind(data_name,Key);
    Range = sscanf(data_name(Index(1) + length(Key):end), '%g', 1);
    
    Key = 'Singles';
    Index = strfind(data_name,Key);
    Sing_Num = sscanf(data_name(Index(1) + length(Key):end), '%g', 1);
    
    new_row = {job_name, data_name, val.hidden_size, Range, Sing_Num, min(val.vloss_list), tloss};
    
    Results = [Results; new_row];
    
end

winners = {};
ranges = unique(Results.Range);
for n = 1:length(ranges)
   subTab = Results(Results.Range == ranges(n), :);
   singles = unique(subTab.Single_Num);
   
   score = {};
   for m = 1:length(singles)
       score(m, :) = {singles(m), mean(subTab{subTab.Single_Num == singles(m), 6})};
   end
   [~,i] = min([score{:,2}]);
   
   wins = subTab{subTab.Single_Num == score{i,1}, :};
   winners{end+1, 1} = wins{1,2};
end

for subDir = 1:length(SaveLocs)
    loc = strcat(SaveDir, SaveLocs(subDir), '/');
    for n = 1:length(winners)
        load_name1 = winners{n,1};
        Key = 'rg';
        Index = strfind(load_name1,Key);
        Range1 = sscanf(load_name1(Index(1) + length(Key):end), '%g', 1);
        
        [Input1, fs] = audioread(strcat(loc, winners{n,1}, '-input.wav'));
        [Target1, fs] = audioread(strcat(loc, winners{n,1}, '-target.wav'));
        
        for m = n+1:length(winners)
            load_name2 = winners{m,1};
            Key = 'rg';
            Index = strfind(load_name2,Key);
            Range2 = sscanf(load_name2(Index(1) + length(Key):end), '%g', 1);
            
            [Input2, fs] = audioread(strcat(loc, winners{m,1}, '-input.wav'));
            [Target2, fs] = audioread(strcat(loc, winners{m,1}, '-target.wav'));
            
            audiowrite(strcat(loc,Code,Pedal, 'Doubles_rg', num2str(Range1), '_rg', num2str(Range2), '-input.wav'), [Input1; Input2], fs);
            audiowrite(strcat(loc,Code,Pedal, 'Doubles_rg', num2str(Range1), '_rg', num2str(Range2), '-target.wav'), [Target1; Target2], fs);
        end
    end
    
    [Input1, fs] = audioread(strcat(loc, winners{1,1}, '-input.wav'));
    [Target1, fs] = audioread(strcat(loc, winners{1,1}, '-target.wav'));
    [Input2, fs] = audioread(strcat(loc, winners{2,1}, '-input.wav'));
    [Target2, fs] = audioread(strcat(loc, winners{2,1}, '-target.wav'));
    [Input3, fs] = audioread(strcat(loc, winners{3,1}, '-input.wav'));
    [Target3, fs] = audioread(strcat(loc, winners{3,1}, '-target.wav'));
    audiowrite(strcat(loc,Code,Pedal,'Triple', '-input.wav'), [Input1; Input2; Input3], fs);
    audiowrite(strcat(loc,Code,Pedal,'Triple', '-target.wav'), [Target1; Target2; Target3], fs);
end