


for data = 1:3
    switch data
        case 1
            dir = 'train';
            speeds = [10:10:100];
        case 2
            dir = 'val';
            speeds = [10:10:100];
        case 3
            dir = 'test';
            speeds = [10:10:100];
    end

[input, fs] = audioread(strcat('../Dataset/Data/', dir,'/input.wav'));

spl_num = length(speeds);
spl_len = floor(length(input)/spl_num);

input = input(1:spl_num*spl_len);

output = zeros(spl_num*spl_len,1);
LFO_cond = zeros(spl_num*spl_len,1);

for split = 0:spl_num-1
    in = input((split * spl_len) + 1:(split+1) * spl_len);
    speed = speeds(split+1);
    
    %Generate a single period of LFO
    LFO_single = generate_LFO(speed, fs, 0);

    %Now we place multiple of the single period LFO-signals after one another
    %to generate long LFO-signal that can be used to model the phaser.
    LFO = LFO_single;

    while length(LFO(:,1)) < length(in)*1.4  %We want to make sure it's long enough
        LFO = [LFO; LFO_single]; 
    end

    signal_phasered = phasing_algorithm(in, fs, 0.5, 0.4, LFO, 1);
    
    LFO_cond((split * spl_len) + 1:(split+1) * spl_len) = LFO(1:length(in),3);
    output((split * spl_len) + 1:(split+1) * spl_len) = signal_phasered(1:length(in));
    
    
end

audiowrite(strcat('../Dataset/Data/', dir,'/SwTo-input.wav'), [input, LFO_cond], fs)
audiowrite(strcat('../Dataset/Data/', dir,'/SwTo-target.wav'), output, fs)
    
    
    
        
end



