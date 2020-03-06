
datas = {'train','val','test'};
speeds = [0:20:100];

[input, fs] = audioread(string(strcat('../Dataset/Data/', datas(1),'/input.wav')));
splits(1) = 0;

for each = 2:length(datas)
    
    splits(each) = length(input);
    [x, fs] = audioread(string(strcat('../Dataset/Data/', datas(each),'/input.wav')));
    input = [input;x];
    
end

splits(end+1) = length(input);

output = zeros(length(speeds)*length(input),1);
LFO_cond = zeros(length(speeds)*length(input),1);
in_real = zeros(length(speeds)*length(input),1);
in_len = length(input);

for each = 1:length(speeds)
    
    speed = speeds(each);
    
    %Generate a single period of LFO
    LFO_single = generate_LFO(speed, fs, 0);

    %Now we place multiple of the single period LFO-signals after one another
    %to generate long LFO-signal that can be used to model the phaser.
    %LFO = LFO_single;

    repeats = ceil(in_len/length(LFO_single)) + 10;
    LFO = repmat(LFO_single,[repeats,1]);

    signal_phasered = phasing_algorithm(input, fs, 0.5, 0.6, LFO, 2);
    
    LFO_cond(((each-1) * in_len) + 1:each * in_len) = LFO(1:in_len,3);
    output(((each-1) * in_len) + 1:each * in_len) = signal_phasered(1:in_len);
    in_real(((each-1) * in_len) + 1:each * in_len) = input;
    
    
end

for each = 1:length(datas)
    
    stsamp = splits(each) + 1;
    ensamp = splits(each+1);
    len = ensamp - stsamp + 1;
    
    in_seg = zeros(length(speeds)*len, 1);
    LFO_seg = zeros(length(speeds)*len, 1);
    out_seg = zeros(length(speeds)*len, 1);
    
    for sps = 1:length(speeds)
        
        in_seg(len*(sps-1) + 1:len*sps) = in_real(splits(end)*(sps-1) + stsamp: splits(end)*(sps-1) + ensamp);
        LFO_seg(len*(sps-1) + 1:len*sps) = LFO_cond(splits(end)*(sps-1) + stsamp: splits(end)*(sps-1) + ensamp);
        out_seg(len*(sps-1) + 1:len*sps) = output(splits(end)*(sps-1) + stsamp: splits(end)*(sps-1) + ensamp);
        
        
    end
    
    audiowrite(strcat('../Dataset/Data/', string(datas(each)),'/SwTo-input.wav'), [in_seg, LFO_seg], fs)
    audiowrite(strcat('../Dataset/Data/', string(datas(each)),'/SwTo-target.wav'), out_seg, fs)
end

    