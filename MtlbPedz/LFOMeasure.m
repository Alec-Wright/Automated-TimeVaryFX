
speeds = [0:20:100];
fs = 44100;
chirps = 600;
spacing = 30;
M = 64;
a = -0.9;

testSig = ChirpMaker(chirps, spacing, M, a, fs);
testSig = testSig';

in_len = length(testSig);

output = zeros(length(speeds)*length(in_len),1);
LFO_cond = zeros(length(speeds)*length(in_len),1);
invChirp = zeros(length(speeds)*length(in_len),1);

for each = 1:length(speeds)
    
    speed = speeds(each);
    
    %Generate a single period of LFO
    LFO_single = generate_LFO(speed, fs, 0);
    
    repeats = ceil(in_len/length(LFO_single)) + 10;
    LFO = repmat(LFO_single,[repeats,1]);
    
    signal_phasered = phasing_algorithm(testSig, fs, 0.5, 0.6, LFO, 2);
    
    inv = ChirpInverter(signal_phasered, M, a);
    inv = flip(inv);
    
    LFO_cond(((each-1) * in_len) + 1:each * in_len) = LFO(1:in_len,3);
    output(((each-1) * in_len) + 1:each * in_len) = signal_phasered(1:in_len);
    invChirp(((each-1) * in_len) + 1:each * in_len) = inv(1:in_len);
 
end

for each 
