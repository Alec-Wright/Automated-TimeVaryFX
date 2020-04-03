clear all
close all

%% 1 - User Parameters

% Create plots (1 for yes, anything else for no)
plots = 0;

% Speeds to test at
speeds = [0:20:100];
% speeds = [100];
fs = 44100;
% Number of chirps to include in the chirp train
chirps = 1200;
% Spacing between chrip starts (in ms)
spacing = 30;
% Number of all-pass filters and coefficient
M = 64;
a = -0.9;

%% 2 - Create chirp train

% Convert from ms to samples
spacing = round(fs*spacing/1000);

% Create chirp train using ChirpMaker function
[testSig, chrpLen] = ChirpMaker(chirps, spacing, M, a);
testSig = testSig';

% Lazy fix so fft is always even in length
if mod(chrpLen,2) == 1
    chrpLen = chrpLen + 1;
end

% Save length of chirp train
in_len = length(testSig);

%output = zeros(length(speeds)*length(in_len),1);
%LFOCond = zeros(length(speeds)*length(in_len),1);
%invChirp = zeros(length(speeds)*length(in_len),1);

if plots == 1
    figure
end

%% Process chirp train for each speed setting

for each = 1:length(speeds)
    
    % Save speed setting and true frequency
    speed = speeds(each);
    speedknob(each).speed = speed;
    speedknob(each).trueFreq = 0.069*exp(0.040*speed);
    
    %Generate a single period of LFO
    LFO_single = generate_LFO(speed, fs, 0);
    %Generate a full length LFO
    repeats = ceil(in_len/length(LFO_single)) + 10;
    LFO = repmat(LFO_single,[repeats,1]);
    
    % Save true single LFO and repeated LFO
    speedknob(each).LFO_single = LFO_single;
    speedknob(each).LFO = LFO;
    
    %Apply phaser to test signal
    signal_phasered = phasing_algorithm(testSig, fs, 0.5, 0.6, LFO, 2);
    %Invert processed chirps to impulse with reverse all-pass chain
    inv = flip(ChirpInverter(signal_phasered, M, a));
    
    % Save phasered signal and response chirps
    speedknob(each).signal_phasered = signal_phasered;
    speedknob(each).imp_resps = inv;
    
    % Create matrix to hold chirps, in frequency domain and time-domain
    speedknob(each).chirpogramMag = zeros(round(chrpLen/2), chirps);
    speedknob(each).chirpogramPhi = zeros(round(chrpLen/2), chirps);
    speedknob(each).chirpstd = zeros(chrpLen, chirps);
    
    % Iterate over all chirps
    for n = 0:chirps-1
        startInd = (each-1)*in_len + 1 + n*spacing;
        endInd = (each-1)*in_len + n*spacing + chrpLen;
        
        % Save time domain chirp
        speedknob(each).chirpstd(:, n+1) = inv(startInd:endInd);
        
        % Save frequency domain chirp
        fd = fft(speedknob(each).chirpstd(:, n+1));
        speedknob(each).chirpogramMag(:, n+1) = abs(fd(1:chrpLen/2));
        speedknob(each).chirpogramPhi(:, n+1) = angle(fd(1:chrpLen/2));
    end
    
    peak_picker(speedknob(each), length(testSig), spacing)
    
end
    
    
    

    
    speedknob(each).chirpFd = fft(speedknob(each).chirpogramMag);

%      B = autocorr2d(speedknob(each).chirpogramMag);
    
%      speedknob(each).chirpCorr = B;
     
     corrLen = length(xcorr(speedknob(each).chirpogramMag(1,:)));
     
     speedknob(each).chirpCorr = zeros(corrLen, chrpLen/2);
     
     
     for bins = 1:chrpLen/2
         speedknob(each).chirpCorr(:,bins) = xcorr(speedknob(each).chirpogramMag(bins,:));
     end
     
     speedknob(each).autoCorr = mean(speedknob(each).chirpCorr,2);
     
     %subplot(3,1,1)
%      surf(B);
%      shading interp;
    if plots == 1
        figure(1)
        subplot(3,4,each*2 - 1)
        imagesc(speedknob(each).chirpogramMag)
    
     
       subplot(3,4,each*2)
       plot(speedknob(each).autoCorr)
%      findpeaks(speedknob(each).autoCorr)
    end
     
     [pks,locs,w,p] = findpeaks(speedknob(each).autoCorr);
     
     prompeaks = pks(p>0.5*max(p));
     promlocs = locs(p>0.5*max(p));
     
     if plots == 1
        figure(2)
        subplot(3,2,each)
        plot(xcorr(output))
     end
     
%      figure(2)
%      subplot(3,2,each)
%      plot(abs(fft(tingy)))
    
%      [vals, inds] = maxk(B,5);
%      [vals2, inds2] = maxk(vals',5);
%  
%end

% % create the video writer with 1 fps
% writerObj = VideoWriter('myVideo.avi');
% writerObj.FrameRate = 10;
% % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% figure(3)
% for i=1:chirps
%     % convert the image to a frame
%     plot(speedknob(each).chirpogramMag(:,i))
%     frame = getframe(3);  
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);


function [predLFO] = peak_picker( x, LFOLen, spacing)

    TargLFO = x.LFO(:,3);
    TargLFO = (TargLFO - min(TargLFO))/(max(TargLFO)- min(TargLFO));

    magSpec = x.chirpogramMag(:,1);
    [pks, locs, wids, proms] = findpeaks(-magSpec);
    peaks_tot = length(pks);
    
    peak_num = wids>0.1*max(wids);
    peak_num = find(peak_num);
    
    tot_ch = size(x.chirpogramMag,2) - 1;
    
    for i = 1:length(peak_num)
        peak(i).bin = zeros(tot_ch,4);
    end

    for chirps = 1:tot_ch
        
        magSpec = x.chirpogramMag(:,chirps);
        [pks, locs, wids, proms] = findpeaks(-magSpec);
        
        %plot(magSpec)
        if length(pks) > peaks_tot
            disp('da_fuq')
        end
        
        for i = 1:length(peak_num)
            peak(i).bin(chirps,1) = locs(peak_num(i));
            peak(i).bin(chirps,2) = pks(peak_num(i));
            peak(i).bin(chirps,3) = wids(peak_num(i));
            peak(i).bin(chirps,4) = proms(peak_num(i));
        end

    end
    
    predLFO = zeros(LFOLen,length(peak_num));
    figure
    
    for i = 1:length(peak_num)
        
        peak(i).max = max(peak(i).bin(:,1));
        peak(i).min = min(peak(i).bin(:,1));
        
        peak(i).bin(:,1) = (peak(i).bin(:,1) - peak(i).min)...
              /(peak(i).max - peak(i).min);
        
        xAx = 1:spacing:spacing*length(peak(i).bin);
        
        peak(i).predLFOcuSp = spline(xAx, peak(i).bin(:,1), 1:LFOLen);
        %peak(i).predLFOcuSp = peak(i).predLFOcuSp/max(peak(i).predLFOcuSp);
        
        t = 0:0.03:(tot_ch - 1)*0.03;
        LFOx = peak(i).bin(:,1);
        
        ts = 0:1/44100:(tot_ch - 1)*0.03;
        
        [Ts,T] = ndgrid(ts,t);
        peak(i).predLFOSinc = sinc(Ts - T)*LFOx;
        
        peak(i).predLFOSinc = (peak(i).predLFOSinc - min(peak(i).predLFOSinc))...
            /(max(peak(i).predLFOSinc) - min(peak(i).predLFOSinc));
        
        subplot(2,2,i)
        
        plot(TargLFO(1:length(peak(i).predLFOcuSp)))
        hold on
        plot(peak(i).predLFOcuSp)
        plot(peak(i).predLFOSinc)
        
        peak(i).error = TargLFO(1:length(peak(i).predLFOcuSp)) - peak(i).predLFOcuSp';
        peak(i).errorSinc = TargLFO(1:length(peak(i).predLFOSinc)) - peak(i).predLFOSinc;
        
        peak(i).meanErr = mean(abs(peak(i).error));
        peak(i).meanErrSinc = mean(abs(peak(i).errorSinc));
        peak(i).meanHeight = mean(abs(peak(i).bin(:,2)));
        peak(i).meanWid = mean(abs(peak(i).bin(:,3)));
        peak(i).meanProm = mean(abs(peak(i).bin(:,4)));
        peak(i).binDiff = peak(i).max - peak(i).min;
        
        plot(abs(peak(i).error))
        plot(abs(peak(i).errorSinc))
        legend('Target','Predicted Spline','Predicted Sinc', 'Error Spline', 'Error Sinc')
        title('PredLFO + error')
        
        
        
%         subplot(3,4,i*4 - 2)
%         plot(peak(i).bin(:,2))
%         title('height')
%         
%         subplot(3,4,i*4 - 1)
%         plot(peak(i).bin(:,3))
%         title('widths')
%         
%         subplot(3,4,i*4)
%         plot(peak(i).bin(:,4))
%         title('prominence')
        
    end
    
    subplot(2,2,4)
    scatter(1:3, [peak(i).meanErr])
    hold on
    scatter(1:3, [peak(i).meanHeight])
    
    
 
    
    
    
end


    
  
