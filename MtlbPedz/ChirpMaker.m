function x = ChirpMaker(chirps, spacing, M, a, fs) 

    sps = round(fs*spacing/1000);
 
    delhf = (1 - a^2)/(1 + 2*a*cos(pi) + a^2);
    dellf = (1 - a^2)/(1 + 2*a*cos(0) + a^2);
    grpdel = max(delhf,dellf);
    
    x = 0;
    for ch = 0:chirps-1
       x(ch*sps + 1) = 1; 
    end
    
    x(end+1:end+(round(M*grpdel)) + 50) = 0;
    
    for n = 1:M
        x = filter([a, 1], [1, a], x);
    end

end