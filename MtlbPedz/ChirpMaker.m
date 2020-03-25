function [x, chrpLen] = ChirpMaker(chirps, spacing, M, a) 

    
 
    delhf = (1 - a^2)/(1 + 2*a*cos(pi) + a^2);
    dellf = (1 - a^2)/(1 + 2*a*cos(0) + a^2);
    grpdel = max(delhf,dellf);
    
    x = 0;
    for ch = 0:chirps-1
       x(ch*spacing + 1) = 1; 
    end
    
    chrpLen = round(M*grpdel);
    
    x(end+1:end+(chrpLen + 50)) = 0;
    
    for n = 1:M
        x = filter([a, 1], [1, a], x);
    end

end