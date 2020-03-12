function InvChirp = ChirpInverter(x, M, a) 

    InvChirp = flip(x);

    for n = 1:M
        InvChirp = filter([a, 1], [1, a], InvChirp);
    end

end