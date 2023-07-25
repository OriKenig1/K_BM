function gamma = calc_gamma(x,GMMStruct)
    
    M = GMMStruct.Order;

    Pr = zeros(length(x),M);

    for m = 1:M
    
        Pr(:,m) = GMMStruct.Alpha(m).*(mvnpdf(x,GMMStruct.Means(:,m).',GMMStruct.Covars(:,:,m)));
    
    end

    gamma = Pr./(sum(Pr,2)*ones(1,M) + 1e-20);


end