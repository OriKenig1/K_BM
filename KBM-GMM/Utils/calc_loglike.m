function diff = calc_loglike(x,N,GMMStruct,GMMStruct2)
   
    f1 = zeros(N,1);
    f2 = zeros(N,1);
    
    for m = 1:GMMStruct.Order

        f1 = f1 + GMMStruct.Alpha(m)* ...
            mvnpdf(x,GMMStruct.Means(:,m)',GMMStruct.Covars(:,:,m));

        f2 = f2 + GMMStruct2.Alpha(m)* ...
            mvnpdf(x,GMMStruct2.Means(:,m)',GMMStruct2.Covars(:,:,m));

    end

    loglike1 = sum(log(f1 + eps));
    loglike2 = sum(log(f2 + eps));

    diff = abs( ( loglike1 - loglike2 ) / loglike2 );

end


