function x = Generate_Obs(GMMStruct,GMMStruct_cont,N,epsilon)

    alphabet = [0,1];
    prob = [1-epsilon,epsilon];
    U = randsrc(1,N,[alphabet; prob]).';

    x_int = [];

    for m = 1:GMMStruct_cont.Order

        x_curr = mvnrnd(GMMStruct_cont.Means(:,m).',GMMStruct_cont.Covars(:,:,m),GMMStruct_cont.Alpha(m)*N);
        x_int = [x_int;x_curr];

    end
    
    x = [];

    for m = 1:GMMStruct.Order

        x_curr = mvnrnd(GMMStruct.Means(:,m).',GMMStruct.Covars(:,:,m),GMMStruct.Alpha(m)*N);
        x = [x;x_curr];

    end
    
    x = (1-U).*x + U.*x_int;

end