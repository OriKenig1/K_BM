function V = calc_V(x,GMMStruct)

    [N,p] = size(x);

    total = 0;
    Order = GMMStruct.Order;
    
    W1 = GMMStruct.Alpha;
    M1 = GMMStruct.Means;
    C1 = GMMStruct.Covars;

    for k = 1:1:Order
        for m = 1:1:Order
            total = total + W1(m)*W1(k)*mvnpdf(M1(:,m)', M1(:,k)', C1(:,:,m) + C1(:,:,k));
        end
    end

V = total - ( 2 * sum(GMMStruct.f)./N ); 