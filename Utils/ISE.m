function total = ISE(GMM1, GMM2)
    total = 0;
    Order = GMM1.Order;
    
    W1 = GMM1.Alpha;
    M1 = GMM1.Means;
    C1 = GMM1.Covars;
    
    W2 = GMM2.Alpha;
    M2 = GMM2.Means;
    C2 = GMM2.Covars;

    for k = 1:1:Order
        for m = 1:1:Order
            total = total + W1(m)*W1(k)*mvnpdf(M1(:,m)', M1(:,k)', C1(:,:,m) + C1(:,:,k));
            total = total - 2 * W1(k)*W2(m)*mvnpdf(M2(:,m)', M1(:,k)', C2(:,:,m) + C1(:,:,k));
            total = total + W2(k)*W2(m)*mvnpdf(M2(:,m)', M2(:,k)', C2(:,:,m) + C2(:,:,k));
        end
    end

    if(isnan(total))
        total = 10^13;
    end

end