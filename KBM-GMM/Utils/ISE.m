function total = ISE(GMM1, GMM2)

    total = 0;
    GMM1_squared = 0;
    GMM1_GMM2 = 0;
    GMM2_squared = 0;

    Order1 = GMM1.Order;
    Order2 = GMM2.Order;
    
    W1 = GMM1.Alpha;
    M1 = GMM1.Means;
    C1 = GMM1.Covars;
    
    W2 = GMM2.Alpha;
    M2 = GMM2.Means;
    C2 = GMM2.Covars;
    
    for k = 1:1:Order1
        for m = 1:1:Order1
            GMM1_squared = GMM1_squared + W1(m)*W1(k)*mvnpdf(M1(:,m)', M1(:,k)', C1(:,:,m) + C1(:,:,k));         
        end
    end
    
    for k = 1:1:Order1
        for m = 1:1:Order2
            GMM1_GMM2 = GMM1_GMM2 + W1(k)*W2(m)*mvnpdf(M2(:,m)', M1(:,k)', C2(:,:,m) + C1(:,:,k));
        end
    end

    for k = 1:1:Order2
        for m = 1:1:Order2
            GMM2_squared = GMM2_squared + W2(m)*W2(k)*mvnpdf(M2(:,m)', M2(:,k)', C2(:,:,m) + C2(:,:,k));         
        end
    end

    total = GMM1_squared - 2 * GMM1_GMM2  + GMM2_squared;

end