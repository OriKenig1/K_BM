function GMMStruct_0 = KMEANS_initGuess(x,M)

[N,p] = size(x);

GMMStruct_0 = CreateGMMStruct(p, M);

[idx,C] = kmeans(x,M,'MaxIter',200);

for m = 1:M

    GMMStruct_0.Alpha(m) = sum(idx == m) / N;

    GMMStruct_0.Means(:,m) = C(m,:)';

    GMMStruct_0.Covars(:,:,m) = eye(p);

end

GMMStruct_0.Alpha = GMMStruct_0.Alpha;