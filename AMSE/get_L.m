function L = get_L(p)
% This function generates the matrix L such that vech(A) = L*vec(A)
% for every symmetric p-dimensional matrix A.

% Initialize the matrix L
L = zeros(p*(p+1)/2, p^2);

e=@(a) [zeros(a-1,1);1;zeros(p-a,1)];

for j = 1:p
    for i = j:p
        E = e(i)*e(j)';
        u = zeros(p*(p+1)/2,1);
        u((j-1)*p + i - j*(j-1)/2) = 1;
        L = L + u*E(:)';
    end
end

end