function [rep,coeff] = interpol(x,y)
N = numel(x)
P = zeros(N,N)
phi = []
g = [];
Y =y';
for k = 1:N % k = ordre max du polynombre
    for i = 1:k
        phi(:,i) = (x.^(i-1))';
    end
    A = inv(phi' * phi) * phi'*Y;
    for m = 1:k
        P(m,k) = A(m)
    end
    g(k,:) = polyval(flip(P(:,k)'),x);
end
rep = g(N,:)
coeff = flip(P(:,N)')
end