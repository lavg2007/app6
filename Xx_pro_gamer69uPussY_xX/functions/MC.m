function [coeff,error] = MC(x,y,num)
if num == 2
    A  = [numel(x) sum(x); sum(x) sum(x.^2)];
    B = [sum(y);sum(y.*x)];
    coeff = inv(A) * B;
    g = coeff(1) + coeff(2) * x;
end
   
% calcul des erreurs
% R2
error.R2 = sum((g-mean(y)).^2)/sum((y-mean(y)).^2);
% RMS
error.rms =sqrt(1/numel(x)*sum((g-y).^2));
end