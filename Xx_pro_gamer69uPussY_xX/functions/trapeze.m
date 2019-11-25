function [trapeze_result,trapeze_error] = trapeze(x,f)
h = mean(diff(x));
trapeze_result = h(1)/2 * (f(1) + f(end) + 2 * sum(f(2:end-1)));

% error
fder1 = diff(f)./diff(x);
trapeze_error = h^2/12 * (fder1(end)-fder1(1));
end