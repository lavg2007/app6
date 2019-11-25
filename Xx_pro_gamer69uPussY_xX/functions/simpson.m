function [simpson_result,simpson_error] = simpson(x,f)
h = mean(diff(x));
simpson_result = h/3 * (f(1) + f(end) + 4*sum(f(2:2:end-1))+2*sum(f(3:2:end-1)));


% error
if numel(x) > 3
    fder3 = diff(diff(diff(f)./diff(x)));
    simpson_error = h^4/180 *(fder3(end)-fder3(1));
else
    simpson_error = NaN;
end
end