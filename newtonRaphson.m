
function [x,i] = newtonRaphson(x,func, dfunc, tol)
     i = 0;
     while 1 
        if abs(feval(func,x)) <= tol
           break
        else
            i = i+1;
            x = x - feval(func,x)/feval(dfunc,x);
        end
     end
end

