function q = NewRaph(phi,jacfunc,q,t,kmax,tol)
% Calculates roots to given function (func) from the jacobian function, and
% from given initial guess values, maximum iteration and termination
% tolerance.
%
% by Mads Rosenh?j, Aarhus University, 2019.
%
% Input:
%   func    = Root search function
%   jacfunc = Jacobian function
%   q       = initial guess
%   K       = Constant structure
%   kmax    = maximum number of iterations
%   tol     = termination tolerance
%
% Output:
%   x       = calculated step size
%   delta2  = improved initial step size, to be used next time LineSearch
%             is called with a similar function as input.



n = 1;
normCorr = 100;
while normCorr > tol && n < kmax
    res = phi(q,t);
    jac = jacfunc(q,t);
    Corr = jac\res;
    q = q - Corr;
    normCorr = norm(Corr);
    n = n + 1;
end



