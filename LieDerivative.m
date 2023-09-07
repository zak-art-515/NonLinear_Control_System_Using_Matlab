function df=LieDerivative(h,x)
% The LieDerivaive MATLAB function is used
% to find the jacobian vector of a given output
% h(x) : Is the output function
% x : The state vector
% df : The jacobian of h along x
    if nargin<2 & nargin==0
        error('not enough input argument');
    end
    df=[];
    n=length(x);
    for ii=1:n
        xx=x(ii);
        dff=diff(h,xx);
        df=[df,dff];
    end
    df;
end