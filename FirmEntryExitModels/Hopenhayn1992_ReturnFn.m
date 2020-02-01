function F=Hopenhayn1992_ReturnFn(n_val,aprime_val, a_val, s_val, p, alpha, cf)
% Note that neither aprime_val nor a_val is actually used for anything in
% this model. But VFI toolkit is not set up to handle that they do not exist.

F=-Inf;

% This example is taken from Chris Edmonds lecture notes.
pi=p*s_val*(n_val^alpha)-n_val-cf;

F=pi;

end