function c=Aiyagari1994_ConsumptionFn(aprime, a, z,alpha,delta,r)
% The return function is essentially the combination of the utility
% function and the constraints.

F=-Inf;
w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));
c=w*z+(1+r)*a-aprime; % Budget Constraint
% c=w l_t+(1+r)a_t-a_{t+1}

end