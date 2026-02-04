function c=Aiyagari1994_ConsumptionFn(aprime, a, z,alpha,delta,r)
% Note that this is essentially just a copy-paste of the ReturnFn up to the
% point where we calculate consumption

w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));
c=w*z+(1+r)*a-aprime; % Budget Constraint
% c=w l_t+(1+r)a_t-a_{t+1}

end