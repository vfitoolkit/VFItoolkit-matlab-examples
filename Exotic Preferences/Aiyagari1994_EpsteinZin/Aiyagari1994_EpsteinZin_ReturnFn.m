function F=Aiyagari1994_EpsteinZin_ReturnFn(aprime_val, a_val, s_val,alpha,delta,gamma,r)
% A list of the parameters to be used
% alpha
% delta
% gamma
% r


F=-Inf;
w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));
c=w*s_val+(1+r)*a_val-aprime_val; 
%c=wl_t+(1+r)a_t-a_{t+1}
if c>0
%     if mu==1
%         F=log(c);
%     else
%         F=(c^(1-mu) -1)/(1-mu);
%     end
    F=c; % As Epstein-Zin
end

end