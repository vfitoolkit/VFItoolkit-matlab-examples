function F=Imrohoroglu1989_ReturnFn(aprime_val, a_val, s_val,sigma,r_l,r_b) %theta,y

F=-Inf;

if a_val>=0 %If lending
    r=r_l;
else %If borrowing
    r=r_b;
end

C=a_val-aprime_val/(1+r)+s_val;
if C>0
    F=(C^(1-sigma))/(1-sigma);
end

end