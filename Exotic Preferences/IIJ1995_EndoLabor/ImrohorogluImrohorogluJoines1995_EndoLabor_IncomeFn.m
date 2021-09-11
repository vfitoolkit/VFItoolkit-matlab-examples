function income=ImrohorogluImrohorogluJoines1995_EndoLabor_IncomeFn(l,kprime,k,z,w,h,zeta, epsilon_j,I_j,SSdivw, Tr_beq,workinglifeincome,g,agej) %i

earnings=z*w*h*l*epsilon_j*I_j;
u=zeta*w*h; % Unemployment benefits are at replacement rate phi;
unemploymentbenefit=(l==0)*u*I_j;
socialsecuritybenefit=w*SSdivw*(1-I_j);
if g>0 % So not using baseline model
    socialsecuritybenefit=w*SSdivw*(1-I_j)*workinglifeincome*(1/((1+g)^agej));
    % Note: the actual social security benefit is constant, but the model
    % has been detrended by (1+g)^t and it is for these reason it here
    % appears to decrease with age.
end
% q=earnings+unemploymentbenefit+socialsecuritybenefit; % This is the notation used by IIJ1995 and is just here to ease comprehension.

% income=r*k+earnings+unemploymentbenefit+socialsecuritybenefit+Tr_beq;

% From Figure 3 it became clear that above commented out formula was not what
% is being plotted as the concept plotted does not include capital income
% but does include the pension (can be seen from the shape of the income
% profile during retirement)    
income=earnings+unemploymentbenefit+socialsecuritybenefit+Tr_beq;

end