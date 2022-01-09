function Phi_aprime=ImrohorogluImrohorogluJoines2003_PhiaprimeFn(l,kprime,bprime,k,ebar,b,z,r,alpha,delta,A,epsilon_j,agej, Jr, emax,ebar,kgridspacing,ebargridspacing,n_a1,n_a2,n_a3)

w=(1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1));

earnings=z*w*l*epsilon_j;

ebarprime=ebar+min(earnings,emax)/agej;

% Change kprime value to the index of the grid point. Following line assumes linear grid.
kprime_c=1+floor(kprime/kgridspacing)+(rem(kprime,kgridspacing)>kgridspacing/2);

% Retired people will stay at ebar, Pre-retirement move to eprime.
ebarprime=b*ebar+(1-b)*ebarprime;
ebarprime_c=1+floor(ebarprime/egridspacing)+(rem(ebarprime,egridspacing)>egridspacing/2);
% Stop ebarprime from going off top of grid
ebarprime_c=n_a2*(ebarprime_c>=n_a2)+ebarprime_c*(ebarprime_c<n_a2);

% Because bprime is just 0 or 1 changing from value to index just involves adding one.
bprime_c=bprime+1;

n_a=[n_a1,n_a2,n_a3]; % Because of how evaluating fns on gpu works I had to split this vector in three, can now put it back together.
Phi_aprime=sub2ind_homemade(n_a,[kprime_c,ebarprime_c,bprime_c]);


end