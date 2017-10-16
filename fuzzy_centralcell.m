function [centralcell_fuzzycomp]=fuzzy_centralcell(a,b,i,j)
% Calculates the fuzzy similarity as the maximum grade of membership to the
% intersection set of two fuzzy sets a and b.

centralcell_fuzzycomp=max(min(a(i,j,1),b(i,j,1)),min(a(i,j,2),b(i,j,2)));

end