function [spatialfuzzy_similarity]=nbh_cc(a,b,c,d,i,j)
% Calculate the neighbourhood+central cell similarity
% if a is neighbourhood fuzzy vector then b is the crisp vector of the
% corresponding comparison set, vice versa for c and d.

A(i,j)=max(min(a(i,j,1),b(i,j,1)),min(a(i,j,2),b(i,j,2)));
B(i,j)=max(min(c(i,j,1),d(i,j,1)),min(c(i,j,2),d(i,j,2)));
spatialfuzzy_similarity=min(A(i,j),B(i,j));
end