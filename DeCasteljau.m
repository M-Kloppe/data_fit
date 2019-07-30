%De Casteljau algorithm for d=2

%M. Kloppe, Juni 2019

%Paramter:
%b1,b2,b3 - barycentric coordinates (columns)
% c - row vector of B-coefficients of a given polynomial
function [value]=DeCasteljau(b1,b2,b3,c)

B=[b1,b2,b3];

c11=c(1,1:3);
c12=c(1,[2,4,5]);
c13=c(1,[3,5,6]);


c21=c11*B';
c22=c12*B';
c23=c13*B';

C=[c21',c22',c23'];
M=C*B';

value=diag(M);

end