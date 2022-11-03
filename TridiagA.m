function [Q] = TridiagA(lambda,mu,N)
%TridiagM creates a tridiagonal matrix of a three-term vector recurrence
%relation.
%   Inputs-->lambda, mu are birth and death rates. Note
%             that these are given as Matlab function
%             handles that depend on "n".
%
%             N: Size of the stochastic matrix
%             This is a parameter we should choose wisely since it
%             represents number of differential equations in the system
%
%   Output--> Q: Tridiagonal matrix with dimensions N x N,
% ______________________________________________________________________

% initializing diagonals
d1 = zeros(N-1,1);
d2 = zeros(N,1);
d3 = zeros(N-1,1);

% element 1,1
d2(1) = -lambda(0);

% element N,N
d2(N) = -mu(N-1);

% bottom diagonal elements
for i=0:N-2
    d1(i+1) = lambda(i);
end

% top diagonal elements
for i=1:N-1
    d3(i) = mu(i);
end

% main diagonal elements
for i=2:N-1
    d2(i) = -lambda(i-1)-mu(i-1);
end

% putting the diagonals together
Q = diag(d1,-1) + diag(d2,0) + diag(d3,1);
end