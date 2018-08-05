% DESCRIPTION
% ===
% Infer couplings of the Ising model by the naive mean-field inversion.
% 
% INPUT
% ===
% S   rows are sequences (possible states: -1 or 1)
% B   number of sequences
% N   number of loci
% 
% OUTPUT
% ===
% J   the matrix corresponding to couplings (symmetric and the diagonal
%     elements are zero)
% 
% REFERENCE
% ===
% - H. J. Kappen and F. B. Rodriguez,
%   "Efficient learning in Boltzmann machines using linear response theory",
%   Neural Comput., 10, 1137 (1998).
% - T. Tanaka,
%   "Mean-field theory of Boltzmann machine learning",
%   Phys. Rev. E, 58, 2302 (1998).
% 
% HISTORY
% ===
% - 2018-05-05  v1.1
%   - add flag for invertibility
% - 2018-04-15  v1.0
%   - initial draft

function [J,invertible] = Ising_NMFI(S,B,N)

% check
[sz1,sz2] = size(S);
if (sz2 ~= N || sz1 ~= B)
  error('`S` should be provided as rows are sequences.')
end

% naive mean-field inversion
m = sum(S,1)/B; % faster than mean
C = (S.'*S)/B - m.'*m;

r = rank(C);

if r >= N
  invertible = true;
  J = -inv(C);
  for k = 1:N
    J(k,k) = 0;
  end
else
  invertible = false;
  J = 0;
  fprintf('The covariance matrix is not invertible.\n')
end

end
