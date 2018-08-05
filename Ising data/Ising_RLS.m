% DESCRIPTION
% ===
% Infer couplings of the Ising model by the regularized least squares.
% 
% INPUT
% ===
% S       rows are sequences (possible states: -1 or 1)
% B       number of sequences
% N       number of loci
% lambda  strength of regularization
% 
% OUTPUT
% ===
% J       the matrix corresponding to couplings (symmetric and the diagonal
%         elements are zero)
% 
% REFERENCE
% ===
% - M. Andreatta and S. Laplagne and S. C. Li and S. Smale, 
%   "Prediction of residue-residue contacts from protein families using similarity kernels and least squares regularization",
%   arXiv:1311.1301v3 [q-bio.BM] (2014).
% 
% HISTORY
% ===
% - 2018-08-05
%   - copied from `script_fabio_data_time_average.m` (with additional comments)

function [J,invertible] = Ising_RLS(S,B,N,lambda)

% check
[sz1,sz2] = size(S);
if (sz2 ~= N || sz1 ~= B)
  error('`S` should be provided as rows are sequences.')
end

if numel(lambda) ~= 1
  error('This implementation only accepts scalar lambda.')
end

% regularized least squares
m = sum(S,1)/B; % faster than mean
C = (S.'*S)/B - m.'*m;
den = C*C;
if lambda > 0
  for k = 1:N
    den(k,k) = den(k,k) + lambda;
  end
end

r = rank(den);

if r >= N
  J = -(C/den);
  for k = 1:N
    J(k,k) = 0;
  end
else
  invertible = false;
  J = 0;
  fprintf('The denominator in the RLS formula is not invertible.\n')
end

end
