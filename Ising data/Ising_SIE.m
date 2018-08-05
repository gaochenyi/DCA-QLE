% DESCRIPTION
% ===
% Infer couplings of the Ising model by the small-interaction expansion.
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
% - R. A. Neher and B. I. Shraiman,
%   "Statistical genetics and evolution of quantitative traits",
%   Rev. Mod. Phys., 83, 4, 1283 (2011).
% 
% HISTORY
% ===
% - 2018-06-21  v1.0
%   - initial draft

function J = Ising_SIE(S,B,N)

% check
[sz1,sz2] = size(S);
if (sz2 ~= N || sz1 ~= B)
  error('`S` should be provided as rows are sequences.')
end

% small interaction expansion
m = mean(S,1);
C = (S.'*S)/B - m.'*m;

J = zeros(N);

for i = 1:N
  for j = i+1:N
    J(j,i) = Ising_SIE_atom(m(i),m(j),C(j,i));
  end
end

J = J + J.';

end

% Do inference by eqs. 20--22 in Rev. Mod. Phys., 83, 4, 1283.
% 
% The formula assumes pairwise couplings are small (small enough
% for expansion).
% 
% chi_i  = <s_i>
% chi_ij = <s_i s_j> - <s_i><s_j>
function J_ij = Ising_SIE_atom(chi_i, chi_j, chi_ij)
J_ij = chi_ij / ((1-chi_i*chi_i)*(1-chi_j*chi_j));
end
