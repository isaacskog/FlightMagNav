function Phi = spatial_rbf_2d(r,centers,lengthScale)
% Isotropic Gaussian RBFs in the horizontal NED plane.

D2 = pdist2_local(r,centers).^2;
Phi = exp(-0.5*D2/lengthScale^2);
end


function D = pdist2_local(A,B)
% Pairwise Euclidean distances.

AA = sum(A.^2,2);
BB = sum(B.^2,2).';

D2 = AA + BB - 2*(A*B.');
D2 = max(D2,0);

D = sqrt(D2);
end