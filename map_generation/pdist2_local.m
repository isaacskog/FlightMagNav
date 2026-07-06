function D = pdist2_local(A,B)
%PDIST2_LOCAL Pairwise Euclidean distances without Statistics Toolbox.

    AA = sum(A.^2,2);
    BB = sum(B.^2,2).';

    D2 = AA + BB - 2*(A*B.');
    D2 = max(D2,0);

    D = sqrt(D2);
end
