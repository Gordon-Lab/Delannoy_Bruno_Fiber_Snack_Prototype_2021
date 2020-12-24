%Compute tensor singular values as a function of instantiating tensor
%components

function [out] = tensor_component_cutoff(data,tensor_rank);

counter = 0;
for i=1:tensor_rank;
    counter = counter + 1;
    M1 = cp_als(data,i);
    sigma = M1.lambda;
    sigma = sigma./sum(sigma); %Convert to variance
    out(i) = min(sigma);
end;

figure; plot(out); %This gives a plot of variance explained as a function of tensor singular value included
title('Variance explained vs tensor singular value included');


    