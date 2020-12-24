%Create 'randomized' tensor

function out = rand_core(data);

%Input: data is an m x n x p tensor
%Output: out is the first value of the core tensor as calculated from
%shuffling the 'data' input

%Note if doing this on the Mirpur cohort, the rows of the matrix need to be
%the patients and the columns need to be taxa

for i=1:length(data(1,1,:)); %Loop through the third dimension
    data_new = data(:,:,i);
    data_new = data_new(randperm(end),randperm(end));
    data_rand(:,:,i) = data_new;
    clear data_new;
end;

    