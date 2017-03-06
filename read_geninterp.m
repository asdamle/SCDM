function [kpts, evals] = read_geninterp(fname)

fid = fopen(fname,'r');
A = textscan(fid,'%d%f%f%f%f','Headerlines',3);

num_eig = length(find(A{1}==1));
num_kpt = max(A{1});
evals = reshape(A{5},num_eig,num_kpt);
temp = [A{2} A{3} A{4}];
kpts = temp(1:num_eig:end,:);
fclose(fid);