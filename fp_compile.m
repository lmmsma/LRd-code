%Use this program to compile all the fixed points into a single matrix

%Reads in files of located fixed point
folder = 'fp_found/';%Need fixed points to be saved in this folder
bcls_tested = [1000:-50:600 590:-10:510 505 500:-10:70];%Vector of bcls
L = length(bcls_tested);

%Compile all fixed points into a single matrix, fp_found
fp_found = zeros(17, L);
for i = 1:L
   fname = [folder 'lrddata_1cell_b' num2str(bcls_tested(i))];
   fp = load(fname, 'Y');
   fp_found(:, i) = fp.Y(:, 1);
end

%Save the matrix to a separate file in the folder specified above
eval(['save ' folder 'compiled_fp.mat fp_found']);

%Display full matrix if desired (uncomment)
%disp(fp_found);
