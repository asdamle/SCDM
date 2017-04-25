function [] = write_amn(Q,nk,fname)

[nband, nwan]= size(Q{1});

disp('Writing amn file...');
fid   = fopen(fname, 'w');
fprintf(fid,'Input from SCDM-k\n');
fprintf(fid,'  %10d  %10d  %10d\n', nband, nk, nwan);
for k= 1 : nk
  for j = 1 : nwan % check
    for i = 1 : nband % check  
      fprintf(fid, '%5d%5d%5d  %16.12f  %16.12f\n', ...
        i, j, k, real(Q{k}(i,j)), imag(Q{k}(i,j)));
    end
  end
end
fclose(fid);