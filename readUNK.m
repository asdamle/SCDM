function [Unk,nx,ny,nz] = readUNK(k)

fname = sprintf('UNK%05d.1',k);
fid = fopen(fname, 'r');
nx = fscanf(fid,'%d',1);
ny = fscanf(fid,'%d',1);
nz = fscanf(fid,'%d',1);
kidx = fscanf(fid,'%d',1);
nband = fscanf(fid,'%d',1);
Unk = fscanf(fid,'%g',2*nx*ny*nz*nband);
Unk = Unk(1:2:end) + 1i*Unk(2:2:end);
Unk = reshape(Unk,[nx*ny*nz,nband]);
for g = 1 : nband
    Unk(:,g) = Unk(:,g) / norm(Unk(:,g));
end
fclose(fid);