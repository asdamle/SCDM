function [] = isolated_wannier(fin,fout)

% nwan must be integer, no threshold option

% Read basic information
fid=fopen(fin,'r');
nband=fscanf(fid,'%d',1);
nk=fscanf(fid,'%d',1);
nwan=fscanf(fid,'%d',1);
kpts = zeros(nk,3);
for k = 1 : nk
  kpts(k,:) = fscanf(fid,'%g',3);
end
fclose(fid);

[check, gamma_idx] = min(sum(abs(kpts),2));
if check ~= 0
    disp('no gamma point found')
end

[Psi,nx,ny,nz] = readUNK(gamma_idx);

disp('Computing columns...');
[~, ~, perm]=qr(Psi','vector');  
cols = perm(1:nwan);

[Xs,Ys,Zs] = ndgrid(0:nx-1,0:ny-1,0:nz-1);
Xs = Xs/nx;
Ys = Ys/ny;
Zs = Zs/nz;
rposFrac = [Xs(:) , Ys(:) , Zs(:)];

% Project to Phi state from the gamma calculation for Amn.
% The phase factor is important for approximating integration in the
% global domain.
disp('Performing projections...');
Amn = cell(nk,1);
for k = 1 : nk
    
    centerFrac = rposFrac(cols,:);
    centerFrac = centerFrac - round(centerFrac);
    phase = exp(1i*2*pi*centerFrac*kpts(k,:)');

    Psi = readUNK(k);
    Psi = diag(phase)*Psi(cols,:);

    Psi = Psi';

    %     Amn{k} = Psi / (sqrtm(Psi'*Psi));
    [Umn, Smn, Vmn] = svd(Psi,0);   %#ok
    Amn{k} = Umn*(Vmn');
end

% Write the overlap matrix .amn that can be called by Wannier90
disp('Writing amn file...');
fid   = fopen(fout, 'w');
fprintf(fid,'Input from SCDM-k\n');
fprintf(fid,'  %10d  %10d  %10d\n', nband, nk, nwan);
for k= 1 : nk
  for j = 1 : nwan % check
    for i = 1 : nband % check  
      fprintf(fid, '%5d%5d%5d  %16.12f  %16.12f\n', ...
        i, j, k, real(Amn{k}(i,j)), imag(Amn{k}(i,j)));
    end
  end
end
fclose(fid);


