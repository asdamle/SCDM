function [] = valance_and_conduction_wannier(fin,fout)

% if nwan >= 1 use as number of wannier functions to choose
% if nwan < 1 use as a relative threshold on diagonal of R

% Read basic information
fid=fopen(fin,'r');
nband=fscanf(fid,'%d',1);
nk=fscanf(fid,'%d',1);
nwan=fscanf(fid,'%g',1);
if nwan >= 1
    nwan = round(nwan);
end
feig = fscanf(fid,'%s',1);
mu = fscanf(fid,'%g',1);
sigma = fscanf(fid,'%g',1);
kpts = zeros(nk,3);
for k = 1 : nk
  kpts(k,:) = fscanf(fid,'%g',3);
end
fclose(fid);

E = load(feig);
eigKpts = reshape(E(:,3),nband,nk);

[check, gamma_idx] = min(sum(abs(kpts),2));
if check ~= 0
    disp('no gamma point found')
end


occ = (1/2)*erfc((eigKpts-mu)/sigma);

[Psi,nx,ny,nz] = readUNK(gamma_idx);

disp('Computing columns...');
[~, ~, perm]=qr((Psi*diag(occ(:,gamma_idx)))','vector');  
if nwan >= 1
    cols = perm(1:nwan);
else
    nwan = find(abs(diag(R))/abs(R(1,1)) < nwan,1);
    cols = perm(1:nwan);
end


[Xs,Ys,Zs] = ndgrid(0:nx-1,0:ny-1,0:nz-1);
Xs = Xs/nx;
Ys = Ys/ny;
Zs = Zs/nz;
rposFrac = [Xs(:) , Ys(:) , Zs(:)];

% The phase factor is important for approximating integration in the
% global domain.
disp('Performing projections...');
Amn = cell(nk,1);
for k = 1 : nk
    
    centerFrac = rposFrac(cols,:);
    centerFrac = centerFrac - round(centerFrac);
    phase = exp(1i*2*pi*centerFrac*kpts(k,:)');

    Psi = readUNK(k);
    Psi = diag(phase)*Psi(cols,:)*diag(occ(:,k));

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


