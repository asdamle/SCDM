clear; close all;
load insul1280.mat
Q = diag(Q);
Psi = V(:,1:32);
P = Psi*(Psi');
[Phi , piv] = scdm_entangled(Psi,eye(32),32);

Psi_small = Psi(:,piv<=160);
Phi_small = Phi(1:160,piv<=160);

ymax = max([Psi_small(:); Phi_small(:)])*1.1;
ymin = min([Psi_small(:); Phi_small(:)])*1.1;

figure
plot(Phi_small(:,1),'k','LineWidth',3)
% ylim([ymin ymax])
axis off
set(gcf, 'Position', [0, 0, 500, 500])
fname = 'isolated_phi.pdf';
save2pdf(fname,gcf,600);

figure
plot(Phi_small,'LineWidth',3)
% ylim([ymin ymax])
axis off
set(gcf, 'Position', [0, 0, 500, 500])
fname = 'isolated_phi_many.pdf';
save2pdf(fname,gcf,600);

figure
plot(Psi_small(:,1),'k','LineWidth',3)
% ylim([ymin ymax])
axis off
set(gcf, 'Position', [0, 0, 1000, 500])
fname = 'isolated_psi.pdf';
save2pdf(fname,gcf,600);

