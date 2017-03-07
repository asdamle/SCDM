clear; close all;
load insul1280.mat
cutoff = 160;
Psi = V(:,1:128);
P = Psi*(Psi');
f = exp((D(1:128)-.2).^2/2);

[Phi , piv] = scdm_entangled(Psi,diag(f),32);

Psi_small = Psi(1:cutoff,piv<=cutoff);
Phi_small = Phi(1:cutoff*2,piv<=cutoff);

ymax = max([Psi_small(:); Phi_small(:)])*1.1;
ymin = min([Psi_small(:); Phi_small(:)])*1.1;

figure
plot(Phi_small(:,1),'k','LineWidth',3)
% ylim([ymin ymax])
axis off
set(gcf, 'Position', [0, 0, 500, 500])
fname = 'Ecase2_phi.pdf';
save2pdf(fname,gcf,600);

figure
plot(Phi_small(1:cutoff,:),'LineWidth',3)
% ylim([ymin ymax])
axis off
set(gcf, 'Position', [0, 0, 500, 500])
fname = 'Ecase2_phi_many.pdf';
save2pdf(fname,gcf,600);

% figure
% plot(Psi_small,'LineWidth',3)
% ylim([ymin ymax])
% axis off
% set(gcf, 'Position', [0, 0, 500, 500])
% fname = 'Ecase2_psi.pdf';
% save2pdf(fname,gcf,600);



