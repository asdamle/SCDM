clear; close all;
load insul1280.mat
cutoff = 160;
Psi = V(:,1:128);
P = Psi*(Psi');
f = (1/2)*erfc((D(1:128)-.2)/.5);

[Phi , piv] = scdm_entangled(Psi,diag(f),64);

Psi_small = Psi(1:cutoff,piv<=cutoff);
Phi_small = Phi(1:cutoff,piv<=cutoff);

Phi_small_many = Phi(1:(cutoff/2),piv<=(cutoff/2));

ymax = max([Psi_small(:); Phi_small(:)])*1.1;
ymin = min([Psi_small(:); Phi_small(:)])*1.1;

figure
plot(Phi_small(:,1),'k','LineWidth',3)
% ylim([ymin ymax])
axis off
set(gcf, 'Position', [0, 0, 500, 500])
fname = 'Ecase1_phi.pdf';
save2pdf(fname,gcf,600);

figure
plot(Phi_small_many,'LineWidth',3)
% ylim([ymin ymax])
axis off
set(gcf, 'Position', [0, 0, 500, 500])
fname = 'Ecase1_phi_many.pdf';
save2pdf(fname,gcf,600);

% figure
% plot(Psi_small,'LineWidth',3)
% % ylim([ymin ymax])
% axis off
% set(gcf, 'Position', [0, 0, 500, 500])
% fname = 'Ecase1_psi.pdf';
% save2pdf(fname,gcf,600);




