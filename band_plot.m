function [] = band_plot(kpts,evals,kpath,kpath_loc,marker,varargin)


nvarargin = length(varargin);
stp = 1;
if( nvarargin >= 1 )
  stp = varargin{1};
end

% plots on current active figure with markers defined by marker string

num_evals = size(evals,1);
num_kpts = size(evals,2);

% first, find fixed points?

% kL = [0.50000  0.50000 0.5000];
% kG = [0.00000  0.00000 0.0000];
% kX = [0.50000 -0.50000 0.0000];
% kK = [0.37500 -0.37500 0.0000];
% 
% locL = find(sum(abs(kpts-repmat(kL,num_kpts,1)),2)<1e-6);
% locG = find(sum(abs(kpts-repmat(kG,num_kpts,1)),2)<1e-6);
% locX = find(sum(abs(kpts-repmat(kX,num_kpts,1)),2)<1e-6);
% locK = find(sum(abs(kpts-repmat(kK,num_kpts,1)),2)<1e-6);

locL = kpath_loc(kpath=='L');
locG = kpath_loc(kpath=='G');
locM = kpath_loc(kpath=='M');
locX = kpath_loc(kpath=='X');
locK = kpath_loc(kpath=='K');
locW = kpath_loc(kpath=='W');



% plot with spacing based on distance between points in space
x = sqrt(sum(diff(kpts).^2,2));
x = [0; cumsum(x)];
for k = 2:length(kpath_loc)-1
    adj = x(kpath_loc(k)+1)-x(kpath_loc(k));
    x((kpath_loc(k)+1):end) = x((kpath_loc(k)+1):end)-adj;
end

hold on
for k = 1:num_evals
  plot(x(1:stp:end),evals(k,1:stp:end),marker,'LineWidth',2)
end
    
% make it look nice
set(gca,'FontSize',28)

ax = gca;

ticks = [x(locL); x(locG); x(locM); x(locX); x(locK); x(locW)]; 

label = cell(1,length(ticks));
idx = 1;
for p = 1:length(locL)
    label{idx} = 'L';
    idx = idx+1;
end
for p = 1:length(locG)
    label{idx} = '\Gamma';
    idx = idx+1;
end
for p = 1:length(locM)
    label{idx} = 'M';
    idx = idx+1;
end
for p = 1:length(locX)
    label{idx} = 'X';
    idx = idx+1;
end
for p = 1:length(locK)
    label{idx} = 'K';
    idx = idx+1;
end
for p = 1:length(locW)
    label{idx} = 'W';
    idx = idx+1;
end
[ticks, idx] = sort(ticks,'ascend');
label = label(idx);
ax.XTick = ticks;
ax.XTickLabel= label;
ax.XGrid = 'on';
ax.GridAlpha = .5;
xlim([x(1) x(end)])
buffer = (max(evals(:))-min(evals(:)))/10;
ylim([min(evals(:))-buffer max(evals(:))+buffer])
box on
ylabel('Energy (eV)');
