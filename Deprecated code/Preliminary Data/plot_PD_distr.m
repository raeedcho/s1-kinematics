function plot_PD_distr(pds,num_bins)
% plots distribution of PDs
pd_bin = linspace(-pi, pi, num_bins);
freq = histc(pds,pd_bin);

freq_norm = freq/max(freq);

freq_norm(end) = freq_norm(1);
pd_bin(end) = pd_bin(1);

% plot polar pd distribution
h_pol = polar(pd_bin(:),freq_norm(:),'b-');
set(h_pol,'LineWidth',3)