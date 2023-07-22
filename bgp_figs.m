% This code generates Figure 3 of our comment.

% Figure 3
%
% First, obtain BGP productivity growth and BGP profit share for each
% discount rate in rvec, exactly as in LMS's code.  Lines 15-30 of bgp_figs.m match
% Lines 3-18 of LMS's main script calibration_EMA_submit.m exactly, without modification.
% (These results are used to generate the black solid line in our Figure 3.)

close all; clear all; n=50; xinit = zeros(1,2*n); xinit(n+1) = 1; flatpi=1;
sig=12; lamb = 1.21; rvec = [6:-0.1:0.1,0.05,0.02];
pivec=compute_pi_fast(sig,lamb,n);
pivec(n+1+flatpi:end) = pivec(n+1+flatpi); pivec(1:n+1-flatpi)=pivec(n+1-flatpi);
c=33.3569^2;
pi=pivec*c;
kap=3.9345;
pishr=pivec(n+1:-1:1)+pivec(n+1:end);

tic;for i=1:length(rvec)
    [xvec, muvec, ~, gvec(i),~,flag] = gen_compute_eqm(lamb,pi,1,kap,rvec(i),xinit);
    if flag>0;xinit=xvec;end
    invcost=[xvec.^2/c,0];
    avgmu(i)=(0:n)*muvec';
    LI(i)=(pishr-invcost(n+1:end)-invcost(n+1:-1:1))*muvec';
end;toc;

black = 'k';

% Second, obtain BGP productivity growth and BGP profit share for each
% discount rate in rvec, assuming that the cost of achieving an 
% investment success rate eta_s is (100/sqrt(2)/33.3569 eta_s)^2. 
% (LMS's code solves for BGP investment success rates 
% using an investment cost (100/sqrt(2)/33.3569 eta_s)^2,
% but uses an investment cost twice as large to obtain BGP profit shares.)
% These results are used to obtain the red dashed line in our Figure 3.

c = 100/sqrt(2)/33.3569;
xinit = zeros(1,2*n); xinit(n+1) = 1/100;

for i=1:length(rvec)
    [xvec_mod, muvec_mod, ~, gvec_mod(i),~,flag_mod] = gen_compute_eqm_correct(lamb,pivec,1,kap/100,rvec(i)/100,xinit,c);
    if flag_mod>0;xinit=xvec_mod;end
    invcost_mod=[c^2 * xvec_mod.^2,0];
    avgmu_mod(i)=(0:n)*muvec_mod';
    LI_mod(i)=(pishr-invcost_mod(n+1:end)-invcost_mod(n+1:-1:1))*muvec_mod';
end

gvec_outer = [gvec;  100*gvec_mod];
LI_outer = 1*[LI;  LI_mod];

[status,msg,msgID] = mkdir('figures_comment');

% Our Figure 3

figure; 
set(gcf, 'PaperUnits', 'inches');
x_width=5;
y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 

ax=  subplot(1,1,1);
p1 = plot(1*(rvec+gvec_outer(2, :)),LI_outer(2, :),'--','LineWidth',2,'Color','r'); 
hold on 
p2 = plot(1*(rvec+gvec_outer(1, :)), LI_outer(1, :), '-k', 'LineWidth', 1.75);
xtickformat(ax, 'percentage');
ylim([0.135 0.2]); 
ax.YGrid = 'on';
xlim([0, 7])
box off; 
set(gca,'ytick', [0.14 0.15 0.16 0.17, 0.18, 0.19, 0.20]); 
xlabel('Interest rate', 'Interpreter', 'latex');
title('Profit share', 'Interpreter', 'latex')


l1 = legend([p1, p2], ...
           ["Correcting inconsistency in LMS code", "As reported in LMS Figure 4"], 'Interpreter', 'latex', ... 
           'Location', 'northeast', 'FontSize', 10);
l1.ItemTokenSize = [10; 5];
set(l1, 'box', 'off')
saveas(gcf, "figures_comment/inv_cost_profit_share.eps", 'epsc');


