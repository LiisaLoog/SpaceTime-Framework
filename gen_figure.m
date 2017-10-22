%--------------------------------------------------------------------------
% Anders Eriksson 2017
%
% Matlab script for plotting the output from SimChangeMigrate
%
%--------------------------------------------------------------------------

clearvars;

%--------------------------------------------------------------------------
% Load simulation output
%--------------------------------------------------------------------------

d = load('sim_output.txt');
t = d(:,1)*0.0025; % Convert time from generations to unit of 1,000 years.
x = d(:,2); y=d(:,3); g = d(:,4:end);
clear d;

g = g - repmat(mean(g),size(g,1),1);

% 
T = ceil(max(t)); L = 4; ts = (0:0.2:T-L)'; nT = length(ts);
ta = ts + L/2; 

n = 101; 
w = linspace(-4,4,n)'; v = atan(10.^w);

iMax = zeros(nT,1);
pCorr = zeros(n,nT);
for i_t = 1:nT
	is = find(ts(i_t) <= t & t <= ts(i_t)+L);
	
	Dg2 = g(is,:)*g(is,:)'/length(is); Cd = diag(Dg2); e = Dg2(:,1)*0+1; Dg2 = Cd*e' + e*Cd' - 2*Dg2;
	Dt2 = (t(is)*e'-e*t(is)').^2; 
	Ds2 = (x(is)*e'-e*x(is)').^2 + (y(is)*e'-e*y(is)').^2;

	Y   = sqrt(abs(Dg2(:)));
	X_S = sqrt(Ds2(:))*100; 
	X_T = sqrt(Dt2(:));

    for k = 1:n
		pCorr(k,i_t) = corr(Y,sqrt(cos(v(k))*X_S+sin(v(k))*X_T)); 
    end
	[~,iMax(i_t)] = max(pCorr(:,i_t));
end

figure(1); clf;
imagesc(ta-T+5,w,pCorr); axis xy; box off; colorbar;
hold on; plot(ta-T+5,w(iMax),'k-','linewidth',2); hold off;
xlabel('Time (years BP)');
ylabel('log_{10} S_{max}');
hold on; h = plot([[-25 -15]',[-15 -10]',[-10 -0]']+5,[1 1]'*(log10([0.0002 0.001 0.005])+3.3),'k--','linewidth',2); hold off;
axis([-17 -2 -1 2])
set(gca,'xtick',(-16:2:-2),'XTickLabel',abs((-16:2:-2)*1000))
set(gca,'fontsize',12)
