% Test #1: Does the implementation reproduce results from the literature? 
% We use the parameters from Pal, M., Ali, M. M., Woo, J. (2006): 
% Exponentiated Weibull distribution. Statistica 66(1), 139-147.
% In Figure 1 the "alpha=2"-curve is used. Note that they use a different
% parameterization.
pd = ExponentiatedWeibull(1/0.5, 2, 2);
x = [0:0.01:6];
f = pd.pdf(x);
fig1 = figure('position', [100 100 450 280]);
plot(x, f);
message = sprintf(['In Pal et al. (2006), Fig. 1 \n' ...
     'the PDF peaks at ~2 m with density of ~0.45.']);
text(2, 0.48, message, 'horizontalalignment', ...
    'left', 'verticalalignment', 'bottom', 'fontsize', 8);add 
ylabel('Density (-)');
xlabel('Significant wave height (m)');
box off
ylim([0 0.6]);

% Test #2: Does the parameter estimation work correctly?
pdTrue = ExponentiatedWeibull(1, 1, 2);
n = 100000;
nOfSamples = 50;
alphaEstimated = nan(nOfSamples, 1);
betaEstimated = nan(nOfSamples, 1);
deltaEstimated = nan(nOfSamples, 1);
pdEstimated = ExponentiatedWeibull.empty(nOfSamples, 0);
for i = 1:nOfSamples
    sample = pdTrue.drawSample(n);
    pdEstimated(i) = ExponentiatedWeibull();
    pdEstimated(i).fitDist(sample, 'WLS');
    alphaEstimated(i) = pdEstimated(i).Alpha;
    betaEstimated(i) = pdEstimated(i).Beta;
    deltaEstimated(i) = pdEstimated(i).Delta;
end

fig2 = figure('position', [100 100 500, 230]);
subplot(1, 3, 1)
hold on
plot([0.5 1.5], [1 1], '-k')
boxplot(alphaEstimated, {'alpha'})
text(1.15, 1, [num2str(mean(alphaEstimated), '%1.3f') '+-' ...
    num2str(std(alphaEstimated), '%1.3f')], 'fontsize', 8, ...
    'verticalalignment', 'bottom'); 
box off

subplot(1, 3, 2)
hold on
plot([0.5 1.5], [1 1], '-k')
boxplot(betaEstimated, {'beta'})
text(1.15, 1, [num2str(mean(betaEstimated), '%1.3f') '+-' ...
    num2str(std(betaEstimated), '%1.3f')], 'fontsize', 8, ...
    'verticalalignment', 'bottom');
box off

subplot(1, 3, 3)
hold on
plot([0.5 1.5], [2 2], '-k')
boxplot(deltaEstimated, {'delta'})
text(1.15, 2, [num2str(mean(deltaEstimated), '%1.3f') '+-' ...
    num2str(std(deltaEstimated), '%1.3f')], 'fontsize', 8, ...
    'verticalalignment', 'bottom'); 
box off
suptitle(['Parameter estimation, true parameters: ' ...
    num2str(pdTrue.Alpha) ', ' num2str(pdTrue.Beta) ', ' num2str(pdTrue.Delta)]);
