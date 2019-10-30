pdTrue = ExponentiatedWeibull(1, 1, 2);
n = 100000;
nOfSamples = 100;

alphaEstimated = nan(nOfSamples, 1);
betaEstimated = nan(nOfSamples, 1);
deltaEstimated = nan(nOfSamples, 1);
for i = 1:nOfSamples
    sample = pdTrue.drawSample(n);
    pdEstimated(i) = ExponentiatedWeibull();
    pdEstimated(i).fitDist(sample, 'WLS');
    alphaEstimated(i) = pdEstimated(i).Alpha;
    betaEstimated(i) = pdEstimated(i).Beta;
    deltaEstimated(i) = pdEstimated(i).Delta;
end

fig = figure('position', [100 100 400, 130]);
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
