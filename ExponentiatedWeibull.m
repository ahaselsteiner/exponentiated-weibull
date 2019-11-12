classdef ExponentiatedWeibull < handle
% The exponentiated Weibull distribution is a 3-parameter 
% probability distribution. See, e.g. https://en.wikipedia.org/wiki/Exponentiated_Weibull_distribution
   properties
      Alpha % Scale parameter.
      Beta % Shape parameter #1.
      Delta % Shape parameter #2.
   end
   
   methods
      function obj = ExponentiatedWeibull(alpha, beta, delta)
         if nargin > 2
            obj.Alpha = alpha; 
            obj.Beta = beta; 
            obj.Delta = delta; 
         end
      end
      
      function [pHat pConfI] = fitDist(this, sample, method)
          if nargin < 2
              method = 'MLE'; % Maximum likelihood estimation.
          end
          if method == 'MLE' % Maximum likelihood estimation.
              start = [1 1 5.001];
              lb = [0 0 0];
              ub = [100 100 100]; % MLE does not not converge for dataset A with limit (inf inf inf).
              [pHat pConfI] = mle(sample, 'pdf', @(x, alpha, beta, delta) ...
                  this.pdf(sample, alpha, beta, delta), ...
                  'start', start, 'lower', lb, 'upper', ub);
          elseif method == 'WLS' % Weighted least squares.
              n = length(sample);
              i = [1:n]';
              pi = (i - 0.5) ./ n;
              xi = sort(sample);
              delta0 = 2;
              [deltaHat, WLSError] = fminsearch(@(delta) ...
                  estimateAlphaBetaWithWLS(delta, xi, pi), delta0);
              [temp, pHat] = estimateAlphaBetaWithWLS(deltaHat, xi, pi);
          else
              error('Error. The input "method" must be either "MLE" or "WLS".')
          end
          this.Alpha = pHat(1);
          this.Beta = pHat(2);
          this.Delta = pHat(3);
      end
      
      function f = pdf(this, x, alpha, beta, delta)
          % Probability density function.
          pdf = @(x, alpha, beta, delta) delta .* beta ./ alpha .* (x ./ alpha).^(beta - 1) ...
                .* (1 - exp(-1 * (x ./ alpha).^beta)).^(delta - 1) .* exp(-1 .* (x ./ alpha).^beta);
          if nargin < 3
              f = pdf(x, this.Alpha, this.Beta, this.Delta);
          else
              f = pdf(x, alpha, beta, delta);
          end              
      end
      
      function F = cdf(this, x)
          % Cumulative distribution function.
          F = (1 - exp( -1 .* (x ./ this.Alpha).^this.Beta)).^this.Delta;
      end
      
      function x = icdf(this, p)
          % Inverse cumulative distribution function.
          x = this.Alpha .* (-1 .* log(1 - p.^(1 ./ this.Delta))).^(1 ./ this.Beta);
      end
      
      function x = drawSample(this, n)
          if n < 2
              n = 1;
          end
          p = rand(n, 1);
          x = this.icdf(p);
      end
      
      function val = negativeloglikelihood(this, x)
          % Negative log-likelihood value (as a metric of goodness of fit).
          val = sum(-log(pdf(x, this.Alpha, this.Beta, this.Delta)));
      end
      
      function mae = meanabsoluteerror(this, sample, pi)
          % Mean absolute error (as a measure of goodness of fit).
          n = length(sample);
          if nargin < 3
              i = [1:n]';
              pi = (i - 0.5) ./ n;
          end
          xi = sort(sample);
          xhati = this.icdf(pi); % The prediction.
          mae = sum(abs(xi - xhati)) / n;
      end
      
      function ax = qqplot(this, sample, qqFig, qqAx, lineColor)
          if nargin > 2
              set(0, 'currentfigure', qqFig);
              set(qqFig, 'currentaxes', qqAx);
          else
              qqFig = figure();
          end
          if nargin < 4
              lineColor = [0 0 0];
          end
          n = length(sample);
          i = [1:n]';
          pi = (i - 0.5) ./ n;
          xi = sort(sample);
          xhati = this.icdf(pi); % The prediction.
          hold on
          plot(xhati, xi, 'kx'); 
          plot(xhati, xhati, 'color', lineColor, 'linewidth', 1.5);
          xlabel('Theoretical quantiles');
          ylabel('Ordered values');
      end
   end
end

function [WLSError, pHat] = estimateAlphaBetaWithWLS(delta, xi, pi) 
    % First, transform xi and pi to get a linear relationship.
    xstar_i = log10(xi);
    pstar_i = log10(-log(1 - pi.^(1 / delta)));

    % Define the weights.
    wi = xi.^2 / sum(xi.^2);

    % Estimate the parameters alphaHat and betaHat.
    pstarbar = sum(wi .* pstar_i);
    xstarbar = sum(wi .* xstar_i);
    bHat = (sum(wi .* pstar_i .* xstar_i) - pstarbar .* xstarbar) / ...
        (sum(wi .* pstar_i.^2) -  pstarbar^2);
    aHat = xstarbar - bHat * pstarbar;
    alphaHat = 10^aHat;
    betaHat = 1 / bHat;
    pHat = [alphaHat, betaHat, delta];
    
    % Compute the weighted least squares error.
    xiHat = alphaHat * (-1 * log(1 - pi.^(1 / delta))).^(1 / betaHat);
    WLSError = sum(wi .* (xi - xiHat).^2);
end
