classdef ExponentiatedWeibull < handle
% The exponentiated Weibull distribution is a 3-parameter 
% probability distribution. See, e.g. https://en.wikipedia.org/wiki/Exponentiated_Weibull_distribution
   properties
      Alpha % Scale parameter.
      Beta % Shape parameter #1.
      ShapeTwo % Shape parameter #2.
   end
   
   methods
      function obj = ExponentiatedWeibull(alpha, beta, shapeTwo)
         if nargin > 2
            obj.Alpha = alpha; 
            obj.Beta = beta; 
            obj.ShapeTwo = shapeTwo; 
         end
      end
      
      function parmHat = fitDist(this, sample, method)
          if nargin < 2
              method = 'MLE'; % Maximum likelihood estimation.
          end
          if method == 'MLE' % Maximum likelihood estimation.
              start = [1 1 5.001];
              lb = [0 0 0];
              ub = [100 100 100]; % MLE does not not converge for dataset A with limit (inf inf inf).
              parmHat = mle(sample, 'pdf', @(x, alpha, beta, shapeTwo) ...
                  this.pdf(sample, alpha, beta, shapeTwo), ...
                  'start', start, 'lower', lb, 'upper', ub);
          elseif method == 'WLS' % Weighted least squares.
              n = length(sample);
              i = [1:n]';
              pi = (i - 0.5) ./ n;
              xi = sort(sample);
              shapeTwo0 = 2;
              [shapeTwoHat SQR_min] = fminsearch(@(shapeTwo) ...
                  estimateAlphaBetaWithWLS(shapeTwo, xi, pi), shapeTwo0);
              [temp parmHat] = estimateAlphaBetaWithWLS(shapeTwoHat, xi, pi);
          else
              error('Error. The input "method" must be either "MLE" or "WLS".')
          end
          this.Alpha = parmHat(1);
          this.Beta = parmHat(2);
          this.ShapeTwo = parmHat(3);
      end
      
      function f = pdf(this, x, alpha, beta, shapeTwo)
          % Probability density function.
          pdf = @(x, alpha, beta, shapeTwo) shapeTwo .* beta ./ alpha .* (x ./ alpha).^(beta - 1) ...
                .* (1 - exp(-1 * (x ./ alpha).^beta)).^(shapeTwo - 1) .* exp(-1 .* (x ./ alpha).^beta);
          if nargin < 3
              f = pdf(x, this.Alpha, this.Beta, this.ShapeTwo);
          else
              f = pdf(x, alpha, beta, shapeTwo);
          end              
      end
      
      function F = cdf(this, x)
          % Cumulative distribution function.
          F = (1 - exp( -1 .* (x ./ this.Alpha).^this.Beta)).^this.ShapeTwo;
      end
      
      function x = icdf(this, p)
          % Inverse cumulative distribution function.
          x = this.Alpha * (-1 * log(1 - p.^(1 ./ this.ShapeTwo))).^(1 ./ this.Beta);
      end
      
      function val = negativeloglikelihood(this, x)
          % Negative log-likelihood value (as a metric of goodness of fit).
          val = sum(-log(pdf(x, this.Alpha, this.Beta, this.ShapeTwo)));
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

function [SQR, parmHat] = estimateAlphaBetaWithWLS(shapeTwo, hsi, pi) 
    % Transform variables to Weibull paper [see, Scholz (2008), but in terms of notation switch "alpha" and "beta"]
    % First, transform the probability pi:
    pstar_i = log10(-log(1 - pi.^(1 / shapeTwo)));
    % Then, transfoorm the significant wave height hsi:
    yi = log10(hsi);

    % Define the weights.
    wi = hsi.^2 / sum(hsi.^2);

    % Estimate the two Weibull parameters, alphahat and betahat
    % shapeTwo
    pstarbar = sum(wi .* pstar_i);
    ybar = sum(wi .* yi);
    bhat = (sum(wi .* pstar_i .* yi) - pstarbar .* ybar) / ...
        (sum(wi .* pstar_i.^2) -  pstarbar^2);
    ahat = ybar - bhat * pstarbar;
    betaHat = 1 / bhat;
    alphaHat = 10^ahat;

    % Create a distribution object and compute estimates hs_hat_WLS
    hs_hat_WLS = alphaHat * (-1 * log(1 - pi.^(1 / shapeTwo))).^(1 / betaHat);

    % Compute the weighted sum of square residuals
    SQR = sum(wi .* (hsi - hs_hat_WLS).^2);
    parmHat = [alphaHat, betaHat, shapeTwo];
end
