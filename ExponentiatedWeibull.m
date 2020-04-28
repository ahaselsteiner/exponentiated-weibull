% This software was written for the publication
% "Predicting wave heights for marine design by prioritizing extreme 
% events in a global model" by A.F. Haselsteiner and K-D. Thoben, see
% https://arxiv.org/pdf/1911.12835.pdf .

classdef ExponentiatedWeibull < handle
% The exponentiated Weibull distribution is a 3-parameter 
% probability distribution. See, e.g. https://en.wikipedia.org/wiki/Exponentiated_Weibull_distribution
   properties
      Alpha % Scale parameter.
      Beta % Shape parameter #1.
      Delta % Shape parameter #2.
      BootstrapParm % Parameters estimated using bootstrap.
      ParameterSE % Parameters' standard error estimated using bootstrapping.
      ParameterCI % Parameters'confidence interval standard error estimated using bootstrapping.
   end
   
   methods
      function obj = ExponentiatedWeibull(alpha, beta, delta)
         if nargin > 2
            obj.Alpha = alpha; 
            obj.Beta = beta; 
            obj.Delta = delta; 
         end
      end
      
      function [pHat, pConfI] = fitDist(this, sample, method, varargin)
          if nargin < 3
              method = 'MLE'; % Maximum likelihood estimation.
          end
          
          p = inputParser;
          addOptional(p,'alpha', nan, @isnumeric);
          addOptional(p,'beta', nan, @isnumeric);
          addOptional(p,'delta', nan, @isnumeric);
          parse(p, varargin{:})
          isFixed(1) = ~isnan(p.Results.alpha);
          isFixed(2) = ~isnan(p.Results.beta);
          isFixed(3) = ~isnan(p.Results.delta);
            
          if method == 'MLE' % Maximum likelihood estimation.
              start = [1 1 5.001];
              lb = [0 0 0];
              ub = [100 100 100]; % MLE does not not converge for dataset A with limit (inf inf inf).
              if sum(isFixed) == 0
              [pHat, pConfI] = mle(sample, 'pdf', @(x, alpha, beta, delta) ...
                  this.pdf(sample, alpha, beta, delta), ...
                  'start', start, 'lower', lb, 'upper', ub);
              elseif sum(isFixed) == 1 || sum(isFixed) == 2
                  error('Error. Not implemented yet.')
              else
                  error('Error. At least one parameter needs to be free to fit it.')
              end
          elseif method == 'WLS' % Weighted least squares.
              n = length(sample);
              i = [1:n]';
              pi = (i - 0.5) ./ n;
              xi = sort(sample);
              delta0 = 2;
              if sum(isFixed) == 0
                  [deltaHat, WLSError] = fminsearch(@(delta) ...
                      estimateAlphaBetaWithWLS(delta, xi, pi), delta0);
                  [temp, pHat] = estimateAlphaBetaWithWLS(deltaHat, xi, pi);
              elseif sum(isFixed) == 1
                  if isFixed(3) == 1
                      deltaHat = p.Results.delta;
                      [temp, pHat] = estimateAlphaBetaWithWLS(deltaHat, xi, pi);
                  else
                      error('Error. Fixing alpha or beta is not implemented yet.')
                  end
              elseif sum(isFixed) == 2
                  error('Error. Not implemented yet.')
              else
                  error('Error. At least one parameter needs to be free to fit it.')
              end
          else
              error('Error. The input "method" must be either "MLE" or "WLS".')
          end
          this.Alpha = pHat(1);
          this.Beta = pHat(2);
          this.Delta = pHat(3);
      end
      
      function [parmHat, pStd, pCi] = fitDistAndBootstrap(this, ...
                sample, method, B, alpha)
          % Estimates the parameters of the distribution and estimate the 
          % parameters' uncertainty using bootstrapping.
          %
          % For bootstrapping, see, for example, "An introduction to the
          % bootstrap" by B. Efron and R. J. Tibshirani (1993).
          if nargin < 5
              alpha = 0.05;
          end
                    bAlphas = nan(B, 1);
          bBetas = nan(B, 1);
          bDeltas = nan(B, 1);
          for i = 1:B
              bSample = datasample(sample, length(sample), 'replace', true);
              bParmHat = this.fitDist(bSample, method);
              bAlphas(i) = bParmHat(1);
              bBetas(i) = bParmHat(2);
              bDeltas(i) = bParmHat(3);
          end
         
          % The index of the interval is chosen as in Efron and Tibshirani
          % (1993), p. 160.
          iLower = floor((B + 1) * (alpha / 2));
          iUpper = B + 1 - iLower;
          
          % Compute the estimators' standard deviations, see Eq. 2.3 in 
          % Efron and Tibshirani (1993).
          pStd = [std(bAlphas), std(bBetas), std(bDeltas)];
          
          % Compute the 1 - alpha intervals based on bootstrap percentiles 
          % (see Efron and Tibshirani (1993), pp. 168). 
          sortedAlphas = sort(bAlphas);
          sortedBetas = sort(bBetas);
          sortedDeltas = sort(bDeltas);
          pCi = ...
              [sortedAlphas(iLower), sortedBetas(iLower), sortedDeltas(iLower);
               sortedAlphas(iUpper), sortedBetas(iUpper), sortedDeltas(iUpper)];
          
          parmHat = this.fitDist(sample, method); % Calling fitDist also sets the class' parameters.
          this.BootstrapParm = [bAlphas, bBetas, bDeltas];
          this.ParameterSE = pStd;
          this.ParameterCI = pCi;
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
          x(x < 0) = 0;
          F = (1 - exp( -1 .* (x ./ this.Alpha).^this.Beta)).^this.Delta;
      end
      
      function x = icdf(this, p)
          % Inverse cumulative distribution function.
          p(p > 1) = NaN;
          p(p < 0) = NaN;
          x = this.Alpha .* (-1 .* log(1 - p.^(1 ./ this.Delta))).^(1 ./ this.Beta);
      end
      
      function x = drawSample(this, n)
          if n < 2
              n = 1;
          end
          p = rand(n, 1);
          x = this.icdf(p);
      end
      
      function mk = kMoment(this, k)
          % kth moment of the distribution.
          xk = @(x, k) this.icdf(x).^k;
          fun = @(k) integral(@(x) xk(x, k), 0, 1);
          mk = fun(k);
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
              figure();
          end
          if nargin < 5
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
          ax = gca;
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
