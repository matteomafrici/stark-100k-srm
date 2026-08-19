function [a, Inc_a, n, Inc_n, R2] = Uncertainty(p, rb)

%
% Original file: Incertezze.m, V 1.02
%
% 04 August 2011 - Dossi
%
% syntax: [a, Inc_a, n, Inc_n, R2] = Uncertainty(p, rb)
%
% rb = a * p^n -> log(rb) = log(a*p^n) = log(a) + n * log(p)
% Y = log(rb); X = log(p); q = log(a); m = n
% Y = m*X + q
%
X = log(p);                                                         % convert pressure values to log scale
N = length(p);                                                      % number of pressure samples
Y = log(rb);                                                        % convert burn-rate values to log scale

delta = N * sum(X.^2)-(sum(X))^2;                                   % denominator for the q and m estimates
m = ( N*sum(X.*Y) - sum(X)*sum(Y) ) / delta;                        % slope of the fit line Y = mX + q, m = n
n = m;                                                              % from the line fit above
q = ( sum(X.^2)*sum(Y) - sum(X)*sum(X.*Y)) / delta;                 % intercept of the fit line Y = mX + q, m = n

Y_eval = m .* X + q;                                                % values evaluated from the fitted line
sig = ( ( sum( ( Y - Y_eval ).^2) ) / (N-2) )^0.5;                  % RMS deviation of the data from the fitted line

Inc_m = sig * (N/delta)^0.5;                                        % uncertainty of the slope m (same as n)
Inc_n = Inc_m;                                                      % from the line fit above
Inc_q = sig * ( sum( X.^2 ) / delta )^0.5;                          % uncertainty of the intercept q

a = exp(q);                                                         % back out a from the intercept
a_new_up = a * exp(Inc_q);                                          % upper bound of the a uncertainty
a_new_down = a / exp(Inc_q);                                        % lower bound of the a uncertainty
Inc_a = ( a_new_up - a_new_down) / 2;                               % symmetric uncertainty of a


M_Y = mean(Y);                                                      % mean of Y
dev_reg = sum ( (Y_eval - M_Y).^2 );                                % regression (explained) deviance
dev_tot = sum ( (Y - M_Y).^2 );                                     % total deviance
R2 = dev_reg / dev_tot;                                             % coefficient of determination

return