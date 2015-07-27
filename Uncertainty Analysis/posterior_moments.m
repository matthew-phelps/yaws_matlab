 function [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(xx,info,mh_conf_sig)
0002 % Computes posterior mean, median, variance, HPD interval, deciles, and density from posterior draws.
0003 %
0004 % INPUTS
0005 %    xx            [double]    Vector of posterior draws (or prior draws ;-)
0006 %    info          [integer]   If equal to one the posterior density is estimated.
0007 %    mh_config_sig [double]    Scalar between 0 and 1 specifying the size of the HPD interval.
0008 %
0009 %
0010 % OUTPUTS
0011 %    post_mean     [double]    Scalar, posterior mean.
0012 %    post_median   [double]    Scalar, posterior median.
0013 %    post_var      [double]    Scalar, posterior variance.
0014 %    hpd_interval  [double]    Vector (1*2), Highest Probability Density interval
0015 %    post_deciles  [double]    Vector (9*1), deciles of the posterior distribution.
0016 %    density       [double]    Matrix (n*2), non parametric estimate of the posterior density. First and second
0017 %                              columns are respectively abscissa and ordinate coordinates.
0018 %
0019 % SPECIAL REQUIREMENTS
0020 %    Other matlab routines distributed with Dynare: mh_optimal_bandwidth.m
0021 %                                                   kernel_density_estimate.m.
0022 %
0023 
0024 % Copyright (C) 2005-2011 Dynare Team
0025 %
0026 % This file is part of Dynare.
0027 %
0028 % Dynare is free software: you can redistribute it and/or modify
0029 % it under the terms of the GNU General Public License as published by
0030 % the Free Software Foundation, either version 3 of the License, or
0031 % (at your option) any later version.
0032 %
0033 % Dynare is distributed in the hope that it will be useful,
0034 % but WITHOUT ANY WARRANTY; without even the implied warranty of
0035 % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
0036 % GNU General Public License for more details.
0037 %
0038 % You should have received a copy of the GNU General Public License
0039 % along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
0040 
0041 xx = xx(:);
0042 xx = sort(xx);
0043 
0044 
0045 post_mean = mean(xx);
0046 post_median = median(xx);
0047 post_var = var(xx);
0048 
0049 number_of_draws = length(xx);
0050 hpd_draws = round((1-mh_conf_sig)*number_of_draws);
0051 
0052 if hpd_draws>2
0053     kk = zeros(hpd_draws,1);
0054     jj = number_of_draws-hpd_draws;
0055     for ii = 1:hpd_draws
0056         kk(ii) = xx(jj)-xx(ii);
0057         jj = jj + 1;
0058     end
0059     [kmin,idx] = min(kk);
0060     hpd_interval = [xx(idx) xx(idx)+kmin];
0061 else
0062     hpd_interval=NaN(1,2);    
0063 end
0064 if length(xx)>9
0065     post_deciles = xx([round(0.1*number_of_draws) ...
0066                        round(0.2*number_of_draws) ...
0067                        round(0.3*number_of_draws) ...
0068                        round(0.4*number_of_draws) ...
0069                        round(0.5*number_of_draws) ...
0070                        round(0.6*number_of_draws) ...
0071                        round(0.7*number_of_draws) ...
0072                        round(0.8*number_of_draws) ...
0073                        round(0.9*number_of_draws)]);
0074 else
0075     post_deciles=NaN(9,1);
0076 end
0077 density = [];
0078 if info
0079     number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
0080     if post_var > 1e-12
0081         bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
0082         kernel_function = 'gaussian';     % Gaussian kernel for Fast Fourrier Transform approximaton.
0083         optimal_bandwidth = mh_optimal_bandwidth(xx,number_of_draws,bandwidth,kernel_function);
0084         [density(:,1),density(:,2)] = kernel_density_estimate(xx,number_of_grid_points,...
0085                                                           number_of_draws,optimal_bandwidth,kernel_function);
0086     else
0087         density = NaN(number_of_grid_points,2);
0088     end
0089 end