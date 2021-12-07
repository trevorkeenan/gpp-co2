function [x] = press(D)
%PRESS Prediction error sum of squares.
% This m-file returns a useful residual scaling, the prediction error sum
% of squares (PRESS). According to Myers and Montgomery (2002, p 46-47), to 
% calculate PRESS, select an observation i. Fit the regression model to the
% remaining n-1 observations and use this equation to predict the withheld
% observation y_i. Denoting this predicted value by ye_(i), we may find the
% prediction error for point i as e_(i)=y_i - ye_(i).
% The prediction error is often called the ith PRESS residual. This
% procedure is repeated for each observation i = 1,2,...,n, producing a set
% of n PRESS residuals e_(1),e_(2),...,e_(n). Then the PRESS statistic is
% defined as the sum of squares of the n PRESS residuals as in,
%
%            PRESS = i_Sum_n e_(i)^2 = i_Sum_n [y_i - ye_(i)]^2
%
% Thus PRESS uses such possible subset of n-1 observations as an estimation
% data set, and every observation in turn is used to form a prediction data
% set. In the construction of this m-file, we use this statistical approach.
% However, as we have seen that calculating PRESS requires fitting n
% different regressions. But, also it is possible to calculate PRESS from
% the results of a single least squares fit to all n observations. It turns
% out that the ith PRESS residual is,
%
%                          e_(i) = e_i/(1 - h_ii)
%
% Thus, because PRESS is just the sum of the squares of the PRESS residuals,
% a simple computing formula is
%
%              PRESS = i_Sum_n [e_i/(1 - h_ii)]^2
%
% It is easy to see that the PRESS residual is just the ordinary residual
% weighted according to the diagonal elements of the hat matrix h_ii. Also,
% for all the interested people, here we just indicate, in an inactive form,
% this statistical approaching.
%
% Data points for which h_ii are large will have large PRESS residuals.
% These observations will generally be high influence points. Generally, a
% large difference between the ordinary residual and the PRESS residual will
% indicate a point where the model fits the data well, but a model built
% without that point predicts poorly. 
%
% This is also known as leave-??one-??out cross-????validation (LOOCV) in linear
% models as a measure of the accuracy. [an anon's suggestion]
% 
% In order to improve the matrix script for avoiding the squares condition
% number in the regression parameter estimation are by using a pivoted QR 
% factorization of X.
%
% Syntax: function x = press(D)
%
% Inputs:
%    D - matrix data (=[X Y]) (last column must be the Y-dependent variable).
%           (X-independent variables).
% Output:
%    x - prediction error sum of squares (PRESS).
%
% Example: 
% From example 2.1 from Myers and Montgomery (2002, p.23) we are interested 
% to calculate the prediction error sum of squares (PRESS). Data are,
%
%                      X1         X2         Y
%                   -----------------------------
%                      -1         -1        1004
%                       1         -1        1626
%                      -1       0.6667       852
%                       1       0.6667      1506
%                       0      -0.4444      1272
%                       0      -0.7222      1270
%                       0       0.6667      1269
%                      -1      -0.1667       903
%                       1      -0.1667      1555
%                       0         -1        1260
%                       0       0.9444      1146
%                       0      -0.1667      1276
%                       0          1        1225
%                    0.1667    -0.1667      1321
%                   -----------------------------
%
% Data matrix must be:
%  D=[-1 -1 1004;1 -1 1636;-1 0.6667 852;1 0.6667 1506;0 -0.4444 1272;
%  0 -0.7222 1270;0 0.6667 1269;-1 -0.1667 903;1 -0.1667 1555;0 -1 1260;
%  0 0.94444 1146;0 -0.1667 1276;0 1 1225;0.1667 -0.1667 1321];
%
% Calling on Matlab the function: 
%    x = press(D)
%
% Answer is:
%
% x = 2.2225e+004   (= 22,225.0)
%
% Created by A. Trujillo-Ortiz, R. Hernandez-Walls, A. Castro-Perez
%            and K. Barba-Rojo
%            Facultad de Ciencias Marinas
%            Universidad Autonoma de Baja California
%            Apdo. Postal 453
%            Ensenada, Baja California
%            Mexico.
%            atrujo@uabc.mx
%
% Copyright (C) April 02, 2007.
%
% ---We deeply thank to Bart and John D'Errico for call us the attention to
% improve the matrix script for avoiding the squares condition number
% in the regression parameter estimation by using a pivoted QR 
% factorization of X (08-26-2013)---
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A., R. Hernandez-Walls, K. Barba-Rojo, and 
%   A. Castro-Perez (2006). press:Prediction error sum of squares.
%   A MATLAB file. [WWW document]. URL http://www.mathworks.com/
%   matlabcentral/fileexchange/loadFile.do?objectId=14564
%
%  Reference:
%  Myers, R. H. and Montgomery, D. C. (2002), Response Surface Methodology:
%          Process and Product Optimization Using Designed Experiments. 2nd.
%          Ed. NY: John Wiley & Sons, Inc.
%

n = size(D,1);
c = [1:n];
I = reshape(c(repmat(1:n,n-1,1)),n,n-1);
I = flipud(I);

e = [];
for i = 1:n;
    idx = I(i,:)';
    DD = D(idx,:);
    [r c] = size(DD);
    n = r; %number of data
    Y = DD(:,c); %response vector
    X = [ones(n,1) DD(:,1:c-1)]; %design matrix
    [Q,R] = qr(X);
    b = R\(Q'*Y); %parameters estimation by using a pivoted QR 
                  %factorization of X
    ye = [1 D(i,1:c-1)]*b;
    ee = (D(i,end)-ye)^2;
    e = [e;ee];
end

x = sum(e); %prediction error sum of squares (PRESS)

return,

%Prediction error sum of squares (PRESS) by using hat matrix
%[r c] = size(D);
%n = r; %number of data
%Y = D(:,c); %response vector
%X = [ones(n,1) D(:,1:c-1)]; %design matrix
%[Q,R] = qr(X);
%b = R\(Q'*Y); %parameters estimation by using a pivoted QR 
                  %factorization of X
%Ye = X*b; %expected response value
%e = Y-Ye; %residual term
%H = X*inv(R'*R)*X'; %hat matrix
%hii = diag(H); %leverage of the i-th observation
%x = sum((e./(1-hii)).^2); %prediction error sum of squares (PRESS)