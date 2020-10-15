% Driver file for local minimization of a real valued function
% using strong Wolfe line search (Algorithm 3.2 in Nocedal and Wright).
%
% Packages Required:
%   1) Matlab Symbolic Toolbox
%
% User must input:
%   1) Symbolic variables
%   2) Objective function
%   3) Initial iterate
%
% Output:
%   1) Final iterate and function value after algorithm terminates.
%   2) Plots of function calls, step sizes, norm of gradient, and function
%   values at each iteration.
   

function main

%Line search parameters for (Strong) Wolfe LS method (Adapted from Alg 3.2 in
%Chapter 3 of Nocedal and Wright)
c1 = 1e-5;
c2 = 0.5;
alpham = 100;

%Termination parameters
eps     =  1.0e-10; %epilson for measuring size of step
epsf    =  1.0e-12; %epsilon for measuring norm of gradient
maxit   =  100000;
iter    =  0;

%%%%%%%%%%%%%
%% Test Case 1
% [vars,obj,x0] = functionFromFile('test1.txt');
%%%%%%%%%%%%%

%%%%%%%%%%%%%
%% Test Case 2
% [vars,obj,x0] = functionFromFile('test2.txt');
%%%%%%%%%%%%%

%%%%%%%%%%%%%%
%% Test Case 3
% generate random convex quadratic
NUM_VARS = 10;
syms vars [NUM_VARS 1] real;
sqrtQ = rand(NUM_VARS,NUM_VARS);
g =  rand(NUM_VARS,1);
Q = sqrtQ*sqrtQ.';
obj = vars'*Q*vars + vars'*g;
x0 = rand(NUM_VARS,1); % initial iterate
%%%%%%%%%%%%%%

% write m files to evaluate objective and its gradient
matlabFunction(obj,'Vars',{vars},'File','objEval');
matlabFunction(jacobian(obj,vars)','Vars',{vars},'File','gradEval');

%Initialization
xc   =  x0;
fc =  objEval(xc);
Dfc = gradEval(xc);
nDfc     =  norm(Dfc);

%Initialize main loop.
ndiff   =  1;
data    =  [iter,nDfc,ndiff,1,fc];

%Are we already at a solution, or should we continue.
if nDfc <= epsf %norm of gradient
    disp('Termination due to small gradient.')
else
    
    
    %The main loop.
    %while norm of gradient is greater than epsf, norm of difference is greater
    %than eps, and iter is less than max iter
    while nDfc > epsf && ndiff > eps && iter < maxit
        if iter==0
            theta = nDfc;
            d=-Dfc/theta;
        else
            s=xn-xc;
            y = Dfcn-Dfc;
            
            if (s'*y)>0
                theta = (y'*y)/(y'*s);
            else
                theta = nDfc;
            end
            d=-Dfcn/theta;
        end
        
        if iter>0
            xc=xn; %set new xn equal to old xc
            Dfc=Dfcn;
        end
        
        DDfnc   =  Dfc'*d; % compute directional derivative
        % run weak Wolfe line search
        [xn, fcall] = wolfeLS(xc,d,fc,Dfc,DDfnc,c1,c2,alpham, eps);
        
        ndiff   =  norm(xn-xc); %size of step
        fc      = objEval(xn);
        Dfcn    = gradEval(xn);
        nDfc    =  norm(Dfcn);
        iter    =  iter + 1;
        data    =  [data;[iter,nDfc,ndiff,fcall,fc]];
        
        %Report reason for termination.
        if nDfc <= epsf
            disp('Termination due to small gradient.')
            break
        elseif ndiff <= eps
            disp('Termination due to small steps.')
            break
        elseif iter == maxit
            disp('Exceeded maximum number of iterations.')
            break
        end
    end
end






%Report some results
[datalength,~]   =  size(data);
onenorm              =  ones(datalength,1);
Tfcall               =  onenorm'*data(:,4);
disp('Final point is '), xc
disp('Final function value is'), fc
disp('Final norm of the gradient is '), nDfc
disp('Total number of function calls is '), Tfcall
disp('Total number of gradient calls is '), iter

figure;
subplot(2,2,1), plot(data(:,4))
xlabel('iteration'),ylabel('function calls')
subplot(2,2,2), semilogy(data(:,1),data(:,2))
xlabel('iteration'),ylabel('norm of the gradient')
subplot(2,2,3), semilogy(data(:,1),data(:,3))
xlabel('iteration'),ylabel('distance between iterates')
subplot(2,2,4), semilogy(data(:,1),data(:,5))
xlabel('iteration'),ylabel('function value')


end




