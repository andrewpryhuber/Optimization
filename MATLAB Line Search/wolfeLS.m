
function [xi, fcalls] = wolfeLS(xc,d,fc,Dfc,DDfnc,c1, c2, alpham,eps)

% Algorithm 3.2 in Nocedal and Wright
%Line search algorithm satisfying strong Wolfe conditions
%The subroutine requires the following input:
%       xc    = the current point,
%       d     = the direction of search,
%       fc    = the current function value,
%       fnc   = a string variable for the function name,
%       Dfc   = derivative of fnc at xc
%       DDfnc = the directional derivative of fnc at xc in the
%               direction d, must have  DDfnc < 0,
%       c1     = the slope modification parameter in (0,1),
%       c2     = the slope modification parameter in (0,1),
%       alpham = the maximum alpha parameter in (0,1),
%

if DDfnc >= 0
    error('Not a direction of descent. Program has been terminated.')
end
if c1<= 0 || c1>= 1
    error('The parameter c1 in the backtracking subroutine is not in (0,1).')
end
if c2<=0 || c2 >=1
    error('The parameter c2 in the backtracking subroutine is not in (0,1).')
end

alpha0 = 0; %write minimum alpha
alpha = 1;
alphai = alpha;

fi_1 = fc;
fcalls=1;

while true
    
    xi = xc + alphai*d;
    fi = objEval(xi);
    gi = gradEval(xi);
    DDfni = gi'*d;
    
    %check to see if the current iterate sufficiently decreases
    if (fi > fc + c1*alphai*DDfnc) || ((fcalls > 1) && (fi >= fi_1))
        [alpha, numcalls] = zoom(xc,d,fc,fi,Dfc,DDfnc, c1, c2, alpha0,alphai);
        xi = xc + alpha*d;
        fcalls=fcalls+numcalls;
        return;
    end
    %current iterate suffciently decreases, so check if it is too close
    %if abs(DDfni) <= -c2*DDfnc %change to this for strong wolfe
    if DDfni>=c2*DDfnc
        alpha = alphai;
        xi = xc + alpha*d;
        return;
    end
    %behind the minimum
    if DDfni >= 0
        [alpha, numcalls] = zoom(xc, d, fc,fi,Dfc,DDfnc,c1,c2,alphai,alpha0);
        xi = xc + alpha*d;
        fcalls=fcalls+numcalls;
        return;
    end
    
    alpha0 = alphai;
    alphai = min(alpham, alphai*3);
    fi_1 = fi;
    fcalls = fcalls+1;
    if norm(d) <= eps
        disp('linesearch step too small')
        if fi >= fc
            xi  =  xc;
        end
        break
    end
end
end


% zoom function as described in Algorithm 3.3 in Nocedal and Wright
% returns parameter alpha and number of internal function calls

function [alpha, count] = zoom(xc,d,fc,f_lo,~,DDfnc,c1, c2, alphaLo,alphaHi)
count=0;
alphap=0;
alpha =1 ;
eps = 1e-8;
while true || abs(alphap-alpha)>eps
    alpha = 1/2*(alphaLo+alphaHi);
    xj = xc + alpha*d;
    fj = objEval(xj);
    Dfcj = gradEval(xj);
    count=count+1; %update count
    
    %if sufficient decrease in point alpha, then set maximum of the
    %feasible interval to alpha
    if ((fj > fc + c1*alpha*DDfnc) || (fj >= f_lo))
        alphaHi = alpha;
    else
        %what if you fulfulled strong wolfe?
        if Dfcj'*d>=c2*DDfnc
            return;
        end
        if Dfcj'*d*(alphaHi-alphaLo) >= 0 %positive slope and alphaHi>alphaLo
            alphaHi = alphaLo;
            alphaLo = alpha;
            f_lo =fj;
        end
        
    end
    
    if alphap==alpha
        return
    end
    alphap=alpha;
    
end
end
