function [xvec, muvec, gvec, g, res,flag] = gen_compute_eqm_correct(lambda,pivec,param,kappa,r,xinit,c)
% See a description of this function in the comment beginning on Line 30 below.    
    n = (length(pivec)-1)/2;
    if isempty(xinit)
        xinit = zeros(1,2*n); xinit(n+1) = 1;
    end
    tol=1e-12;
    %opt = optimoptions('fsolve','Display','off','FunctionTolerance',tol,'OptimalityTolerance',tol,'MaxFunctionEvaluations',50000,'MaxIterations',1000);
    opt = optimoptions('fsolve','Display','off','FunctionTolerance',tol,'OptimalityTolerance',tol);
    [xvec,~,flag] = fsolve(@(xvec)gen_eqm_eqns_correct(xvec,pivec,param,kappa,r,c),xinit,opt);xvec=abs(xvec);
    res = sum(abs(gen_eqm_eqns_correct(xvec,pivec,param,kappa,r,c)));
    %xvec=xvec/r;
    muvec = zeros(1,n+1);
    muvec(1) = 1; muvec(2) = 2*xvec(n+1) / (xvec(n)+kappa);
    xtmp = [xvec,0];
    for i=2:n
        muvec(i+1) = muvec(i)* xtmp(n+i) / (xtmp(n+1-i)+kappa);
    end
    muvec = muvec./sum(muvec);
    
    [g,gvec]=gen_compute_g(muvec,xvec,lambda,kappa);
    if 1==0
    gvec(1) = 2*muvec(1)*xvec(n+1);
    gvec(2:n+1) = muvec(2:n+1).*xtmp(n+2:(2*n+1));
    gvec = gvec.*log(lambda);
    g = sum(gvec);
    end
end

% This function (gen_compute_eqm_correct.m) is used to obtain the red dashed line in
% in our Figure 3. 
%
% This function is identical to gen_compute_eqm.m in LMS's replication
% code, except that:
%
% (i) Line 12 is commented out.  Commenting out Line 12 does
% not change the equations solved or the interpretation of the function's
% output, because Line 3 of gen_eqm_eqns_correct is also commented out.
%
% (ii) Line 10 calls gen_eqm_eqns_correct.m, not LMS's gen_eqm_eqns.m
% gen_eqm_eqns_correct.m solves the model using the investment cost functional form
% provided in LMS's paper: (c^2 eta_s)^2.  LMS state (p. 213) that, in their
% quantitative model, the cost of achieving investment success rate eta is c(eta_s)=(c * eta_s)^2.
% This functional form assumption is used here.  
%
% When calling this function to produce the red line in our Figure 3,
% we choose a value for c (c = 100/sqrt(2)/33.3569 )
% such that the BGP investment success rates match those in LMS's code.  
% (This value for c is different from the value reported in LMS Table 1, which reports 
% c = 33.4.  Using c = 100/sqrt(2)/33.3569 and the functional form assumption for
% investment costs in their paper reproduces LMS's BGP investment success rates.
% Thus, their code effectively imposes a much lower investment cost than in the model
% as described in their paper.  When the model described in their paper is solved
% correctly with c = 33.4, there is no inverted-U relation between productivity growth and the interest
% rate. In this economy, productivity growth and the profit share are
% very low (as the discount rate falls from 6% toward 0, productivity growth rises from
% about 0.01 percent to about 0.03 percent, or from 1 to 3 basis points), because the investment cost
% in their paper is (c * eta_s)^2 and a value of c = 33.4 corresponds to a very high investment cost.)
%
% Overview of investment cost assumptions in LMS's paper and code:
% LMS state (p. 213) that, in their quantitative model, the cost of achieving
% investment success rate eta_s is c(eta_s) = (c * eta_s)^2, with c = 33.4 (their Table 1).
% Their code instead solves for investment success rates and productivity growth
% as if the investment cost were (100 / sqrt(2) / 33.4 * eta_s)^2.  Their code 
% calculates profit shares as if the investment cost were (100 / 33.4 * eta_s)^2, 
% or twice as large. 
