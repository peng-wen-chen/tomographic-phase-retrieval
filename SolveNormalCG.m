function [x, err] = SolveNormalCG( b, x,input) %p0,T_hat,n)
    tol=1e-20;
    
    Ax= CG_LHS(x,input.T,input.p,input.n);
    r = b - Ax;
    p = r;
    rsold = r(:)' * r(:);
    for ii = 1: input.CGiter
    Ap= CG_LHS(p,input.T, input.p,input.n);
         alpha = rsold / (p(:)' * Ap(:));
        x = x + alpha * p;
        if (rsold)< tol
            break
        end
        r = r - alpha * Ap;
        rsnew = r(:)' * r(:);
        if sqrt(rsnew) < tol
            
              break
        end
         p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
    err= sqrt(rsnew)/norm(b(:));
end