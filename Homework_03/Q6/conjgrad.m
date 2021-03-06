function [x, rsnew] = conjgrad(A, b, x)
    r = b - 2*A * x;
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = 2*A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end