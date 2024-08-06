function [out] = MyMax(Matrix)

[N, M] = size(Matrix);

IndexMatrix = zeros(N, M);

for n = 2:(N-1)
    for m = 1:M
        A      = Matrix(n, m);
        Aup    = Matrix(n-1, m);
        Adown  = Matrix(n+1, m);
        if (m==1) Aleft  = 0; else Aleft   = Matrix(n, m-1); end
        if (m==M) Aright = 0; else Aright  = Matrix(n, m+1); end
        if (A>Aup)&&(A>Adown)&&(A>Aleft)&&(A>Aright)
            IndexMatrix(n, m) = 1;
        else
            IndexMatrix(n, m) = 0;
        end
    end
end

if sum(sum(IndexMatrix))==0
    out = min(min(Matrix));
else
    out    = max(Matrix(logical(IndexMatrix)));
end