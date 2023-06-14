function [Mx,My,Mz] = multiplication_theorems_real(N)
% multiplication theorems for real-valued SHs of maximum order N
% Franz Zotter and Thomas Deppisch, 2020

    Mx=zeros((N)^2,(N+1)^2);
    My=zeros((N)^2,(N+1)^2);
    Mz=zeros((N)^2,(N+1)^2);
    
    w = @(n,m) sqrt((n+m-1)*(n+m)/((2*n-1)*(2*n+1)));
    v = @(n,m) sqrt((n-m)*(n+m)/((2*n-1)*(2*n+1)));
    s = @(n,m) sqrt(2-(m==0))*sqrt(2-(n==0));
    sgn = @(m) (2*(m>=0)-1);
    nm = @(n,m) n.^2+n+m+1;
    
    for n=1:N-1
        for m=[-n+2:-1, 1:n] %anm, a(n,0)=0
            Mx(nm(n, m),nm(n-1,m-1)) = -sgn(-m)*w(n  , m  )/s(m,m-1);
            % anm
            My(nm(n,-m),nm(n-1,m-1)) = -sgn(-m)*w(n  , m  )/s(m,m-1);
        end
        if n>1
            m=0;
            % anm
            My(nm(n,-m),nm(n-1,m-1)) = -sgn(-m)*w(n  , m  )/s(m,m-1);
        end
    end
    for n=0:N-1
        for m=[-n:-1,1:n] %bnm, b(n,0)=0
            Mx(nm(n, m),nm(n+1,m-1)) = sgn(-m)*w(n+1,-m+1)/s(m,m-1);
            % bnm
            My(nm(n,-m),nm(n+1,m-1)) = sgn(-m)*w(n+1,-m+1)/s(m,m-1);
        end
        m=0; %bnm
        My(nm(n,-m),nm(n+1,m-1)) = sgn(-m)*w(n+1,-m+1)/s(m,m-1);    
    end
    for n=1:N-1
        for m=[-n:-2, 1:n-2] %cnm, c(n,-1)=0
            Mx(nm(n, m),nm(n-1,m+1)) = -sgn(m)*w(n  ,-m  )/s(m,m+1);
            % bnm, c(n,-1)=0, c(n,0)=0
            My(nm(n,-m),nm(n-1,m+1)) =  sgn(m)*w(n  ,-m  )/s(m,m+1);
        end
        m=0; % cnm
        Mx(nm(n, m),nm(n-1,m+1)) = -sgn(m)*w(n  ,-m  )/s(m,m+1);
    end

    for n=0:N-1
        for m=[-n:-2, 1:n] %dnm, d(n,-1)=0
            Mx(nm(n, m),nm(n+1,m+1)) =  sgn(m)*w(n+1, m+1)/s(m+1,m);
            %dnm, d(n,-1)=d(n,0)=0
            My(nm(n,-m),nm(n+1,m+1)) = -sgn(m)*w(n+1, m+1)/s(m+1,m);
        end
        m=0; %dnm
        Mx(nm(n, m),nm(n+1,m+1)) =  sgn(m)*w(n+1, m+1)/s(m+1,m);
    end
    
    for n=1:N-1
        for m=-n+1:n-1
            Mz(nm(n,m),nm(n-1,m  )) = v(n,m);
        end
    end
    for n=0:N-1
        for m=-n:n
            Mz(nm(n,m),nm(n+1,m  )) = v(n+1,m);
        end
    end
end

