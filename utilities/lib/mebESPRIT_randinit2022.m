% Simplified MEB-ESPRIT core function 
%
% function est = mebESPRIT_randinit2022(R, M0, C, Qa, arr)
%
% - R: covariance cell matrix of size (num_array x num_array)
%     Each cell contains a 2D matrix of size 
%     (num_spherical_harmonics x num_spherical_harmonics) 
%
% - M0, C{1}, C{2}, C{3}: EB-ESPRIT coefficients 
% - Qa: assumed number of sources 
% - arr: array structure (only array positions in arr.xyz are used here)
% 
% For the construction of coefficients, please refer to:
% J. -W. Choi, F. Zotter, B. Jo and J. -H. Yoo, 
% "Multiarray Eigenbeam-ESPRIT for 3D Sound Source Localization 
% With Multiple Spherical Microphone Arrays," 
% in IEEE/ACM Transactions on Audio, Speech, and Language Processing, 
% vol. 30, pp. 2310-2325, 2022, doi: 10.1109/TASLP.2022.3183929
%
% Author: Jung-Woo Choi (KAIST) and Franz Zotter (IEM, Graz)
% jwoo@kaist.ac.kr, zotter@iem.at


function est = mebESPRIT_randinit2022(R, M0, C, Qa, arr)
% Perform mebESPRIT
  run common.m

   [L1, L] = size(M0); 
   narr = size(R,1);
   nxyz = 3; %x,y,z
   
   RT = cell2mat(R);       % convert to 2D mtx (narr*nm, narr*nm)
   [UT,~] = eigs(RT,Qa);   % total signal subspace (narr*nm, Qa)
    
   for ia=1:narr % for each array
       Ul{ia} = UT((1:L)+(ia-1)*L, :); % partitioning subspace for each array
       MU{ia} = M0*Ul{ia}; %{L}(H1,Qa)
             
       for ixyz=1:nxyz
           % RHS matrix of all arrays
           CU{ia}(:,:,ixyz) = C{ixyz}*Ul{ia}; %{L}(H1,Qa,3) 
       end
   end

    %% %%%%%%%%%%%%%%%%%%%%%%% jointevd
    [ma, pos, prob] = jointevd(CU, MU, 100, arr); % merged with Zotter's gjsd
    E = real(ma); % (Qa, xyz, narr)

    %% %%%%%%%%%%%%%%%%%%%%%%  Build est structure  
    est.ev  = E;                         % (Qa, xyz, narr)
    est.evnorm = reshape(vecnorm(E,2,2),Qa,narr); % calc unitnorm
    est.doa = xyz2phith(E);
    est.xyz = pos; % (Qa, xyz)
    est.prob = prob;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Joint EVD (iter)
function [ma, pos, b] = jointevd(A, B, maxit, arr)
bdebug = 0;

A0=A;
B0=B;

tol = 1e-7;
tol2= 1e-1; %1e-4;
reg = 1e-6;
maxreinit=4;

%- algorithmic settings
l = numel(A);
[n,m,r] = size(A{1});
n0 = n;
m0 = m;
r0 = r;

%- initialize
ma = zeros(m,r,l);
b = zeros(m,l);
pos = zeros(m,3);
q = 1; 
reinit=0;

while (m>1) % for each source

    x = randn(m,1);   % ones(m,1); % rand init 

    x = x(:)/norm(x); % normalize right eigenvec
    mm = zeros(r,l);  % randn(r,l);
    y  = zeros(n,l);  % left eigenvec
    it = 1;
    while (it<maxit)
        xold = x; 
        %1 left eigenvector x
        for li=1:l
            y(:,li) = B{li} * x;       
            b(q,li) = norm(y(:,li));
            
        end
        
        %2 eigenvalues
        for li = 1:l
            for ri=1:r
                mm(ri,li) = real(y(:,li)'*A{li}(:,:,ri)*x); 
            end
            mm(:,li)=mm(:,li)/norm(mm(:,li));
        end
        [mm, pos_] = pos3d(mm,arr.xyz,b(q,:));

        ma(q,:,:)=mm;
        pos(q,:) = pos_(:)';
        
        %3 right eigenvector
        Py = zeros(n,n,l);
        for li=1:l
            Py(:,:,li) = eye(n) - y(:,li)*y(:,li)'/(b(q,li))^2;
        end
        Pm = zeros(r,r,l);
        for li=1:l
            Pm(:,:,li) = eye(r) - mm(:,li)*mm(:,li)'; %*norm(mm(:,li))^2;
        end
        C = reg*eye(m,m);
        
        for li=1:l
            for ri=1:r
                for ri2=1:r
                    C = C + A{li}(:,:,ri)'*Py(:,:,li)*A{li}(:,:,ri2)*Pm(ri,ri2,li); 
                end
            end
        end
        x = C \ x;
        x = x/norm(x);
        
        %4 converged?
        e = 1-abs(xold'*x);
       
        if e<=tol
           break
        end % end if e<=tol
        it = it +1;
    end % end while it

    x = xold;     

    Ah=zeros(n,m,r,l);
    Bh=zeros(n,m,l);
    
    [ux, betax]= gallery('house',x(:));
    
    for li=1:l
        [uy, betay]= gallery('house',y(:,li));
        Bh(:,:,li) = B{li} - betay*uy*(uy'*B{li}) - betax*(B{li}*ux)*ux' + betax*betay*(uy'*B{li}*ux)*(uy*ux');
        for ri=1:r
            Atmp = A{li}(:,:,ri); 
            Ah(:,:,ri,li)= Atmp - betay*uy*(uy'*Atmp) - betax*(Atmp*ux)*ux' + betax*betay*(uy'*Atmp*ux)*(uy*ux');
            mc(ri,li) = squeeze(Ah(1,1,ri,li))./squeeze(Bh(1,1,li));
        end
       
    end
    
    %- valid estimation?
    if (norm(squeeze(Ah(2:end,1,:)),'Fro')+...
            norm(squeeze(Bh(2:end,1,:)),'Fro')< ...
        (norm(squeeze(Ah(1,1,:)),'Fro')+norm(squeeze(Bh(1,1,:)),'Fro'))*tol2)||(reinit>=maxreinit) 
        if bdebug
            fprintf('Ah for q = %d \n', q);
            fprintf([repmat([repmat(' %2.3f',1,m) ' |'], 1, r*l) '\n'],reshape(abs(Ah),size(Ah,1),[])');
            fprintf('Bh for q = %d \n', q);
            fprintf([repmat([repmat(' %2.3f',1,m) ' |'], 1, l) '\n'],reshape(abs(Bh),size(Bh,1),[])');
            fprintf('lambda (mc) for q = %d \n', q);
            fprintf([repmat([repmat(' %2.3f',1,r) ' |'], 1, l) '\n'],mc.');
            fprintf('lambda (mm) for q = %d \n', q);
            fprintf([repmat([repmat(' %2.3f',1,r) ' |'], 1, l) '\n'],mm.');
            fprintf('norm(y) for q = %d \n', q);
            fprintf([repmat(' %2.3f |', 1, l) '\n'],b(q,:));
        end
        m=m-1;
        n=n-1;
        q=q+1;
        
        for li=1:l
            A{li}=reshape(Ah(2:end,2:end,:,li),[n,m,r]);
            B{li}=reshape(Bh(2:end,2:end,li),[n,m]);
        end

        reinit=0;
    else
     
        reinit=reinit+1;
%         disp(['reinit:' num2str(q) '_' num2str(reinit)])
    end
end % end while m>1

% last source: matrices become vectors
mm = zeros(r,l);
for li=1:l
    x = 1; % right eigenvector is a unit scalar
    y = B{li}(:,1)*x;
    b(q,li) = norm(y);
    y = y / b(q,li);
    mm(:,li) = real(y'*reshape(A{li}(:,1,:),size(A{1},[1 3]))*x);
    mm(:,li)= mm(:,li)/norm(mm(:,li));
end
[ma(q,:,:), pos_] = pos3d(mm,arr.xyz,b(q,:));
pos(q,:) = pos_(:)';

% 
if bdebug
    fprintf('lambda for q = %d \n', q);
    fprintf([repmat([repmat(' %2.3f',1,r) ' |'], 1, l) '\n'],mm.');
end
%     
end % end jointevd

function H=house(v) % complex version was worse

[u, beta]= gallery('house',v(:));
H = eye(numel(v)) - beta*(u*u');

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D position estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (for each source)
function [evout, Pos] = pos3d(ev, arrxyz, b) 
    
    % ev: xyz, arr
    % arrxyz: arr, xyz
    % Pos: 3,1 (xyz)
    
    axyz = arrxyz'; % xyz, arr

    w = b/rms(b);
    w = w.^(1/3);
    narr = size(axyz,2);
    
    E = zeros(narr*3,3); D=zeros(narr*3,1);
    for ia=1:narr
        E((ia-1)*3+(1:3),:) = w(ia)*(eye(3) - ev(:,ia)*ev(:,ia)');
        D((ia-1)*3+(1:3)) = E((ia-1)*3+(1:3),:)*axyz(:,ia); 
    end

    Pos = (E'*E)\E'*D; %E\D; % (xyz, 1) % global position (pseudo inv is important)
    apos = Pos - axyz; % xyz, narr (relative position from each array)
    evout = apos./vecnorm(apos,2,1);
    
    
end % end function pos3D

