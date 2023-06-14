% common.m containing anonymous functions 

vect = @(x) x(:); % reshape to vector
unvect = @(x,A) reshape(x,size(A)); % unvectorize to size of A
dB10 = @(x) 10*log10(abs(x)); % 10 dB
dB20 = @(x) 2*dB10(x); % 20 dB
sq   = @(x) squeeze(x);

% matlab elevation to zenith
aziElev2aziZenith = @(azel) [azel(:,1,:) pi/2-azel(:,2,:)];

%- xyz to ev (ONLY for complex)
xyz2ev = @(xyz) [xyz(:,1,:) + 1i*xyz(:,2,:), xyz(:,1,:) - 1i*xyz(:,2,:), xyz(:,3,:)]./ ...
                   sqrt(sum(abs(xyz).^2,2)); % nsrc, 3
ev2xyz = @(ev) real([ (ev(:,1) + ev(:,2))/2, (ev(:,1) - ev(:,2))/2i, ev(:,3)]); % nsrc,3

%- xyz to phith
xyz2phith = @(xyz) [atan2(xyz(:,2,:),xyz(:,1,:)), atan2(sqrt(xyz(:,1,:).^2+xyz(:,2,:).^2),xyz(:,3,:))]; % (phi,theta),3
phith2xyz = @(phith) [cos(phith(:,1,:)).*sin(phith(:,2,:)) sin(phith(:,1,:)).*sin(phith(:,2,:)) cos(phith(:,2,:))];

%- xyz to spherical 
xyz2sph = @(xyz) [sqrt(sum(abs(xyz).^2,2)) atan2(xyz(:,2,:),xyz(:,1,:)), atan2(sqrt(xyz(:,1,:).^2+xyz(:,2,:).^2),xyz(:,3,:))]; % (r, phi,theta),3
sph2xyz = @(sph) rphith(:,1,:).*[cos(rphith(:,2,:)).*sin(rphith(:,3,:)) sin(rphith(:,2,:)).*sin(rphith(:,3,:)) cos(rphith(:,3,:))];

%- angular distance
adist = @(phith1,phith2) acos( cos(phith1(:,2,:)).*cos(phith2(:,2,:)) + sin(phith1(:,2,:)).*sin(phith2(:,2,:))... 
                        .*cos(phith1(:,1,:)-phith2(:,1,:))); %%elementwise
                    
adistm = @(phith1,phith2) acos( mtimesx(cos(phith1(:,2,:)),cos(phith2(:,2,:)),'t') + ...
                        mtimesx(sin(phith1(:,2,:)),sin(phith2(:,2,:)),'t')... 
                        .* cos(phith1(:,1,:)- permute(phith2(:,1,:),[2 1 3])) ); %% all 
                    
%- position distance
distm = @(xyz1,xyz2) vecnorm( permute(xyz1,[1 3 2]) - permute(xyz2,[3 1 2]) , 2, 3);

%- off diagonal
offdiag = @(X) X - diag(diag(X));
