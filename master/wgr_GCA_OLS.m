function [cgc] = wgr_GCA_OLS(y,x,z,order)
%% y->x condtion on z. based on covariance matrix
%% x: N*nx, y: N*ny, z:N*nz.
%% 
[N,nx]=size(x);
%now
X = x(order+1:end,:);
%past
past_ind = repmat([1:order],N-order,1) + repmat([0:N-order-1]',1,order);
xz = [x z];
XZ_past = reshape(xz(past_ind,:),N-order,order*size(xz,2));
Y_past = reshape(y(past_ind,:),N-order,order*size(y,2));
XZY_past = [XZ_past Y_past];

% Remove mean
xzyc = bsxfun(@minus,XZY_past,sum(XZY_past,1)/(N-order)); 
Xc = bsxfun(@minus,X,sum(X,1)/(N-order)); 
xzc = bsxfun(@minus,XZ_past,sum(XZ_past,1)/(N-order)); 

%Covariance matrix
cov_X = (Xc' * Xc) / (N-order-1); 
cov_xz = (xzc' * xzc) / (N-order-1);
cov_xzy = (xzyc' * xzyc) / (N-order-1);
cov_X_xz = (Xc' * xzc) / (N-order-1); 
cov_X_xzy = (Xc' * xzyc) / (N-order-1);

%Partial cross-covariance
cov_X_xz = cov_X - cov_X_xz/cov_xz*cov_X_xz';
cov_X_xyz = cov_X - cov_X_xzy/cov_xzy*cov_X_xzy';
cgc = log(det(cov_X_xz)/det(cov_X_xyz));
end