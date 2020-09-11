function [R]=findRotation(A,B,varargin)
%FindRotation is adopted from absor.m by Matt Jacobson [*]
%and finds the rotation that best maps one collection of 3d point coordinates to
%another, both with centroids at origo. It is based on Horn's quaternion-based method. 

%DESCRIPTION:
%
%As input data, one has
%
%  A: 3xN matrix whos columns are the coordinates of N source points.
%  B: 3xN matrix whos columns are the coordinates of N target points.
%
%The basic syntax
%
%     [regParams,Bfit,ErrorStats]=absor(A,B)
%
%solves the unweighted/unscaled registration problem
%
%           min. sum_i ||R*A(:,i) + t - B(:,i)||^2
%
%for unknown rotation matrix R and unknown translation vector t.

%OUTPUTS:
% R:   The estimated rotation matrix, R

% [*] Matt J (2020). Absolute Orientation - Horn's method
% (https://www.mathworks.com/matlabcentral/fileexchange/26186-absolute-orientation-horn-s-method), MATLAB Central File Exchange. Retrieved May 27, 2020.


dimension=size(A,1);

if dimension~=size(B,1)
    error 'The number of points to be registered must be the same'
end


%%Centering
left=A; right=B;

M=left*right.';


%%Compute rotation matrix
Mflat = (M(:));
Sxx = Mflat(1);
Syx = Mflat(2);
Szx = Mflat(3);

Sxy = Mflat(4);
Syy = Mflat(5);
Szy = Mflat(6);

Sxz = Mflat(7);
Syz = Mflat(8);
Szz = Mflat(9);

N=[(Sxx+Syy+Szz)  (Syz-Szy)      (Szx-Sxz)      (Sxy-Syx);...
    (Syz-Szy)      (Sxx-Syy-Szz)  (Sxy+Syx)      (Szx+Sxz);...
    (Szx-Sxz)      (Sxy+Syx)     (-Sxx+Syy-Szz)  (Syz+Szy);...
    (Sxy-Syx)      (Szx+Sxz)      (Syz+Szy)      (-Sxx-Syy+Szz)];

[V,D]=eig(N);

[~,emax]=max(real(  diag(D)  )); emax=emax(1);

q=V(:,emax); %Gets eigenvector corresponding to maximum eigenvalue
q=real(q);   %Get rid of imaginary part caused by numerical error

[~,ii]=max(abs(q)); sgn=sign(q(ii(1)));
q=q*sgn; %Sign ambiguity

%map to orthogonal matrix

quat=q(:);
nrm=norm(quat);
if ~nrm
    disp 'Quaternion distribution is 0'
end

quat=quat./norm(quat);

q0=quat(1);
qx=quat(2);
qy=quat(3);
qz=quat(4);
v =quat(2:4);


Z=[q0 -qz qy;...
    qz q0 -qx;...
    -qy qx  q0 ];

R=v*v.' + Z^2;
