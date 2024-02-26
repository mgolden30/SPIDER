%{
Generate a feature matrix G for vectors using simulated QED data
%}
%% Load data
clear;
addpath("../SPIDER_functions/");
%load in vorticity, flow velocity, ....
%load('../../temp.mat')
load('fullData1.mat');
%dt=t(2)-t(1);
%load("InitialCondition1.mat";

o=w;
N = size(o, 1);
dx = 2*pi/N;
omega = o;
wL    = o;
%uL    = zeros( [size(o),2] );
%[Us, Vs, Ws, Ul, Vl, Wl] = separate_scales( omega );
%%return;
%Integrate library
%addpath("../SPIDER_functions/");
number_of_library_terms = 3;   %under-estimate this
number_of_windows       = 128; %number of domains we integrate over
degrees_of_freedom      = 1;   %vectors have one degree of freedom
dimension               = 3;   %how many dimensions does our data have?
envelope_power          = 4;   %weight is (1-x^2)^power
size_vec                = [ 32, 32, 32 ]; %how many gridpoints should we use per integration?
buffer                  = 3; %Don't use points this close to boundary
%define shorthand notation
nl = number_of_library_terms;
nw = number_of_windows;
dof= degrees_of_freedom;
%Make important objects for integration
pol      = envelope_pol(envelope_power, dimension);
G        = zeros( dof*nw, nl );
labels   = cell(nl, 1);
scales   = zeros(1,nl);
%we also need pol to be odd in time. multiply it by t
pol0 = pol; %save for normalization
pol = multiply_pol( pol, legendre_pol( [0,0,1], dimension ) );

U = o;
size_of_data = size(U, 1:dimension);
corners = pick_subdomains( size_of_data, size_vec, buffer, nw );
% Nondimensionalize!
% TODO
% Create grid variable
x = (1:size(U,1))*dx;
y = (1:size(U,2))*dx;
grid = {x,y,t};


%length and time scales
L = x(end) - x(1); %Use size of domain
T = t(end) - t(1); %Use size of domain
%{
clear;
N = 128;
[x,y] = meshgrid( (0:N-1)/N*2*pi );
w = exp( -cos(x).^2 ./(1 + 0.5*sin(3*y)) );
 
w2 = zeros(N,N,2);
w2(:,:,1) = w;
w2(:,:,2) = w;
imagesc(w)
colormap jet
%}

k = -N/2:N/2-1;
k = fftshift(k);
%If you need a third axis for time
kx = reshape(k, [1,N,1]);
ky = reshape(k, [N,1,1]);
N = size(w, 1);
dx = 2*pi/N;
omega = w;
U = 0*omega;
V = U;
omega_x = U;
omega_y = U;
for t = 1:size(omega,3)
  k = 0:N-1;
  k(k>N/2) = k(k>N/2) - N;
  k_sq = k.^2 + k'.^2;
  k_sq(1,1) = 1;
  o = fft2( omega(:,:,t) );
  V(:,:,t) = -real(ifft2( o.*1i.*k./k_sq ));
  U(:,:,t) =  real(ifft2( k'.*o.*1i./k_sq ));
  omega_x(:,:,t) =  real(ifft2( o.*1i.*k ));
  omega_y(:,:,t) =  real(ifft2( k'.*o.*1i ));
end
  
dU_dx = ifft2( (1i*k)  .* fft2(U), "symmetric" );
dV_dx = ifft2( (1i*k)  .* fft2(V), "symmetric" );
dU_dy = ifft2(  1i*k'  .* fft2(U), "symmetric" );
dV_dy = ifft2(  1i*k'  .* fft2(V), "symmetric" );
deriv = @(w,nx,ny) ifft2( (1i*kx).^nx .* (1i*ky).^ny .* fft2(w), "symmetric" );

%dw_dx = deriv(w2,1,1);
%imagesc(squeeze( dw_dx(:,:,1)))


%% END OF BOILERPLATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO TOUCH PAST HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h=t(2)-t(1);
% dtWl=(S(:,:,2:end)-S(:,:,1:end-1))/h;
%dtuL=(uL(:,:,:,dof)-uL(:,:,:,j))/h; %what goes inside parenthesis?
%{
dWldx=ifft2(Dx.*fft2(Wl(:,:,j)),'symmetric');
dWldy=ifft2(Dx.*fft2(Wl(:,:,j)),'symmetric');
Lap_Wl = domain.nu*ifft2(D2.*fft2(Wl(:,:,j)),'symmetric');
dWldx2=ifft2(Dx.*fft2(Wl),'symmetric');
dWldy2=ifft2(Dy.*fft2(Wl),'symmetric');
duLdx=ifft2(Dx.*fft2(uL(:,:,1:j)),'symmetric');
duLdy=ifft2(Dy.*fft2(uL(:,:,1,j)),'symmetric');
duLdx2=ifft2(Dx.*fft2(uL(:,:,1:j),'symmetric'));
duLdy2=ifft2(Dy.*fft2(uL(:,:,1:j),'symmetric'));
%second column is dim,
%}
S = o;
phi2 = o;
uL = o;
uS = o;
wL = o;
wS = o;
advS = o;
advL = o;
Ws = o;
%str = {"S", "S^2", "\phi2", "uL", "uS", "wL", "Ws", "advS", "adL"};
%list_of_terms  ={S,(S).^2,phi2,uL,uS,wL,wS,advS,advL};
%
% str = {"S"};
% list_of_terms = {S};
%
% list_of_vectors={wL,Ws};
% list_of_scalars={uL,uS};
%
% a = 1;
% for p=1:size(list_of_terms)
%     for j=1:3 %first partial derivatives
%         labels{a}= "\partial_" + j + str(p) %does p extract the current term
%         G(:,a)   = SPIDER_integrate( list_of_terms{p}, [j], grid, corners, size_vec, pol );
%         a = a+1;
%     end
%
%     for n=1:3 % all combinations of second partial derivatives
%         for l=1:3
%             labels{a}= "\partial_"+ n + "\partial_" + l +str(p);
%             G(:,a)    = SPIDER_integrate(list_of_terms{p}, [n,l], grid, corners, size_vec, pol );
%             a = a+1;
%         end
%     end
% end
%
% return
% for k=1:size(list_of_vectors)%product of every element with every other element
%     for x=1:size(list_of_vectors)
%        % labels{a}= str(list_of_vectors(k))+ "times" str(list_of_terms(x)) %idk how to fix this error
%         G(:,a)    = SPIDER_integrate(list_of_vectors{k}.*list_of_vectors{x}, [], grid, corners, size_vec, pol );
%     end
% end
a = 1; %running index over library space
%is Wl dx defined somewhere?, if not using fft func? omega or Wl?
%curl?

labels{a} = "u_i_dt";
G(:,a)    = SPIDER_integrate_curl( U, V, [3], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;
 
%{
labels{a} = "di dj u_j";
G(:,a)    =    SPIDER_integrate_curl(  U, 0*U, [1,1], grid, corners, size_vec, pol ) ...
           +   SPIDER_integrate_curl(  V,   U, [1,2], grid, corners, size_vec, pol )...
           +   SPIDER_integrate_curl(0*U,   V, [2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;
%} 

labels{a} = "dj^2 u_i";
G(:,a)    = SPIDER_integrate_curl(U, V, [1,1], grid, corners, size_vec, pol ) ... 
          + SPIDER_integrate_curl(U, V, [2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

labels{a} = "u_j dj u_i";
G(:,a)    = SPIDER_integrate_curl(U.*dU_dx + V.*dU_dy, U.*dV_dx + V.*dV_dy, [], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;

return;

labels{a} = "di dj^2 dk u_k";
G(:,a)    = SPIDER_integrate( U,[1,1,1,1], grid, corners, size_vec, pol ) ...
+ SPIDER_integrate( V,[2,2,2,2], grid, corners, size_vec, pol )+SPIDER_integrate( U,V,[1,2,2,2], grid, corners, size_vec, pol )...
+SPIDER_integrate( U,V,[1,1,2,2], grid, corners, size_vec, pol )+SPIDER_integrate( U,V,[1,2,2,21], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;
 
labels{a} = "dj^2 dk^2 u_i";
G(:,a)    = SPIDER_integrate(U,V, [1,1,1,1], grid, corners, size_vec, pol )+2*SPIDER_integrate(U,V, [1,1,2,2], grid, corners, size_vec, pol )...
     + SPIDER_integrate(U,V, [2,2,2,2], grid, corners, size_vec, pol );
scales(a) = 1;
a         = a+1;
labels{a} = "u_i * dj u_j";
G(:,a)    =(SPIDER_integrate((U+V).*(DU/dx+dV/dy), grid, corners, size_vec, pol )+SPIDER_integrate(V,2, grid, corners, size_vec, pol ));
scales(a) = 1
a         = a+1;
 
labels{a} = "u_j * di u_j";
G(:,a)    = SPIDER_integrate((U+V).*(dU/dx+dV/dx+dU/dy+dV/dy),1, grid, corners, size_vec, pol);
scales(a) = 1;
a         = a+1;
  
labels{a} = "u_j * dj u_i"; %%how this differs from the other data set
  G(:,a)    = SPIDER_integrate( U.*U,V.*V,1, grid, corners, size_vec, pol ).*SPIDER_integrate( U.*U,V.*V, 2,grid, corners, size_vec, pol )%(gradU^2)
  
  scales(a) = 1;
  a         = a+1;
%
%
% phi2 = 0*Wl;
% for t = 1:size(Wl,3)
%   phi2(:,:,t) = phi;
% end
% %
  labels{a} = "di u_j * dj dk u_k";
  G(:,a)    = (SPIDER_integrate((dU/dx+dV/dx+dU/dy+dV/dy).*deriv(U,2,0)+deriv(V,0,2)+deriv(U+V,1,1),1, grid, corners, size_vec, pol ))
%  scales(a) = 1;
%  a         = a+1;
%
% % {
% labels{a} = "di u_j * dk^2 u_j";
%  G(:,a)    = (SPIDER_integrate(U,V,1, grid, corners, size_vec, pol )+ SPIDER_integrate(U,V, 2, grid, corners, size_vec, pol )).*...
%      (SPIDER_integrate(U,V, [1,1], grid, corners, size_vec, pol )+SPIDER_integrate(U,V, [2,2], grid, corners, size_vec, pol ));
%  scales(a) = 1;
%  a         = a+1;
%
%  labels{a} = "dj u_i * dj dk u_k";
%  G(:,a)    = (SPIDER_integrate( U,V, 1, grid, corners, size_vec, pol )+SPIDER_integrate( U,V, 2, grid, corners, size_vec, pol )).*...
%      (SPIDER_integrate( U, [1,1], grid, corners, size_vec, pol)+SPIDER_integrate( U,V, [1,2], grid, corners, size_vec, pol )+SPIDER_integrate( V, [2,2], grid, corners, size_vec, pol ));
%  scales(a) = 1;
%  a         = a+1;
%
%  labels{a} = "dj u_i * dk^2 u_j";
%  G(:,a)    = (SPIDER_integrate( U,V,1, grid, corners, size_vec, pol )+ SPIDER_integrate( U,V,2, grid, corners, size_vec, pol )).*...
%       (SPIDER_integrate( U,V,[1,1], grid, corners, size_vec, pol )+SPIDER_integrate( U,V,[2,2], grid, corners, size_vec, pol ));
% scales(a) = 1;
%  a         = a+1;
%
%  labels{a} = "dj u_j * di dk u_k";
%  G(:,a)    = (SPIDER_integrate( U, 1, grid, corners, size_vec, pol )+SPIDER_integrate( V, 2, grid, corners, size_vec, pol )).*...
%      (SPIDER_integrate( U, [1,1], grid, corners, size_vec, pol )+SPIDER_integrate( U,V,[1,2], grid, corners, size_vec, pol )+...
%      SPIDER_integrate( V, [2,2], grid, corners, size_vec, pol ));
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "dj u_k * di dj u_k";
%  G(:,a)    = (SPIDER_integrate( U,V,1, grid, corners, size_vec, pol )+SPIDER_integrate( U,V,2, grid, corners, size_vec, pol)).*...
%      SPIDER_integrate( U,V,[1,1], grid, corners, size_vec, pol )+2.*SPIDER_integrate( U,V,[1,2], grid, corners, size_vec, pol )+...
%      SPIDER_integrate( U,V,[2,2], grid, corners, size_vec, pol );
%   scales(a) = 1;
%   a         = a+1;
%
%   labels{a} = "dj u_k * di dk u_j"; %%how is this different from previous
%   G(:,a)    = SPIDER_integrate( U, [], grid, corners, size_vec, pol );
%   scales(a) = 1;
%   a         = a+1;
%
%    labels{a} = "dj u_j * dk^2 u_i";
%    G(:,a)    = (SPIDER_integrate( U, 1, grid, corners, size_vec, pol )+ SPIDER_integrate( V, 2, grid, corners, size_vec, pol )).*...
%         SPIDER_integrate( U,V, [1,1], grid, corners, size_vec, pol )+ SPIDER_integrate( U,V, [2,2], grid, corners, size_vec, pol );
%      scales(a) = 1;
%  a         = a+1;
%
%   labels{a} = "dj u_k * dj dk u_i, u_i * dj^2 dk u_k";
%   G(:,a)    = (SPIDER_integrate( U,V,1, grid, corners, size_vec, pol )+SPIDER_integrate( U,V,2, [], grid, corners, size_vec, pol )).*...
%       (SPIDER_integrate( U,V,[1,1], grid, corners, size_vec, pol )+2.*SPIDER_integrate( U,V,[1,2], grid, corners, size_vec, pol )+SPIDER_integrate( U,V,[2,2], grid, corners, size_vec, pol )).*...
%       (SPIDER_integrate( U,[1,1,1], grid, corners, size_vec, pol )+SPIDER_integrate( V,[1,1,2], grid, corners, size_vec, pol )+...
%       SPIDER_integrate(U,[2,2,1], grid, corners, size_vec, pol )+SPIDER_integrate( V,[2,2,2], grid, corners, size_vec, pol ));
%   scales(a) = 1;
%  a         = a+1;
%
%   labels{a} = "u_j * di dj dk u_k, u_j * di dk^2 u_j";
%  G(:,a)    = SPIDER_integrate( dwLdy, [], grid, corners, size_vec, pol );
%  scales(a) = 1;
%  a         = a+1;
%
%  labels{a} = "dt_large scale vorticity";
%  G(:,a)    = SPIDER_integrate( dtwL, [], grid, corners, size_vec, pol);
%  scales(a) = 1;
%   a         = a+1;
%
% labels{a} = "dx2_large scale vorticity";
% G(:,a)    = SPIDER_integrate( dwLdx2, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "dy2_large scale vorticity";
% G(:,a)    = SPIDER_integrate( dwLdy2, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "dt2_wL large scale vorticity";
% G(:,a)    = SPIDER_integrate(dt2Wl, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "large scale u";
% G(:,a)    = SPIDER_integrate( uLx, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
%
% labels{a} = "dx large scale u";
% G(:,a)    = SPIDER_integrate( duLdx, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "dy large scale u";
% G(:,a)    = SPIDER_integrate( duLdy, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "d2x large scale u";
% G(:,a)    = SPIDER_integrate( duLdx2, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
%  labels{a} = "dy2 large scale u";
%  G(:,a)    = SPIDER_integrate(duLdy2, [], grid, corners, size_vec, pol );
%  scales(a) = 1;
%  a         = a+1;
%
%  labels{a} = "dt large scale u";
%  G(:,a)    = SPIDER_integrate(dtuL, [], grid, corners, size_vec, pol );
%  scales(a) = 1;
%  a         = a+1;
%
%  labels{a} = "Wl times wL";
%  G(:,a)    = SPIDER_integrate(Wl.*wL, [], grid, corners, size_vec, pol );
%  scales(a) = 1;
%  a         = a+1;
%
% labels{a} = "Wl times wS";
%  G(:,a)    = SPIDER_integrate(Wl.*wS, [], grid, corners, size_vec, pol );
%  scales(a) = 1;
%  a         = a+1;
%
%
%
% labels{a} = "Wl times uL";
% G(:,a)    = SPIDER_integrate( Wl.*uL, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "Wl times uS";
% G(:,a)    = SPIDER_integrate( Wl.*uS, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "d/dx(Wl*wL)";
% G(:,a)    = SPIDER_integrate(dx_Wl_wL, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "d/dx(Wl*wS)"; %derivative of a scalar product
% G(:,a)    = SPIDER_integrate(dx_Wl_wS, [], grid, corners, size_vec, pol );
% scales(a) = 1;
%
%
% labels{a} = "dt2_Wl";
% G(:,a)    = SPIDER_integrate(dt2Wl , [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "Dt1_Wl";
% G(:,a)    = SPIDER_integrate( dtWl, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "u\dwSdx"; %not possible
% G(:,a)    = SPIDER_integrate( U.*dwSdx, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "dx(Wl*wL)";
% G(:,a)    = SPIDER_integrate( dx_Wl_times_wL, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
% labels{a} = "dy(W*wL)";
% G(:,a)    = SPIDER_integrate( dy_Wl_times_wL, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "dx2(W*wL)";
% G(:,a)    = SPIDER_integrate(dx2_Wl_times_wL, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "dy2(W*wL)";
% G(:,a)    = SPIDER_integrate(dy2_Wl_times_wL, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%
% labels{a} = "v\dwSdx";
% G(:,a)    = SPIDER_integrate( V.*dwSdx, [], grid, corners, size_vec, pol );
% scales(a) = 1;
% a         = a+1;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DON't TOUCH PAST HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize
norm_vec = SPIDER_integrate( 0*U(:,:,:,1) + 1, [], grid, corners, size_vec, pol0 );      
norm_vec = repmat( norm_vec, [dof,1] );
G = G./norm_vec;
G = G./scales;

function vals = SPIDER_integrate_curl( U, V, derivs, grid, corners, size_vec, pol )
  vals = SPIDER_integrate( V, [derivs, 2], grid, corners, size_vec, pol ) ...
       - SPIDER_integrate( U, [derivs, 1], grid, corners, size_vec, pol );
end
