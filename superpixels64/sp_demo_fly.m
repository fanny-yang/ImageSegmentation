% sp_demo.m
%
% See instructions in README.
clear all;
cd '~/Google Drive/registration/superpixel/superpixels64/';

if ~exist('cncut')
    addpath('yu_imncut');
end

addpath('~/Google Drive/registration/superpixel/segbench/Detectors/');
addpath('~/Google Drive/registration/superpixel/segbench/Gradients/');
addpath('~/Google Drive/registration/superpixel/segbench/Util/');
addpath('~/Google Drive/registration/superpixel/segbench/Textons/');
addpath('~/Google Drive/registration/superpixel/segbench/Filters/');
addpath('~/Google Drive/registration/superpixel/superpixels64/yu_imncut/');
%addpath('~/Documents/fruitfly/siqi_201206');
addpath('~/Documents/fruitfly/image2TI');
%addpath('~/Documents/fruitfly/Visualization');
%addpath(genpath('~/Documents/fruitfly/Utility'));
%addpath('~/Documents/fruitfly/morphdet/');

inimg = imread('insitu76293.jpg');
inimg = imresize(inimg, .5);
[x y score] = fembryo(inimg, 3, 0);

figure;
imshow(inimg); hold on; 
plot(x,y,'-');
hold off; 

bwmat = zeros(size(inimg,1),size(inimg,2));
for l = 1:length(x)
    bwmat(y(l),x(l)) = 1;
end
bwmat_int = bwconvhull(bwmat);
phi = double((bwmat_int ==0).*bwdist(bwmat_int) - (bwmat_int ~=0).*bwdist(1 - bwmat_int));
phi1 = bwdist(bwmat);int_crude = phi1 <abs(phi);
%%


Ip = im2double(histeq(rgb2gray(inimg)));
I = zeros(size(inimg));
I(:,:,1) = Ip;I(:,:,3) = Ip;I(:,:,3) = Ip;

N = size(I,1);
M = size(I,2);

% Number of superpixels coarse/fine.
N_sp=40;
%N_sp2=1000;
% Number of eigenvectors.
N_ev=40;

% ncut parameters for superpixel computation
diag_length = sqrt(N*N + M*M);
par = imncut_sp;
par.int=0;
par.pb_ic=1;
par.sig_pb_ic=0.05;
par.sig_p=ceil(diag_length/50);
par.verbose=0;
par.nb_r=ceil(diag_length/60);
par.rep = -0.005;  % stability?  or proximity?
par.sample_rate=0.2;
par.nv = N_ev;
par.sp = N_sp;
par.mask1 = int_crude;

% Intervening contour using mfm-pb
fprintf('running PB\n');
[emag,ephase] = pbWrapper(I,par.pb_timing);
emag = pbThicken(emag);
par.pb_emag = emag;
par.pb_ephase = ephase;
clear emag ephase;

st=clock;
fprintf('Ncutting...');
[Seg, V,S,W] = imncut_sp_mod1(I,par);
fprintf(' took %.2f minutes\n',etime(clock,st)/60);


% st=clock;
% fprintf('Fine scale superpixel computation...');
% Sp2 = clusterLocations(Sp,ceil(N*M/N_sp2));
% fprintf(' took %.2f minutes\n',etime(clock,st)/60);


I = im2double(inimg);
%I_sp = segImage(I,Sp);
%I_sp2 = segImage(I,Sp2);
I_seg = segImage(I,Seg);
figure;
subplot(1,2,1);
imshow(I);
subplot(1,2,2);
imshow(I_seg);
% subplot(1,4,3);
% imshow(I_sp);
% subplot(1,4,4);
% imshow(I_sp2);


