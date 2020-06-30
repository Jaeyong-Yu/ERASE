clear all; close all;clc
% tic

% path(pathdef);
% addpath 'utils'
%% data path

fid_folder = 'D:\JaeyongYu\data\20200615_MB\meas_MID00276_FID35057_ERASE_RS2x2x1_ss30_R82_TP23_12ms_TE30_FOV2.dat'; % data location
params.rvalue = 82;

%%%%%%%%%%%%%%%%%%%%%%%%% dicom direc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savenumber = 1010;
save_dir_name = '20200423ERASE_R146_full';
direc = 'D:\JaeyongYu\data\20200423\postrecon\';
fn = fid_folder;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% flag %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = 0; % Save data
dicom_flag = 0; % dicom import
stimulate_flag = 0; % stimulate program data import

% reconstruction flag
params.QV = 1; % Quadratic phase map, 0: High resolution P matirx, 1: simple P matrix (default)
params.WV = 1; % Weighting function, 0: center quadratic phase, 1: gaussian (defalt)
params.alpha = 1.5; % a lager value of 'alpha' produces a more params.Narrow in WV.1, P matrix weighting factor
params.Posi_FM_dir = 1; % FM positive direction
params.Nega_FM_dir = 0; % FM negative direction
params.PI = 1; % acceleration factor of PI

% reconstructed image selection [1: SUM & CH sum, 2: reSR & CH sum]
data_select = 1; % for dicom & save .mat file & image plot & total map 

% total map
image_flag = 1; % Total image 
tSNR_flag = 0; % Total tSNR image 
indexI =1; % Total image index
indexT =1; % Total tSNR index
nx = 8; 
ny = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reading raw data from dat files

twix_obj = mapVBVD(fid_folder,'removeOS');
A = max(size(twix_obj));
if params.PI ~= 1
    if A == 1
        DATA_ms = twix_obj.image.dataSize;
        rawData = twix_obj.image();
        % rawData_params.Navi = twix_obj.phasecor();
        rawData_PI = twix_obj.refscan();
    else
        DATA_ms = twix_obj{1, A}.image.dataSize;
        rawData = twix_obj{1, A}.image();
        rawData_PI = twix_obj{1, A}.refscan();
    end
else
    if A == 1
        DATA_ms = twix_obj.image.dataSize;
        rawData = twix_obj.image();
        % rawData_params.Navi = twix_obj.phasecor();
    else
        DATA_ms = twix_obj{1, A}.image.dataSize;
        rawData = twix_obj{1, A}.image();
    end
end

params.Nx = DATA_ms(1); % RO
params.Nc = DATA_ms(2); % chanparams.Nal
params.Ny = DATA_ms(3); % PE
params.Nz = DATA_ms(4); % SPEN dir
params.Na = DATA_ms(6); % average
params.Nr = DATA_ms(9); % repetition
params.Nset = DATA_ms(10); % sets
params.Ns = DATA_ms(11); % segments

raw = squeeze(double(permute(rawData, [3 1 2 4 9 10 11 6 5 7 8])));
% PE(params.Ny) RO(params.Nx) C(params.Nc) SPEN(params.Nz) Reptition(params.Nr) ... 
% sets(params.params.Nset) segment(params.Ns) average(params.Na) slice phase contrast

if params.PI ~= 1
    raw_PI = squeeze(double(permute(rawData_PI, [3 1 2 4 9 10 11 6 5 7 8])));
    for p = 1:(params.PI-1)
        raw(params.Ny+p,:,:,:,:) = zeros(size(raw(1,:,:,:,:)));
    end
end
% center_line RO(params.Nx) C(params.Nc) SPEN(params.Nz) Reptition(params.Nr) ... 
% sets(params.params.Nset) segment(params.Ns) average(params.Na) slice phase contrast
%% In-plane GRAPPA reconstruction (PI)
if params.PI ~= 1
    params.Ny = params.Ny + params.PI - 1;
res = zeros(size(raw));
% 
% if params.Nr == 1
% 
% for z = 1 : params.Nz
%     res(:,:,:,z) = GRAPPA(raw(:,:,:,z),raw_PI(:,:,:,z),[5 5],0.1);
%     z
% end
% end
% 
% if params.Nr ~= 1
% 
% for t = 1 : params.Nr
% for z = 1 : params.Nz
%     res(:,:,:,z,t) = GRAPPA(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
%     t
%     z
% end
% end
% end

tic;
aaa = 1;
for t = 1 :  params.Nr
    for z = 1 : params.Nz
        para(1,aaa) = t;
        para(2,aaa) = z;
        aaa = aaa +1;
    end
end

% NumPool = 4; % Number of Parallel computing
% isOpen = parpool('local',NumPool) > 0;

parfor tt = 1 : 10
    t = para(1,tt);
    z = para(2,tt);
    res(:,:,:,tt) = GRAPPA(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
%     res(:,:,:,tt) = GRAPPA1(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
end

parfor tt = 11 : 20
    t = para(1,tt);
    z = para(2,tt);
    res(:,:,:,tt) = GRAPPA(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
%     res(:,:,:,tt) = GRAPPA1(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
end

parfor tt = 21 : length(para)
    t = para(1,tt);
    z = para(2,tt);
    res(:,:,:,tt) = GRAPPA(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
%     res(:,:,:,tt) = GRAPPA1(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
end

% parfor tt = 31 : length(para)
%     t = para(1,tt);
%     z = para(2,tt);
%     res(:,:,:,tt) = GRAPPA(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
% %     res(:,:,:,tt) = GRAPPA1(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
% end
toc;
% parfor tt = 1 : length(para)
%     t = para(1,tt);
%     z = para(2,tt);

%     res(:,:,:,tt) = GRAPPA(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
% %     res(:,:,:,tt) = GRAPPA1(raw(:,:,:,z,t),raw_PI(:,:,:,z,t),[5 5],0.1);
% end

Recon_ksp = permute(res,[2 4 1 3 5]);

else
    Recon_ksp = permute(raw,[2 4 1 3 5]);
end

% [ low_raw ] = LowResKernel(raw, 36);
% Recon_ksp = permute(low_raw,[2 4 1 3 5]);
% Recon_ksp = Recon_ksp(:,:,:,:,6);
% params.Nr=1;
%% SR Reconstruction

[final_SR_tot final_reSR_tot final_SR_Sum] = SR_reconstruction(Recon_ksp,params);

%% data selection 1: SUM, 2: reSR, 3: SUM & CH sum, 4: reSR & CH sum

    clear final_SR

if data_select == 1
final_SR = final_SR_tot;
end

if data_select == 2
final_SR = final_reSR_tot;
end

%% reconstruction image plot

%%%%%%%%%%%%%%%%%%%%%% plot,total SR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, imagesc(abs(squeeze(final_SR_tot(ceil(params.Nx/2),:,:,1,1)))), colormap(gray), title('SR images (SPEN & PE)'), set(gcf, 'color', [1,1,1]), axis('equal')
figure, imagesc(abs(squeeze(final_SR_tot(:,:,ceil(params.Nz/2),1,1)))), colormap(gray), title('SR images (RO & PE)'), set(gcf, 'color', [1,1,1]), axis('equal')
figure, imagesc(abs(squeeze(final_SR_tot(:,ceil(params.Ny/2),:,1,1)))), colormap(gray), title('SR images (RO & SPEN)'), set(gcf, 'color', [1,1,1]), axis('equal')
% figure, imagesc(abs(squeeze(final_SR_tot(:,:,3,1,1)))), colormap(gray), title('SR images (RO & PE) 3'), set(gcf, 'color', [1,1,1]), axis('equal')
% figure, imagesc(abs(squeeze(final_SR_tot(:,:,4,1,1)))), colormap(gray), title('SR images (RO & PE) 4'), set(gcf, 'color', [1,1,1]), axis('equal')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% plot,total SR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, imagesc(abs(squeeze(final_reSR_tot(ceil(params.Nx/2),:,:,1,1)))), colormap(gray), title('SR images (SPEN & PE)'), set(gcf, 'color', [1,1,1]), axis('equal')
% figure, imagesc(abs(squeeze(final_reSR_tot(:,:,ceil(params.Nz/2),1,1)))), colormap(gray), title('SR images (RO & PE)'), set(gcf, 'color', [1,1,1]), axis('equal')
% figure, imagesc(abs(squeeze(final_reSR_tot(:,ceil(params.Ny/2),:,1,1)))), colormap(gray), title('SR images (RO & SPEN)'), set(gcf, 'color', [1,1,1]), axis('equal')
% figure, imagesc(abs(squeeze(final_reSR_tot(:,:,3,1,1)))), colormap(gray), title('SR images (RO & PE) 3'), set(gcf, 'color', [1,1,1]), axis('equal')
% figure, imagesc(abs(squeeze(final_reSR_tot(:,:,4,1,1)))), colormap(gray), title('SR images (RO & PE) 4'), set(gcf, 'color', [1,1,1]), axis('equal')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prenormalized 
% addpath 'C:\Users\Work1\Documents\MATLAB\LIBRARY'
% ssize = size(final_SR);
% mask = zeros(ssize(1),ssize(2),ssize(3));
% mask(find(final_SR > 0.08*max(final_SR(:)))) = 1;
% figure, mosaic(mask,6,6,1)
% 
% th = 15;
% for sl=1:ssize(3)
%     fit1 = polyfit2D_NthOrder(final_SR(:,:,sl),double(mask(:,:,sl)),2);
%     fit_new = 1./fit1;
%     fit_final = fit_new;
%     fit_final(fit_new>th) = th;
%     disp_new(:,:,sl) = final_SR(:,:,sl).*fit_final;
% end


%% tSNR map

if tSNR_flag == 1

FF11 = abs(squeeze(final_SR(:,:,:,:)));
ssize = size(FF11,4);

final_mean = mean(FF11,4);
final_sdt = std(FF11,0,4);
final_tSNR = final_mean./final_sdt;

sx = size(final_tSNR,1);
sy = size(final_tSNR,2);

tSNRmap = zeros(sy*nx, sx*ny);
for iy = 1:nx
    for ix = 1:ny
        tSNRmap((iy-1)*sy+1:iy*sy, (ix-1)*sx+1:ix*sx) = final_tSNR(:,:,indexT).';
        indexT = indexT+1;
    end
end
    
figure, imagesc(abs(tSNRmap)); axis image; colormap 'jet'; caxis([0 15]);

end
%% image map

if image_flag == 1

final_SR_map = final_SR(:,:,:,1);
% final_SR_map1 = disp_new(:,:,:,1);

sx = size(final_SR_map,1);
sy = size(final_SR_map,2);

RASEmap = zeros(sy*nx, sx*ny);
% RASEmap_newdisplay = zeros(sy*nx, sx*ny);

for iy = 1:nx
    for ix = 1:ny
        RASEmap((iy-1)*sy+1:iy*sy, (ix-1)*sx+1:ix*sx) = final_SR_map(:,:,indexI).';
%         RASEmap_newdisplay((iy-1)*sy+1:iy*sy, (ix-1)*sx+1:ix*sx) = final_SR_map1(:,:,indexI).';
        indexI = indexI+1;
    end
end
    
figure, imagesc(abs(RASEmap)); axis image; colormap(gray), title('RASE');
% figure, imagesc(abs(RASEmap_newdisplay)); axis image; colormap(gray); caxis([0 4.5]);
% figure, imagesc(abs(RASEmap_newdisplay)); axis image; colormap(gray), title('RASE after Intensity correction'); 
end

%% Save data

if S==1
    
savename = ['final_SR' num2str(savenumber)];

eval(['final_SR' num2str(savenumber) '= final_SR;'])

save([savename '.mat'],['final_SR' num2str(savenumber)],'-v7.3');

end

%% Stimulate file

if stimulate_flag==1
    
    IM = final_SR;
    flag = 0;
    b2sdt_SPENjk( fn , IM, flag)

end
%% dicom import

if dicom_flag==1
    
    ssize = size(final_SR);
    import_data = abs(final_SR);
    
if params.Nr == 1    

import_data = import_data./max(max(max(import_data)));
reshape_import_data = reshape(import_data,[ssize(1),ssize(2), ssize(3)*1]);

    mkdir(direc,save_dir_name)
    clear j


for n=1:ssize(3)*1;
    
     j(:,:)=squeeze(reshape_import_data(:,:,n));

     dicomwrite(j,num2str(n));
 
%      wh=what;
%      curr_path=wh.path;
%      new_path=cat(2,curr_path,'\',save_dir_name);
new_path = [direc,save_dir_name];
     
     movefile(num2str(n),new_path)
       
end

else
    
import_data = import_data./max(max(max(max(import_data))));
reshape_import_data = reshape(import_data,[ssize(1),ssize(2), ssize(3)*ssize(4)]);

    mkdir(direc,save_dir_name)
    clear j

for n=1:ssize(3)*ssize(4);
    
     j(:,:)=squeeze(reshape_import_data(:,:,n));

     dicomwrite(j,num2str(n));

    new_path = [direc,save_dir_name];
     
     movefile(num2str(n),new_path)
       
end

end


end

toc