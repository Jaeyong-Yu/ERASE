function [final_SR_tot,final_reSR_tot,final_SR_Sum] = SR_reconstruction(Recon_ksp, params)

rawData_reshape = reshape(Recon_ksp,[params.Nx, params.Nz*params.Ny*params.Nc*params.Nr*params.Nset*params.Ns*params.Na]);

Im_odd = zeros(size(rawData_reshape));
Im_even = zeros(size(rawData_reshape));

Im_odd(:,1:2:end) = rawData_reshape(:,1:2:end);
Im_even(:,2:2:end) = rawData_reshape(:,2:2:end);

Im_odd = reshape(Im_odd,[params.Nx,params.Nz,params.Ny,params.Nc,params.Nr*params.Nset*params.Ns*params.Na]);
Im_even = reshape(Im_even,[params.Nx,params.Nz,params.Ny,params.Nc,params.Nr*params.Nset*params.Ns*params.Na]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, imagesc(abs(squeeze(Im_odd(:,(params.Nz/2),:,1)))), title('k-space odd data')
% figure, imagesc(abs(squeeze(Im_even(:,(params.Nz/2),:,1)))), title('k-space even data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Im_odd_x = ifftshift(ifft(fftshift(Im_odd,1),[],1),1);
Im_odd_xy = ifftshift(ifft(fftshift(Im_odd_x,3),[],3),3);

Im_even_x = ifftshift(ifft(fftshift(Im_even,1),[],1),1);
Im_even_xy = ifftshift(ifft(fftshift(Im_even_x,3),[],3),3);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, imagesc(abs(squeeze(Im_odd_xy(:,(params.Nz/2),:,1)))), title('k-space odd data')
% figure, imagesc(abs(squeeze(Im_even_xy(:,(params.Nz/2),:,1)))), title('k-space even data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Im_raw = Im_odd_xy+Im_even_xy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, imagesc(abs(squeeze(Im_raw(:,:,(params.Ny/2),1)))), title('image odd+even data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Im_final = zeros(params.Nx, params.Nz, params.Ny, params.Nr);
% for n=1:params.Nc
%     Im_raw1=Im_raw(:,:,:,n, :);
%     Im_raw2=abs(Im_raw1).^2;
%     Im_final = Im_final+Im_raw2;
% end
% Im_final_total = sqrt(Im_final);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, imagesc(squeeze(Im_final_total(:,:,params.Ny/2))), colormap(gray), title('image odd+even'), axis image
% figure, imagesc(squeeze(Im_final_total(:,params.Nz/2,:))), colormap(gray), title('image odd+even'), axis image
% figure, imagesc(squeeze(Im_final_total(params.Nx/2,:,:))), colormap(gray), title('image odd+even'), axis image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SR, P matrix

%% set the spatial axes (SPEN direction)
kz = 1/params.Nz*(-1/2*(params.Nz-1):1/2*(params.Nz-1));

%% set matrix P
% modified from pFT paper, Eq[1]
pha = zeros(params.Nz,params.Nz);
for index= 1:params.Nz
    pha(index,:) = -kz.^2/2+kz(index)*kz;
end

if params.QV == 1
    
    P = exp(sqrt(-1)*2*pi*params.rvalue*pha);                                         % P matrix, pFT paper, Eq[14]

else
    
    Intfactor = linspace(1,params.Nz,params.rvalue);
    for index = 1:params.Nz
    FL = griddedInterpolant(1:params.Nz, pha(index,:));
    PHI(index, :) = FL(Intfactor);
    end
%     figure, imagesc(pha)
%     figure, imagesc(PHI)

P = exp(sqrt(-1)*2*pi*params.rvalue*PHI);                                         % P matrix, pFT paper, Eq[14]

end

%% Windowing, local k-space points

%%%%%%%%%%%%%%%%%%%%%%%% partial FT points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M1 = params.Nz^2/params.rvalue;                                    % Number of partial points, pFT paper. Eq[12]
M1 = round(M1);

if params.WV ~= 1
if mod(M1, 2) == 0
    M1 = M1-1;
else
    M1;
end
end

disp(['Partial points = ', num2str(M1)             ])

%%%%%%%%%%%%%%%%%%%%%%%%%%% windowing matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = zeros(params.Nz, params.Nz);

for iic=1:params.Nz
    for iir=1:params.Nz
        if (iir-((M1-1)/2)<=iic & iic<=iir+((M1-1)/2))
            w(iir, iic) = 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, imagesc(w), title('windowing matrix')           % pFT paper, Eq[20]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Weighting 1: Gaussian, others: center quadratic phase shape 

if params.WV == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Window(gauss) version %%%%%%%%%%%%%%%%%%%%%%
    w_length = M1;
    w_g = gausswin(w_length, params.alpha);                          % alpha default 2.5
    % 'alpha = 3'. refereparams.Nce by Chen paper
    padding_w = [zeros(size(w)); w; zeros(size(w));];
    w_points = find(diag(w));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure, imagesc(padding_w), title('zeros padding : weighting matrix')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ii=1:params.Nz
        w_col = padding_w(:,ii);
        center_point = params.Nz+w_points(ii);
        if mod(M1,2) == 1
            w_col(center_point-floor(M1/2):center_point+floor(M1/2)) = w_g;
        end
        if mod(M1,2) == 0
            w_col(center_point-floor(M1/2):center_point+floor(M1/2)-1) = w_g;
        end
        W1(:,ii) = w_col;
    end
    W = W1(params.Nz+1:params.Nz+params.Nz,:);

    W_maxreix = W;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure, imagesc(W), title('weighting matrix and windowing')
%     figure, plot(W(16,:))
%     figure, plot(W(:,16))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    W = P.*w; % not Weighting
    w_g1 = unwrap(angle(W(params.Nz/2,:)));
%     figure, plot(w_g1)
    w_g1 = w_g1./max(w_g1);
    w_g = w_g1(find(w_g1~=0));
%     figure, plot(w_g)
    
    padding_w = [zeros(size(w)); w; zeros(size(w));];
    w_points = find(diag(w));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure, imagesc(padding_w), title('zeros padding : weighting matrix')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii=1:params.Nz
        w_col = padding_w(:,ii);
        center_point = params.Nz+w_points(ii);
        if mod(M1,2) == 1
            w_col(center_point-floor(M1/2):center_point+floor(M1/2)) = w_g;
        end
        if mod(M1,2) == 0
            w_col(center_point-floor(M1/2):center_point+floor(M1/2)-1) = w_g;
        end
        W1(:,ii) = w_col;
    end
    W = W1(params.Nz+1:params.Nz+params.Nz,:);
    W_maxreix = W;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure, imagesc(W), title('weighting matrix and windowing')
%     figure, plot(W(16,:))
%     figure, plot(W(:,16))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



%% FM & gradient direction 

if params.Posi_FM_dir ==1
% windoowing and weighting matrix
W2 = fliplr(W);   % Spin-echo type : single 180 pulse, GRE type : no rephase type
% weighting matrix
W3 = fliplr(w);   % Spin-echo type : single 180 pulse, GRE type : no rephase type
end
if params.Nega_FM_dir ==1
% windoowing and weighting matrix
W2 = W;           % Spin-echo type : double 180 pulse, GRE type : rephase type
% weighting matrix
W3 = w;           % Spin-echo type : double 180 pulse, GRE type : rephase type
end 

Pw = W2.*P ; % weighting & window
% Pw = W3.*P ; % only weighting


% figure, imagesc(abs(Pw)), title('windowing & weighted matrix (Pw)')   

%% SR 
% lambda = 0.05

Im_odd_xy = double(permute(Im_odd_xy, [2 1 3 4 5])); % SPEN RO PE C R
Im_even_xy = double(permute(Im_even_xy, [2 1 3 4 5])); % SPEN RO PE C R

F_odd = zeros(size(Im_odd_xy));
F_even = zeros(size(Im_even_xy));

for iir = 1:params.Nr
    for iic = 1:params.Nc
        for iiy = 1:params.Ny
            F_odd(:,:,iiy,iic,iir) = Pw'*squeeze(Im_odd_xy(:,:,iiy,iic,iir));
            F_even(:,:,iiy,iic,iir) = Pw'*squeeze(Im_even_xy(:,:,iiy,iic,iir));
%             F_odd(:,:,iiy,iic,iir)  = (Pw+lambda*eye(size(Pw,1)))'*squeeze(Im_odd_xy(:,:,iiy,iic,iir));
%             F_even(:,:,iiy,iic,iir) = (Pw+lambda*eye(size(Pw,1)))'*squeeze(Im_even_xy(:,:,iiy,iic,iir));  
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%
% Im_odd_xy1 = reshape(Im_odd_xy, [size(Im_odd_xy,1), size(Im_odd_xy,2)*size(Im_odd_xy,3)*size(Im_odd_xy,4)*size(Im_odd_xy,5)]);
% Im_even_xy1 = reshape(Im_even_xy, [size(Im_even_xy,1), size(Im_even_xy,2)*size(Im_even_xy,3)*size(Im_even_xy,4)*size(Im_even_xy,5)]);
% 
% F_odd1 = zeros(size(Im_odd_xy1));
% F_even1 = zeros(size(Im_even_xy1));
% 
% F_odd1 = Pw'*squeeze(Im_odd_xy1);
% F_even1 = Pw'*squeeze(Im_odd_xy1);
% 
% F_odd = reshape(F_odd1, [size(Im_odd_xy,1), size(Im_odd_xy,2), size(Im_odd_xy,3), size(Im_odd_xy,4), size(Im_odd_xy,5)]);
% F_even = reshape(F_even1, [size(Im_even_xy,1), size(Im_even_xy,2), size(Im_even_xy,3), size(Im_even_xy,4), size(Im_even_xy,5)]);
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, imagesc(abs(F_odd(:,:,params.Ny/2,1,1))), colormap(gray), title('F(ODD) : RO & SPEN')
% figure, imagesc(abs(squeeze(F_odd(:,params.Nx/2,:,1,1)))), colormap(gray), title('F(ODD) : SPEN & PE')
% figure, imagesc(abs(squeeze(F_odd(params.Nz/2,:,:,1,1)))), colormap(gray), title('F(ODD) : RO & PE')

% figure, imagesc(abs(F_even(:,:,params.Ny/2,1,1))), colormap(gray), title('F(EVEN) : RO & SPEN')
% figure, imagesc(abs(squeeze(F_even(:,params.Nx/2,:,1,1)))), colormap(gray), title('F(EVEN) : SPEN & PE')
% figure, imagesc(abs(squeeze(F_even(params.Nz/2,:,:,1,1)))), colormap(gray), title('F(EVEN) : RO & PE')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Even Odd direction, Phase correction & re SR

% [ final_SR_Sum, final_reSR ] = Phase_corr_SR( F_odd, F_even, params, Pw );
[ final_SR_Sum ] = Phase_corr_SR1( F_odd, F_even, params, Pw );
% matirx = RO PE SPEN channel and volumes

%%%%%%%%%%%%%%%%%%%%%%%% plot,Sum SR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, imagesc(abs(squeeze(final_SR_Sum(ceil(x/2),:,:,1,1)))), colormap(gray), title('final Sum SR(SPEN & PE)'), set(gcf, 'color', [1,1,1]),axis off
% figure, imagesc(abs(squeeze(final_SR_Sum(:,:,ceil(y/2),1,1)))), colormap(gray), title('final Sum SR(RO & SPEN)'), set(gcf, 'color', [1,1,1]),axis off
% figure, imagesc(abs(squeeze(final_SR_Sum(:,ceil(z/2),:,1,1)))), colormap(gray), title('final Sum SR(RO & PE)'), set(gcf, 'color', [1,1,1]),axis off
% figure, imagesc(abs(squeeze(final_SR_Sum(:,:,10,1,1)))), colormap(gray), title('final Sum SR 3'), set(gcf, 'color', [1,1,1])
% figure, imagesc(abs(squeeze(final_SR_Sum(:,:,4,1,1)))), colormap(gray), title('final Sum SR 4'), set(gcf, 'color', [1,1,1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% plot,total SR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure, imagesc(abs(squeeze(final_reSR(ceil(x/2),:,:,1,1)))), colormap(gray), title('F total (SPEN & PE)'), set(gcf, 'color', [1,1,1])
% figure, imagesc(abs(squeeze(final_reSR(:,:,ceil(y/2),1,1)))), colormap(gray), title('F total (RO & SPEN)'), set(gcf, 'color', [1,1,1])
% figure, imagesc(abs(squeeze(final_reSR(:,ceil(z/2),:,1,1)))), colormap(gray), title('F total (RO & PE)'), set(gcf, 'color', [1,1,1])
% figure, imagesc(abs(squeeze(final_reSR(:,:,3,1,1)))), colormap(gray), title('F total 3'), set(gcf, 'color', [1,1,1])
% figure, imagesc(abs(squeeze(final_reSR(:,:,4,1,1)))), colormap(gray), title('F total 4'), set(gcf, 'color', [1,1,1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% channel sum (after SR)

if params.Nc ~= 1 & params.Nr == 1   
    
final_SR_ch = zeros(params.Nx, params.Ny, params.Nz);
% final_reSR_ch = zeros(params.Nx, params.Ny, params.Nz);

for n=1:params.Nc
    F_total1a=final_SR_Sum(:,:,:,n);
%     F_total1b=final_reSR(:,:,:,n);    
    F_total2a=abs(F_total1a).^2;
%     F_total2b=abs(F_total1b).^2;    
    final_SR_ch = final_SR_ch+F_total2a;
%     final_reSR_ch = final_reSR_ch+F_total2b;    
end

final_SR_tot = sqrt(final_SR_ch);
% final_reSR_tot = sqrt(final_reSR_ch);

end

if params.Nc ~= 1 & params.Nr ~= 1  
    
final_SR_ch = zeros(params.Nx, params.Ny, params.Nz, params.Nr);
% final_reSR_ch = zeros(params.Nx, params.Ny, params.Nz, params.Nr);

for n=1:params.Nc
    F_total1a=squeeze(final_SR_Sum(:,:,:,n,:));
%     F_total1b=squeeze(final_reSR(:,:,:,n,:));    
    F_total2a=abs(F_total1a).^2;
%     F_total2b=abs(F_total1b).^2;    
    final_SR_ch = final_SR_ch+F_total2a;
%     final_reSR_ch = final_reSR_ch+F_total2b;    
end

final_SR_tot = sqrt(final_SR_ch);
% final_reSR_tot = sqrt(final_reSR_ch);

end

final_SR_tot;
% final_reSR_tot;
final_reSR_tot=1;



end

