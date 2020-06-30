function [ final_SR_Sum ] = Phase_corr_SR1( F_odd, F_even, params, Pw )
% phase correction, even and odd line 
% After SR data

%% Phase map correction

phase_diff1 = zeros(size(F_even));


for iiz = 1:params.Nz
phase_diff1(iiz,:,:,:,:) = F_even(iiz,:,:,:,:)./F_odd(iiz,:,:,:,:);
end

phase_diff = angle(phase_diff1);

F_even_corr = zeros(size(F_even));

for iiz = 1:params.Nz
F_even_corr(iiz,:,:,:,:) = exp(-i*phase_diff(iiz,:,:,:,:)).*F_even(iiz,:,:,:,:);
end

final_SR = F_odd+F_even_corr;
final_SR_Sum = permute(final_SR,[2,3,1,4,5]);

% Im_odd1 = zeros(size(F_odd));
% Im_even1 = zeros(size(F_odd));
% 
% %%%%%%%%%% re-inverse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for iir = 1:params.Nr
%     for iic = 1:params.Nc
%         for iiy = 1:params.Ny
% Im_odd1(:,:,iiy,iic,iir) = Pw*squeeze(F_odd(:,:,iiy,iic,iir));
% Im_even1(:,:,iiy,iic,iir) = Pw*squeeze(F_even_corr(:,:,iiy,iic,iir));
%         end
%     end
% end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ImT = Im_odd1+Im_even1;
% 
% final_reSR = zeros(size(ImT));
% for iir = 1:params.Nr
%     for iic = 1:params.Nc
%         for iiy = 1:params.Ny
% final_reSR(:,:,iiy,iic,iir) = Pw'*squeeze(ImT(:,:,iiy,iic,iir));
%         end
%     end
% end
% 
% final_SR_Sum = permute(final_reSR,[2,3,1,4,5]);

end

