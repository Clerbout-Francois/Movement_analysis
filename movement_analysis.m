clear all
close all
clc

%% For the lecture of sequence of seq1
WORKDIR = pwd
SEQPREF = 'seq1'; % name of the directory
SEQDIR = strcat(WORKDIR,'\',SEQPREF);


nb_im = 50; % number of images in the directory 


%% Step needed in order to select all the images from the seq1 directory (from 01 to 50)
for i=1:nb_im
    n = mod(i,10); % remainder after division
    m = floor(i/10); % round toward negative infinity
    numero(i).str=[int2str(m),int2str(n)]; % as images are named seq1_xx, we want to have two numbers for images 1 to 9, that's why we use mod() and floor()
end

for i=1:nb_im
    ima(i).name = strcat(SEQDIR,'\',SEQPREF,'_',numero(i).str,'.jpg'); % we extract each image
    im1 = imread(ima(i).name); 
    ima(i).im=double(im1(:,:,1));
end

%% Start of the analysis

i0 = 13; % number of the image on which we will focalize our study

% Initialization : Estimation of noise standard deviation
diff0 = ima(i0).im - ima(i0-1).im; % difference image
diff0_carre = diff0 .* diff0;
sig_noise = mean(diff0_carre(:)) ^ 0.5

% The noise standard deviation in the difference image is estimated by robust estimation
s=mean(diff0_carre(:))*16;
nb_iter=10;
for i = 1:nb_iter
   sig_noise = sum((diff0_carre(:) <= s) .* diff0_carre(:)) / sum(diff0_carre(:) <= s);
   s = 16*sig_noise;
   sqrt(sig_noise)
end;
sig_noise = sqrt(sig_noise)




%% 1st part: segmentation based on difference information

% we want to determine the difference between ima(i0) and the previous image
% Threshold : k_th    
k_th = 4;
delta_i = 5; % decrement or increment on the frame number for the difference

diff1 = ima(i0).im - ima(i0-delta_i).im; % In our case, difference between i0=13 and image 8
diff1_carre = diff1 .* diff1;
diff1th = abs(diff1) > (k_th * sig_noise);
figure(11); colormap(gray(256));image(ima(i0).im);
title(sprintf('Image n° %d',i0));
saveas(11,sprintf('resultats/11_image_num_%d.jpg',i0));
figure(12); colormap(gray(256));image(ima(i0-delta_i).im);
title(sprintf('Image n° %d - %d',i0,delta_i));
saveas(12,sprintf('resultats/12_image_num_%d-%d.jpg',i0,delta_i));
figure(13); colormap(gray(256));image(abs(diff1+128));
title(sprintf('Image difference, delta-i = %d',delta_i));
saveas(13,sprintf('resultats/13_diff_image___delta_i_%d.jpg',delta_i));
figure(14); colormap(gray(256));imagesc(diff1th);
title(sprintf('Thresholded image difference, i0 =%d, delta-i = %d, k-th = %.1f',i0,delta_i,k_th));
saveas(14,sprintf('resultats/14_th_diff_image___delta_i_%d___k_th_%.1f.jpg',delta_i,k_th));



% Same operation between ima(i0) and the following image.
% then "logical and" between the two images
diff2 = ima(i0).im - ima(i0+delta_i).im;
diff2_carre = diff2 .* diff2;
diff2th = abs(diff2) > (k_th * sig_noise);
diffand = diff1th .* diff2th;

figure(21); colormap(gray(256));image([ima(i0-delta_i).im, ima(i0).im, ima(i0+delta_i).im]);
title(sprintf('Images n° %d-%d ,%d , %d+%d',i0,delta_i,i0,i0,delta_i));
saveas(21,sprintf('resultats/21_images_num_%d_%d_%d.jpg',i0-delta_i,i0,i0+delta_i));
figure(22); colormap(gray(256));image([ abs(diff1+128), abs(diff2+128)]);
title(sprintf('Images difference  ima%d-ima%d ,ima%d-ima%d ',i0,i0-delta_i,i0,i0+delta_i));
saveas(22,sprintf('resultats/22_diff_images___delta_i_%d.jpg',delta_i));
figure(23); colormap(gray(256));imagesc([diff1th, diff2th, diffand] );
title(sprintf('Thresholded images difference (k-th=%.1f) :  diff1=ima%d-ima%d ,diff2=ima%d-ima%d , C=diff1&diff2',k_th,i0,i0-delta_i,i0,i0+delta_i));
saveas(23,sprintf('resultats/23_th_diff_images___delta_i_%d.jpg',delta_i));
figure(24); colormap(gray(256));imagesc(diffand);
title(sprintf('Logical AND of thresholded images difference, i0 =%d, delta-i = %d, k-th = %.1f',i0,delta_i,k_th));
saveas(24,sprintf('resultats/24_th_diffAND_image___delta_i_%d___k_th_%.1f.jpg',delta_i,k_th));


%% 2nde part : determination of a background image
% 1- Calculation of the mean image and the standard deviation image.
for i=1:nb_im
    if (i == 1) imbgd = ima(i).im; imbgdvar = ima(i).im .* ima(i).im;
    else imbgd = imbgd + ima(i).im; imbgdvar = imbgdvar + ima(i).im .* ima(i).im;
    end
end
imbgd = imbgd / nb_im; imbgdvar = imbgdvar / nb_im;
imbgdsigma = (imbgdvar- imbgd .* imbgd ).^0.5;
figure(101);image(imbgd); colormap(gray(256));
title('Image of averages');
saveas(101,'resultats/101_Image_of_averages.jpg');
figure(102);imagesc(imbgdsigma); colormap(gray(256));
title('Image of standard deviations');
saveas(102,'resultats/102_Image_of_standard_deviations.jpg');

% Background re-estimation by recursive median
imbgd0 = imbgd;
for i=1:nb_im
    imbgd = imbgd + (ima(i).im > imbgd) - (ima(i).im < imbgd); 
end
figure(103);image(imbgd); colormap(gray(256));
saveas(103,'resultats/103_Image_medians.jpg');

% 2nd iteration of the recursive median
imbgd1 = imbgd;
for i=1:nb_im
    imbgd = imbgd + (ima(i).im > imbgd) - (ima(i).im < imbgd); 
end
figure(104);image(imbgd); colormap(gray(256));
title('Image of recursive medians');
saveas(104,'resultats/104_Image_medians.jpg');

% difference image - estimated background
figure(111);image(abs(ima(i0).im - imbgd0)); colormap(gray(256));
title('Image minus average');
saveas(111,'resultats/111_Image_minus_average.jpg');
figure(112);image(abs(ima(i0).im - imbgd)); colormap(gray(256));
title('Image minus median');
saveas(112,'resultats/112_Image_minus_median.jpg');
figure(113);image(abs(diff1)); colormap(gray(256));
title(sprintf('Image difference (Absolute value), delta-i = %d',delta_i));
saveas(113,sprintf('resultats/113_abs_diff_image___delta_i_%d.jpg',delta_i));

% thresholded image difference

figure(121);imagesc(abs(ima(i0).im - imbgd0)>(k_th * sig_noise)); colormap(gray(256));
title(sprintf('(Image - Average) thresholded, k-th = %.1f',k_th));
saveas(121,sprintf('resultats/121_th_diff_image_minus_average___k_th_%.1f.jpg',k_th));
figure(122);imagesc(abs(ima(i0).im - imbgd)>(k_th * sig_noise)); colormap(gray(256));
title(sprintf('(Image - Mediann) thresholded, k-th = %.1f',k_th));
saveas(122,sprintf('resultats/122_th_diff_image_minus_median___k_th_%.1f.jpg',k_th));
figure(123);imagesc(abs(diff1)>(k_th * sig_noise)); colormap(gray(256));

figure(124);imagesc(diffand); colormap(gray(256));

%% Alternate filtering
% Opening then closing
ima1 = abs(ima(i0).im - imbgd)>(k_th * sig_noise);
el_struct = ones(3,3);
ima1 = imerode(ima1,el_struct);
ima1 = imdilate(ima1,el_struct);
ima1 = imdilate(ima1,el_struct);
ima1 = imerode(ima1,el_struct);
figure(211); imagesc(ima1); colormap(gray(256));
title(sprintf('Opening then closing ((Image - Median) thresholded, k-th = %.1f',k_th));
saveas(211,sprintf('resultats/211_Ouv_Ferm_th_diff_image_minus_median___k_th_%.1f.jpg',k_th));

% Closing then opening
ima2 = abs(ima(i0).im - imbgd)>(k_th * sig_noise);
el_struct = ones(3,3);
ima2 = imdilate(ima2,el_struct);
ima2 = imerode(ima2,el_struct);
ima2 = imerode(ima2,el_struct);
ima2 = imdilate(ima2,el_struct);
figure(212); imagesc(ima2); colormap(gray(256));
title(sprintf('Closing then opening ((Image - Median) thresholded), k-th = %.1f',k_th));
saveas(212,sprintf('resultats/212_Ferm_Ouv_th_diff_image_minus_median___k_th_%.1f.jpg',k_th));

% Segmentation

ima24 = myCompCon(ima2,4,0);
figure(213); imagesc(ima24); colormap(gray(256));
title(sprintf('RCC4 after Closing then opening ((Image - Median) thresholded), k-th = %.1f',k_th));
saveas(213,sprintf('resultats/213_RCC4_Ferm_Ouv_th_diff_image_minus_median___k_th_%.1f.jpg',k_th));
ima28 = myCompCon(ima2,8,0);
figure(214); imagesc(ima28); colormap(gray(256));
title(sprintf('RCC4 after Closing then opening ((Image - Median) thresholded), k-th = %.1f',k_th));
saveas(214,sprintf('resultats/214_RCC8_Ferm_Ouv_th_diff_image_minus_median___k_th_%.1f.jpg',k_th));

ima3 = abs(ima(i0).im - imbgd)>(k_th * sig_noise);
ima34 = myCompCon(ima3,4,0);
figure(215); imagesc(ima34); colormap(gray(256));
title(sprintf('RCC4 on ((Image - Median) thresholded), k-th = %.1f',k_th));
saveas(215,sprintf('resultats/215_RCC4_th_diff_image_minus_median___k_th_%.1f.jpg',k_th));
ima38 = myCompCon(ima3,8,0);
figure(216); imagesc(ima38); colormap(gray(256));
title(sprintf('RCC8 on ((Image - Median) thresholded), k-th = %.1f',k_th));
saveas(216,sprintf('resultats/216_RCC8_th_diff_image_minus_median___k_th_%.1f.jpg',k_th));







