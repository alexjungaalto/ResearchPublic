close all; 
clear all; 

%N= 10 ; % dataset contains 10 images 
X = []; 

for iter=1:10
    
filename = sprintf('apple_%d.jpg', iter);
A = double(imread(filename));
  %  imshow(A); 
I = rgb2gray(A/255);
[rows, cols]= size(I); 
X = [X reshape(I,rows*cols,1)]; 
figure(iter);
imagesc(I)
colormap('gray');

end

W = pca (X'); 

d = 2 ; 
W = W(:,1:d); 

hat_X = (W*W')*X ; 
figure(11); 

imagesc(reshape(hat_X(:,5),rows,cols)); 
title('reconstruction of first image')
colormap('gray');