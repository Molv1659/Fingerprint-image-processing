clear all;
I = imread('23_2.bmp');
[M,N] = size(I);
I = double(I);
%Normalization
Mean = sum(I(:))/(size(I,1)*size(I,2));
I2 = I.^2;
Mean2 = sum(I2(:))/(size(I,1)*size(I,2));
VAR = Mean2 - Mean*Mean;
Mean0 = 100;
VAR0 = 2000;
for i=1:M
    for j=1:N  
        if(I(i,j)>Mean)
            I(i,j) = Mean0 + sqrt(VAR0*(I(i,j)-Mean)^2/VAR);
        else
            I(i,j) = Mean0 - sqrt(VAR0*(I(i,j)-Mean)^2/VAR);
        end
    end
end
%处理成8的倍数边长，四边扩张12，
%Orientation Image
I = imresize(I,[M-mod(M,8),N-mod(N,8)]);
enhance_I = zeros(size(I));
expand_I = padarray(I,[12,12],'symmetric');
[M,N] = size(expand_I);
Gx = [-1 0 1; -2 0 2; -1  0  1 ];
Gy = [ 1 2 1;  0 0 0; -1 -2 -1 ];
diffx = conv2(expand_I,Gx,'same');
diffy = conv2(expand_I,Gy,'same');
V = fast_local_mean((diffy+1i*diffx).^2,16);
%要取块中心
V = V(17:8:M-15,17:8:N-15);
Cita = 0.5*angle(V);
Cita0 = Cita; 
%Smoothen
temp = cos(2*Cita)+1i*sin(2*Cita);
temp = imfilter(temp,fspecial('gaussian',5));
Cita = 0.5*angle(temp);

figure(1);
subplot(421);imshow(I,[]);
subplot(423);quiver(cos(Cita0),sin(Cita0),0.7,'MaxHeadSize',0);
axis ij; axis image; axis([0 size(Cita,2) 0 size(Cita,1)]);  title('before smoothen Orientation Image');
subplot(424);quiver(cos(Cita),sin(Cita),0.7,'MaxHeadSize',0);
axis ij; axis image; axis([0 size(Cita,2) 0 size(Cita,1)]);  title('Orientation Image');

%frequency image
f = zeros(size(Cita,1),size(Cita,2));
p_maxf = zeros(size(Cita,1),size(Cita,2));
for i =1:size(Cita,1)
    for j=1:size(Cita,2)
        f_pic=abs(fftshift(fft2(expand_I(8*i-7:8*i+24,8*j-7:8*j+24))));
        f_pic(17,17)=0;
        [maxf,order] = max(f_pic(:));
        x = floor(order/32)+1;  
        y = mod(order,32);
        beta = angle((x-17)+1i*(y-17)) - Cita(i,j);
        while(beta>pi)
            beta = beta-pi;
        end
        if(beta<0)
            beta = -beta;
        end
        if(maxf>1000 && beta>pi/3 && beta<2*pi/3) %配合着方向图看f，只有主频率和方向图符合的才要
            f(i,j) = sqrt( (x-17)^2+(y-17)^2 );
        else
            f(i,j) = 0;
        end
        p_maxf(i,j) = maxf;
    end
end
f0 = f;
fun1 = @(x) median(x(:));
f = nlfilter(f,[4 4] ,fun1);
f = imfilter(f,fspecial('gaussian',5));

f = im2double(uint8(f));
subplot(425),imshow(f0,[]);title('before smoothen frequency image');
subplot(426),imshow(f,[]);title('frequency image');


%region mask
mask = zeros(size(Cita,1),size(Cita,2));
mask = imbinarize(mask);
for i=1:size(Cita,1)
    for j=1:size(Cita,2)
        if(f(i,j)>0.0079 && f(i,j)<0.02)
            mask(i,j) = 1;
        end
    end
end
mask = imresize(mask,size(I));
mask0 = mask;
%先腐蚀散着的白点，在膨胀保证指纹都进来
se = strel('disk',10);        
mask = imerode(mask,se);
se = strel('disk',30);
mask = imdilate(mask,se);

citad = Cita*180/pi;
for i = 1:size(Cita,1)
    for j =1:size(Cita,2)
        if(ne(f(i,j), 0))
            gama = Cita(i,j);
            gama = gama + pi/2;
            gama = -gama;
            gama = gama*180/pi;
            [mag,phase] = imgaborfilt(I(8*i-7:8*i,8*(j-1)+1:8*j),0.1/f(i,j),gama);
            enhance_I(8*i-7:8*i,8*(j-1)+1:8*j) = mag.*cos(phase);
        end        
    end
end
enhance_I = im2double(uint8(enhance_I));
enhance_I = imbinarize(enhance_I,0.2);
I = im2double(uint8(I));
enhance_I = enhance_I.*mask + I.*(1-mask);
for time=1:7
    se = strel('disk',1);        
    enhance_I = imerode(enhance_I,se);
    se = strel('disk',1);
    enhance_I = imdilate(enhance_I,se);
end


subplot(427),imshow(mask0);title('before mask');
subplot(428),imshow(mask);title('mask');
subplot(422),imshow(enhance_I,[]);









