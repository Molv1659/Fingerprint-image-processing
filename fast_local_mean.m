

function I1 = fast_local_mean(img,r)
N = (2*r+1)^2;
I1 = fast_local_sum(img,r)/N;
end

function imgcul = img_integral(img)
imgcul=cumsum(img,1);
imgcul=cumsum(imgcul,2);
end

function I = fast_local_sum(img,r)
[m,n] = size(img);
I = zeros(m,n);
img = padarray(img,[r+1 r+1],'symmetric');
imgcul = img_integral(img);
for i=1:m
    for j=1:n
        i1 = i+r+1;
        j1 = j+r+1;
        I(i,j) = imgcul(i1+r,j1+r) - imgcul(i1+r,j1-r-1) - imgcul(i1-r-1,j1+r) + imgcul(i1-r-1,j1-r-1);        
    end
end

end


function I2 = fast_local_var(img,r)
N = (2*r+1)^2;
imgmean = fast_local_sum(img,r)/N;
squareimgmean = fast_local_sum(img.*img , r)/N;
I2 = squareimgmean - imgmean.*imgmean;
end