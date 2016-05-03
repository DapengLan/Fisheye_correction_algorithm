clc;clear;
%I0=extraction_fish('d2.jpg');
tic;
I0=imread('23.jpg');   %提取轮廓
R=(size(rgb2gray(I0),1))/2;
oriwidth=2*R;
oriheigh=2*R;
Newheigh=2*R;
Newwidth=2*R;
I1=uint8(zeros(Newheigh,Newwidth,3));
xc=R;
yc=R;

for  j=1:1:Newwidth  
    for i=1:1:Newheigh
       
       yi=j-xc;             %平移坐标系
       xh=i-yc;
       xk=(xh*sqrt(R^2-yi^2))/R ; %坐标变换
        
        xft=round(xk+xc);                        
        yft=round(yi+yc);
        if (xft<=0)||(xft>oriwidth)||(yft<=0)||(yft>oriheigh)
            I1(i,j,:)=0;
        elseif (xft > 1) && (xft <oriwidth) && (yft> 1) && (yft <oriheigh)
             xx= xk+xc;
             yy=yi+yc;
             a = double(uint16(xx)); 
             b = double(uint16(yy)); 
              x11 = double(I0(a,b,:));                     % x11 <- origImg(a,b)
            x12 = double(I0(a,b+1,:));                   % x12 <- origImg(a,b+1)
            x21 = double(I0(a+1,b,:));                   % x21 <- origImg(a+1,b)
            x22 = double(I0(a+1,b+1,:));                 % x22 <- origImg(a+1,b+1)
            I1(i,j,:) = uint8( (b+1-yy) * ((xx-a)*x21 + (a+1-xx)*x11) + (yy-b) * ((xx-a)*x22 +(a+1-xx) * x12) ); %双线性插值
            else                                                        
            I1(i,j,:) = I0(xft,yft,:);     % else,apply nearest interpolate
        end
    end
end
subplot(1,2,1),imshow(I0);
title('校正前');
subplot(1,2,2),imshow(I1);
title('校正后');
toc;
time=toc;