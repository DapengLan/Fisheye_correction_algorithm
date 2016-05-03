clc;clear all;
tic;
I0=imread('23.jpg');   %提取轮廓
R=round((size(rgb2gray(I0),1))/2);
[heigh,width]=size(rgb2gray(I0));
u0=R;
v0=R;
%theta=2*pi/36;
nwidth=2*R;
nheigh=2*R;
x0=round(nwidth/2);
y0=round(nheigh/2);
I1=uint8(zeros(nheigh,nwidth,3));
for i=1:nwidth
    for j=1:nheigh
        x=i-x0;
        y=j-y0;
        u=R*x/(sqrt(x^2+y^2+R^2));
        v=R*y/(sqrt(x^2+y^2+R^2));
       xft=round(u+u0);                        
        yft=round(v+v0);
        if (xft<=0)||(xft>width)||(yft<=0)||(yft>heigh)
            I1(i,j,:)=0;
        elseif (xft > 1) && (xft <width) && (yft> 1) && (yft <heigh)
             xx= u+u0;
             yy=v+v0;
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
figure;
subplot(1,2,1),imshow(I0);
subplot(1,2,2),imshow(I1);
toc;
time=toc;