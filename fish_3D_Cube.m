clc;
clear;
I0=imread('f.jpg');
I0=uint8(rgb2gray(I0));  
[heigh,width]=size(I0);
r=floor((min(heigh,width))/2);   %原图像半径
oriwidth=2*r;                    %原图像宽长
oriheigh=2*r;
nwidth=4*r;                      %定义校正后图像大小
nheigh=4*r;
%ewidth1=round(0.707*r);
%ewidth2=round(1.293*r);
%eheigh1=round(0.707*r);
%eheigh2=round(1.293*r);
I1=uint8(zeros(nheigh,nwidth));
p1=r;
m1=r;
n1=r;
xc=r;
yc=r;
u0=r;
v0=r;
%-------------立方体顶面校正-------%

for i=1:2*r  
    for j=1:2*r
        m=i-xc;
        n=j-yc;
        u=(r*m)/(sqrt(m^2+n^2+p1^2));
        v=(r*n)/(sqrt(m^2+n^2+p1^2));
        xft=round(u+u0);                        
        yft=round(v+v0);
        if (xft<=0)||(xft>oriwidth)||(yft<=0)||(yft>oriheigh)
            I1(r+i,r+j,:)=0;
        elseif (xft >0) && (xft <oriwidth) && (yft>0) && (yft <oriheigh)
             xx= u+u0;
             yy=v+v0;
             a = double(uint16(xx)); 
             b = double(uint16(yy)); 
              x11 = double(I0(a,b));                     % x11 <- origImg(a,b)
            x12 = double(I0(a,b+1));                   % x12 <- origImg(a,b+1)
            x21 = double(I0(a+1,b));                   % x21 <- origImg(a+1,b)
            x22 = double(I0(a+1,b+1));                 % x22 <- origImg(a+1,b+1)
            I1(r+i,r+j) = uint8( (b+1-yy) * ((xx-a)*x21 + (a+1-xx)*x11) + (yy-b) * ((xx-a)*x22 +(a+1-xx) * x12) ); %双线性插值
            else                                                        
            I1(r+i,r+j) = I0(xft,yft);     % else,apply nearest interpolate
        end      
    end  
end

%%-----------立方体左侧面校正------------------%%
for i=1:2*r
    for j=-r+1:0 
        p=j+xc;
        m=i-yc;
        u=(r*m)/(sqrt(m^2+n1^2+p^2));
        v=(r*(-n1))/(sqrt(m^2+n1^2+p^2));
        xft=round(u+u0);                        
        yft=round(v+v0);
        if (xft<=0)||(xft>oriwidth)||(yft<=0)||(yft>oriheigh)
            I1(r+i,r+j)=0;
        elseif (xft >0) && (xft <oriwidth) && (yft>0) && (yft <oriheigh)
             xx= u+u0;
             yy=v+v0;
             a = double(uint16(xx)); 
             b = double(uint16(yy)); 
              x11 = double(I0(a,b));                     % x11 <- origImg(a,b)
            x12 = double(I0(a,b+1));                   % x12 <- origImg(a,b+1)
            x21 = double(I0(a+1,b));                   % x21 <- origImg(a+1,b)
            x22 = double(I0(a+1,b+1));                 % x22 <- origImg(a+1,b+1)
            I1(r+i,j+r) = uint8( (b+1-yy) * ((xx-a)*x21 + (a+1-xx)*x11) + (yy-b) * ((xx-a)*x22 +(a+1-xx) * x12) ); %双线性插值
            else                                                        
            I1(r+i,j+r) = I0(xft,yft);     % else,apply nearest interpolate
        end
      
    end
   
end

%%-----------立方体右侧面校正------------------%%
for i=1:2*r  
    for j=1:r
        p=j-xc;
        m=i-yc;
        u=(r*m)/(sqrt(m^2+n1^2+p^2));
        v=(r*(n1))/(sqrt(m^2+n1^2+p^2));
        xft=round(u+u0);                        
        yft=round(v+v0);
        if (xft<=0)||(xft>oriwidth)||(yft<=0)||(yft>oriheigh)
            I1(r+i,3*r+j)=0;
        elseif (xft >0) && (xft <oriwidth) && (yft>0) && (yft <oriheigh)
             xx= u+u0;
             yy=v+v0;
             a = double(uint16(xx)); 
             b = double(uint16(yy)); 
              x11 = double(I0(a,b));                     % x11 <- origImg(a,b)
            x12 = double(I0(a,b+1));                   % x12 <- origImg(a,b+1)
            x21 = double(I0(a+1,b));                   % x21 <- origImg(a+1,b)
            x22 = double(I0(a+1,b+1));                 % x22 <- origImg(a+1,b+1)
            I1(r+i,3*r+j) = uint8( (b+1-yy) * ((xx-a)*x21 + (a+1-xx)*x11) + (yy-b) * ((xx-a)*x22 +(a+1-xx) * x12) ); %双线性插值
            else                                                        
            I1(r+i,3*r+j) = I0(xft,yft);     % else,apply nearest interpolate
        end
      
    end
   
end

%%-----------立方体上侧面校正------------------%%
for i=-r+1:0  
    for j=1:2*r
        p=i+xc;
        n=j-yc;
        u=(r*(-m1))/(sqrt(m1^2+n^2+p^2));
        v=(r*(n))/(sqrt(m1^2+n^2+p^2));
        xft=round(u+u0);                        
        yft=round(v+v0);
        if (xft<=0)||(xft>oriwidth)||(yft<=0)||(yft>oriheigh)
            I1(r+i,r+j)=0;
        elseif (xft >0) && (xft <oriwidth) && (yft>0) && (yft <oriheigh)
             xx= u+u0;
             yy=v+v0;
             a = double(uint16(xx)); 
             b = double(uint16(yy)); 
              x11 = double(I0(a,b));                     % x11 <- origImg(a,b)
            x12 = double(I0(a,b+1));                   % x12 <- origImg(a,b+1)
            x21 = double(I0(a+1,b));                   % x21 <- origImg(a+1,b)
            x22 = double(I0(a+1,b+1));                 % x22 <- origImg(a+1,b+1)
            I1(r+i,r+j) = uint8( (b+1-yy) * ((xx-a)*x21 + (a+1-xx)*x11) + (yy-b) * ((xx-a)*x22 +(a+1-xx) * x12) ); %双线性插值
            else                                                        
            I1(r+i,r+j) = I0(xft,yft);     % else,apply nearest interpolate
        end
      
    end
   
end

%%-----------立方体下侧面校正------------------%%
for i=1:r  
    for j=1:2*r
        p=i-xc;
        n=j-yc;
        u=(r*(m1))/(sqrt(m1^2+n^2+p^2));
        v=(r*(n))/(sqrt(m1^2+n^2+p^2));
        xft=round(u+u0);                        
        yft=round(v+v0);
        if (xft<=0)||(xft>oriwidth)||(yft<=0)||(yft>oriheigh)
            I1(3*r+i,r+j)=0;
        elseif (xft >0) && (xft <oriwidth) && (yft>0) && (yft <oriheigh)
             xx= u+u0;
             yy=v+v0;
             a = double(uint16(xx)); 
             b = double(uint16(yy)); 
              x11 = double(I0(a,b));                     % x11 <- origImg(a,b)
            x12 = double(I0(a,b+1));                   % x12 <- origImg(a,b+1)
            x21 = double(I0(a+1,b));                   % x21 <- origImg(a+1,b)
            x22 = double(I0(a+1,b+1));                 % x22 <- origImg(a+1,b+1)
            I1(3*r+i,r+j) = uint8( (b+1-yy) * ((xx-a)*x21 + (a+1-xx)*x11) + (yy-b) * ((xx-a)*x22 +(a+1-xx) * x12) ); %双线性插值
            else                                                        
            I1(3*r+i,r+j) = I0(xft,yft);     % else,apply nearest interpolate
        end
      
    end
   
end



%figure,imshow(I1);
figure;
subplot(1,2,1),imshow(I0);
subplot(1,2,2),imshow(I1);
