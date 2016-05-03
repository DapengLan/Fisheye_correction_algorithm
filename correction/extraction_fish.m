function New_image=extraction_fish(pic)

I0=imread(pic);
I0=uint8(I0);
r=I0(:,:,1);
g=I0(:,:,2);
b=I0(:,:,3);
I=0.59.*r+0.11.*g+0.3.*b; %像素亮度计算公式
I=uint8(I);
[Height,Width]=size(I);
Thre=30;  %预设阈值

for Row1=1:(Height/2);  %循环寻找圆形区域上边界  

    CurRow_Bright=I(Row1,:); 
    Max=max(CurRow_Bright); %求取最大亮度值
    Min=min(CurRow_Bright); %求取最小亮度值
    Lim=Max-Min;  %该扫描线的极限亮度差
     if (Lim>Thre)
         Ytop=Row1;
         break;
         
     end
     
end
     
for Row2=Height:-1:(Height/2);  %循环寻找圆形区域下边界  
    
    CurRow_Bright=I(Row2,:);
    Max=max(CurRow_Bright); %求取最大亮度值
    Min=min(CurRow_Bright); %求取最小亮度值
    Lim=Max-Min;  %该扫描线的极限亮度差
     if (Lim>Thre)
         Ybot=Row2;
         break;
     end
     
end

for Col1=1:(Width/2);  %循环寻找圆形区域左边界  
   
    CurCol_Bright = I(:,Col1);
    Max=max(CurCol_Bright); %求取最大亮度值
    Min=min(CurCol_Bright); %求取最小亮度值
    Lim=Max-Min;  %该扫描线的极限亮度差
     if (Lim>Thre),
         Xleft=Col1;
         break;
     end
     
end

for Col2=Width:-1:Width/2;  %循环寻找圆形区域右边界  
   
    CurCol_Bright = I(:,Col2);
    Max=max(CurCol_Bright); %求取最大亮度值
    Min=min(CurCol_Bright); %求取最小亮度值
    Lim=Max-Min;  %该扫描线的极限亮度差
     if (Lim>Thre),
         Xrig=Col2;
         break;
     end
         
end

%x0=fix(Xleft+ Xrig)/2;
%y0=fix(Ytop+Ybot)/2;
Rx=floor((Xrig-Xleft)/2);
Ry=floor((Ytop-Ybot)/2);
R=max(Rx,Ry);
New_image=uint8(zeros(2*R,2*R,3));
for i= 1:2*R
    for j=1:2*R
        New_image(i,j,:)=I0(Ytop+i,Xleft+j,:);
       
    end
end
figure,imshow(New_image);
%figure,imshow(I0);
%figure,imshow(New_image);

