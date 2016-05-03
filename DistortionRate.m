function In=DistortionRate(G)

pi=3.141592653589793;
G2=rgb2gray(G);
[height,width]=size(G2);
In=zeros(height,width,3);
%===============
% r_m=sqrt((height^2+width^2)/4);
r_m =sqrt((height/2)^2 + (width/2)^2);

theta_m = 6.35*pi/18.0;%(66.3/90)*(pi/2);
h = r_m/tan(theta_m);
Xc=round(height/2);
Yc=round(width/2);
% Xc=488+480;
% Yc=670+640;
for i=1:height
    for j=1:width
        r=sqrt(((i-Xc)^2+(j-Yc)^2)/1.0);
        theta=atan(r/h);
        a=tan(theta);
        if abs(a)>0.00000001;
            D=(theta-a)/a;
            
            x1=ceil((1+D)*(i-Xc));
            y1=ceil((1+D)*(j-Yc));
            x2=floor((1+D)*(i-Xc));
            y2=floor((1+D)*(j-Yc));
            
            if x1+Xc<height && y1+Yc<width && x2+Xc<height && y2+Yc<width
                In(i,j,:)=(G(x1+Xc,y1+Yc,:) + G(x2+Xc,y1+Yc,:) + G(x1+Xc,y2+Yc,:) + G(x2+Xc,y2+Yc,:)) / 4;
            end            
        else
            x=i;
            y=j;   
            if x+Xc<height && y+Yc<width
                In(i,j,:)=G(x+Xc,y+Yc,:);
            end
        end
%         if x+Xc<height && y+Yc<width
%             In(i,j)=G(x+Xc,y+Yc);
%        end
    end
end



