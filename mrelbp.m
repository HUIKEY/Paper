function mrelbp=mrelbp(image,radius,wr,neighbors,mapping)
H1=mrelbp_f(image,radius,wr,neighbors,mapping); 
H2=mrelbp_f(image,radius+1,wr,neighbors,mapping); 
H3=mrelbp_f(image,radius+3,wr+2,neighbors,mapping); 
H4=mrelbp_f(image,radius+7,wr+6,neighbors,mapping); 
mrelbp=[H1 H2 H3 H4];
end
%实现mrelbp_ci,mrelbp_ni,mrelbp_rd三个算子对局部纹理信息的提取
function mrelbp = mrelbp_f(image,radius,wr,neighbors,mapping)
a = 2*pi/neighbors;
d_image=double(image);
spoints=zeros(8,2);
spoints_1=zeros(8,2);
for i = 1:neighbors
    spoints(i,1) = -radius*sin((i-1)*a);
    spoints(i,2) = radius*cos((i-1)*a);
end
% Determine the dimensions of the input image.%确定输入图像的尺寸
[ysize, xsize] = size(image);
miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));
% Block size, each LBP code is computed within a block of size bsizey*bsizex每个LBP代码在大小为bsizey * bsizex的块内计算
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;%ceil：朝正无穷方向取整
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;
%在block里中心点的坐标
% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));
%检查block和img的大小
% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizex || ysize < bsizey)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end
% Calculate dx and dy;计算中心像素需要移动的距离
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.所有可以作为模板中心点的像素集合
C = image(origy:origy+dy,origx:origx+dx);%显示原图像中模板中心点的像素集



%mrelbp_ci
mrelbp_ci=zeros(dy+1,dx+1);
wc=3;                  %滤波半径=领域半径
B_1=medfilt2(C,[wc,wc]);
ave=mean2(B_1);             % 一幅图像的均值
[m,n]=size(B_1);
for i=1:m
    for j=1:n
        if(B_1(i,j)>=ave)
            mrelbp_ci(i,j)=1;
        else
            mrelbp_ci(i,j)=0;
        end
    end
end



%mrelbp_ni
ave_rp=zeros(dy+1,dx+1);          %领域的均值ave_rp                    
for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  %round四舍五入 ceil进一 round退一
  fy = round(y);cy = ceil(y); ry = round(y);  
  fx = round(x);cx = ceil(x); rx = round(x);
  %判断是否需要插值
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
      N = image(ry:ry+dy,rx:rx+dx);
  else
     ty = y - fy;
     tx = x - fx;
     w1 = roundn((1 - tx) * (1 - ty),-6);
     w2 = roundn(tx * (1 - ty),-6);
     w3 = roundn((1 - tx) * ty,-6) ;
     w4 = roundn(1 - w1 - w2 - w3, -6);
     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
     N = roundn(N,-4);
  end   
  B_2=medfilt2(N,[wr,wr]);
  B_2=double(B_2);
  ave_rp= ave_rp + B_2;
end
ave_rp=ave_rp/neighbors;


mrelbp_ni=zeros(dy+1,dx+1);
for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
 
  fy = round(y);cy = ceil(y); ry = round(y);  
  fx = round(x);cx = ceil(x); rx = round(x);
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
      N = image(ry:ry+dy,rx:rx+dx);
  else
     ty = y - fy;
     tx = x - fx;
     w1 = roundn((1 - tx) * (1 - ty),-6);
     w2 = roundn(tx * (1 - ty),-6);
     w3 = roundn((1 - tx) * ty,-6) ;
     w4 = roundn(1 - w1 - w2 - w3, -6);
     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
     N = roundn(N,-4);
  end
  B_2=medfilt2(N,[wr,wr]);
  D = B_2>= ave_rp; 
  % Update the result matrix.更新结果矩阵
  v = 2^(i-1);
  mrelbp_ni = mrelbp_ni + v*D;
end

%mrelbp_rd
r_1=wr;       %大圆半径滤波=领域半径
r_2=wr-1;     %小圆半径滤波=领域半径-1
mrelbp_rd=zeros(dy+1,dx+1);
radius_1=radius-1;    %小圆半径=领域半径-1  大圆半径=领域半径
for i = 1:neighbors
    spoints_1(i,1) = -radius_1*sin((i-1)*a);
    spoints_1(i,2) = radius_1*cos((i-1)*a);
end
for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  y_1 = spoints_1(i,1)+origy;
  x_1 = spoints_1(i,2)+origx;
  
  fy = floor(y); cy = ceil(y); ry = round(y); 
  fx = floor(x); cx = ceil(x); rx = round(x);
  fy_1 = floor(y_1); cy_1 = ceil(y_1); ry_1 = round(y_1);  
  fx_1 = floor(x_1); cx_1 = ceil(x_1); rx_1 = round(x_1);
  %判断大圈是否需要插值
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
      N = image(ry:ry+dy,rx:rx+dx);
  else
     ty = y - fy;
     tx = x - fx;
     w1 = roundn((1 - tx) * (1 - ty),-6);
     w2 = roundn(tx * (1 - ty),-6);
     w3 = roundn((1 - tx) * ty,-6) ;
     w4 = roundn(1 - w1 - w2 - w3, -6);
     N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
         w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
     N = roundn(N,-4);
  end
  
  %判断小圈是否需要插值
  if (abs(x_1 - rx_1) < 1e-6) && (abs(y_1 - ry_1) < 1e-6)
      N_1 = image(ry_1:ry_1+dy,rx_1:rx_1+dx);
  else
     ty_1 = y_1 - fy_1;
     tx_1 = x_1 - fx_1;
     w1 = roundn((1 - tx_1) * (1 - ty_1),-6);
     w2 = roundn(tx_1 * (1 - ty_1),-6);
     w3 = roundn((1 - tx_1) * ty_1,-6) ;
     w4 = roundn(1 - w1 - w2 - w3, -6);
     N_1 = w1*d_image(fy_1:fy_1+dy,fx_1:fx_1+dx) + w2*d_image(fy_1:fy_1+dy,cx_1:cx_1+dx) + ...
         w3*d_image(cy_1:cy_1+dy,fx_1:fx_1+dx) + w4*d_image(cy_1:cy_1+dy,cx_1:cx_1+dx);
     N_1 = roundn(N_1,-4);
  end
 
  %滤波器
  B_3_1=medfilt2(N,[r_1,r_1]);
  B_3_2=medfilt2(N_1,[r_2,r_2]);
  D = B_3_1 >= B_3_2; 
 
  % Update the result matrix.更新结果矩阵
  v = 2^(i-1);
  mrelbp_rd = mrelbp_rd + v*D;
end





%应用映射
bins = mapping.num;
for i = 1:size(mrelbp_rd,1)%r=size(A,1)该语句返回的时矩阵A的行数， c=size(A,2) 该语句返回的时矩阵A的列数。
    for j = 1:size(mrelbp_rd,2)
        mrelbp_ni(i,j) = mapping.table(mrelbp_ni(i,j)+1);
        mrelbp_rd(i,j) = mapping.table(mrelbp_rd(i,j)+1);
    end
end
%绘制直方图
mrelbp_ci=hist(mrelbp_ci(:),2);
mrelbp_ci=mrelbp_ci/sum(mrelbp_ci);
mrelbp_ni=hist(mrelbp_ni(:),0:(bins-1));
mrelbp_ni=mrelbp_ni/sum(mrelbp_ni);
mrelbp_rd=hist(mrelbp_rd(:),0:(bins-1));
mrelbp_rd=mrelbp_rd/sum(mrelbp_rd);
mrelbp=[mrelbp_ci mrelbp_ni mrelbp_rd];
end






 








            




