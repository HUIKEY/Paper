function mrelbp=mrelbp(image,radius,wr,neighbors,mapping)
H1=mrelbp_f(image,radius,wr,neighbors,mapping); 
H2=mrelbp_f(image,radius+1,wr,neighbors,mapping); 
H3=mrelbp_f(image,radius+3,wr+2,neighbors,mapping); 
H4=mrelbp_f(image,radius+7,wr+6,neighbors,mapping); 
mrelbp=[H1 H2 H3 H4];
end
%ʵ��mrelbp_ci,mrelbp_ni,mrelbp_rd�������ӶԾֲ�������Ϣ����ȡ
function mrelbp = mrelbp_f(image,radius,wr,neighbors,mapping)
a = 2*pi/neighbors;
d_image=double(image);
spoints=zeros(8,2);
spoints_1=zeros(8,2);
for i = 1:neighbors
    spoints(i,1) = -radius*sin((i-1)*a);
    spoints(i,2) = radius*cos((i-1)*a);
end
% Determine the dimensions of the input image.%ȷ������ͼ��ĳߴ�
[ysize, xsize] = size(image);
miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));
% Block size, each LBP code is computed within a block of size bsizey*bsizexÿ��LBP�����ڴ�СΪbsizey * bsizex�Ŀ��ڼ���
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;%ceil�����������ȡ��
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;
%��block�����ĵ������
% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));
%���block��img�Ĵ�С
% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizex || ysize < bsizey)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end
% Calculate dx and dy;��������������Ҫ�ƶ��ľ���
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.���п�����Ϊģ�����ĵ�����ؼ���
C = image(origy:origy+dy,origx:origx+dx);%��ʾԭͼ����ģ�����ĵ�����ؼ�



%mrelbp_ci
mrelbp_ci=zeros(dy+1,dx+1);
wc=3;                  %�˲��뾶=����뾶
B_1=medfilt2(C,[wc,wc]);
ave=mean2(B_1);             % һ��ͼ��ľ�ֵ
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
ave_rp=zeros(dy+1,dx+1);          %����ľ�ֵave_rp                    
for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  %round�������� ceil��һ round��һ
  fy = round(y);cy = ceil(y); ry = round(y);  
  fx = round(x);cx = ceil(x); rx = round(x);
  %�ж��Ƿ���Ҫ��ֵ
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
  % Update the result matrix.���½������
  v = 2^(i-1);
  mrelbp_ni = mrelbp_ni + v*D;
end

%mrelbp_rd
r_1=wr;       %��Բ�뾶�˲�=����뾶
r_2=wr-1;     %СԲ�뾶�˲�=����뾶-1
mrelbp_rd=zeros(dy+1,dx+1);
radius_1=radius-1;    %СԲ�뾶=����뾶-1  ��Բ�뾶=����뾶
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
  %�жϴ�Ȧ�Ƿ���Ҫ��ֵ
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
  
  %�ж�СȦ�Ƿ���Ҫ��ֵ
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
 
  %�˲���
  B_3_1=medfilt2(N,[r_1,r_1]);
  B_3_2=medfilt2(N_1,[r_2,r_2]);
  D = B_3_1 >= B_3_2; 
 
  % Update the result matrix.���½������
  v = 2^(i-1);
  mrelbp_rd = mrelbp_rd + v*D;
end





%Ӧ��ӳ��
bins = mapping.num;
for i = 1:size(mrelbp_rd,1)%r=size(A,1)����䷵�ص�ʱ����A�������� c=size(A,2) ����䷵�ص�ʱ����A��������
    for j = 1:size(mrelbp_rd,2)
        mrelbp_ni(i,j) = mapping.table(mrelbp_ni(i,j)+1);
        mrelbp_rd(i,j) = mapping.table(mrelbp_rd(i,j)+1);
    end
end
%����ֱ��ͼ
mrelbp_ci=hist(mrelbp_ci(:),2);
mrelbp_ci=mrelbp_ci/sum(mrelbp_ci);
mrelbp_ni=hist(mrelbp_ni(:),0:(bins-1));
mrelbp_ni=mrelbp_ni/sum(mrelbp_ni);
mrelbp_rd=hist(mrelbp_rd(:),0:(bins-1));
mrelbp_rd=mrelbp_rd/sum(mrelbp_rd);
mrelbp=[mrelbp_ci mrelbp_ni mrelbp_rd];
end






 








            




