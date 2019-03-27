function SM = GetSysMat(option)
%SM 系统矩阵
%sp0 源点坐标
%dp0 探测器横纵坐标 2*Ndetector
%grid 网格参数
%theta 投影角度
os = option.Source2Center;
od = option.Source2Detector-os;
nd = option.detecNum;
H = option.detecLen/2;
mode = option.mode;
[sp0, dp0] = getgeoargs(os, od, nd, H, mode);
theta = (0:option.angIncrement:(option.angNum-1)*option.angIncrement)/360*2*pi;
grid = getgrid(option.imgSize,option.imgPysicSize/2);


SM.lsect=cell(1,1);
SM.isect=cell(1,1);
%区分扇束和平行束
if size(sp0,2)==1%扇束
    SM.mode='fan';
    for t=1:length(theta)
%          if mod(t,10)==0
%              disp(strcat('angle:',num2str(t)));
%          end;    
         sp.x=sp0(1)*cos(theta(t))+sp0(2)*sin(theta(t));
         sp.y=sp0(2)*cos(theta(t))-sp0(1)*sin(theta(t));
         for k=1:length(dp0)
             dp.x=dp0(1,k)*cos(theta(t))+dp0(2,k)*sin(theta(t));
             dp.y=dp0(2,k)*cos(theta(t))-dp0(1,k)*sin(theta(t));
             p=linesect(sp,dp,grid);
             if ~isempty(p)
                 if length(p.x)>=2
                         [lsect,isect]=sectinfo(p,grid);
                         SM.lsect{t}{k}=lsect;
                         SM.isect{t}{k}=isect;
                 else
                         SM.lsect{t}{k}=[];
                         SM.isect{t}{k}=[];
                 end;
              else
                 SM.lsect{t}{k}=[];
                 SM.isect{t}{k}=[];
              end;
         end
    end
    
else %平行束 仅在于源坐标的计算有区别
    SM.mode='par';
    for t=1:length(theta)
%          if mod(t,10)==0
%              disp(strcat('angle:',num2str(t)));
%          end;    
         for k=1:length(dp0)
             sp.x=sp0(1,k)*cos(theta(t))+sp0(2,k)*sin(theta(t));
             sp.y=sp0(2,k)*cos(theta(t))-sp0(1,k)*sin(theta(t));
             dp.x=dp0(1,k)*cos(theta(t))+dp0(2,k)*sin(theta(t));
             dp.y=dp0(2,k)*cos(theta(t))-dp0(1,k)*sin(theta(t));
             p=linesect(sp,dp,grid);
             if ~isempty(p)
                 if length(p.x)>=2
                         [lsect,isect]=sectinfo(p,grid);                        
                         SM.lsect{t}{k}=lsect;
                         SM.isect{t}{k}=isect;
                 else
                         SM.lsect{t}{k}=[];
                         SM.isect{t}{k}=[];
                 end;
              else
                 SM.lsect{t}{k}=[];
                 SM.isect{t}{k}=[];
             end;
         end
    end
end
end
function [sp0,dp0]=getgeoargs(varargin)
%getgeoargs(os,od,nd)
%getgeoargs(os,od,nd,H)
%getgeoargs(os,od,nd,H,mode)
%os:图像中心到源距离
%od:图像中心到探测器距离
%nd:探测器个数
%H:探测器半长
%mode:'fan' or 'par',缺省为扇束

if nargin>=3
    os=varargin{1};
    od=varargin{2};
    nd=varargin{3};
end
if nargin>=4
    H=varargin{4};
else
    H=(os+od)/(os-1);
end
if nargin>=5 
    mode=varargin{5};
else
    mode='fan';
end
if nargin>6  error('Invalid arguments!');end;


dL=2*H/nd;
if strcmp(mode,'fan') sp0=[-os;0];
else if strcmp(mode,'par') sp0=[-os*ones(1,nd);(-H+dL/2):dL:(H-dL/2)];
    else error('Invalid mode!');end
end;
dp0=[od*ones(1,nd);(-H+dL/2):dL:(H-dL/2)];
% disp(strcat('distance from image''s center to origin:',num2str(os),'cm'));
% disp(strcat('distance from image''s center to detector:',num2str(od),'cm'));
% disp(strcat('width of detecor:',num2str(2*H),'cm'));
end
function grid=getgrid(varargin)
%网格参数
%N 网格大小
%L 网格物理半长 默认为1
if nargin==1;
    N=varargin{1};
    grid.L=1;
else
    N=varargin{1};
    grid.L=varargin{2};
end;

% disp(strcat('size of grid:',num2str(2*grid.L),'cm'));
grid.N=N;
grid.xg=grid.L*linspace(-1,1,N+1);
grid.yg=grid.L*linspace(-1,1,N+1);
end
function p=linesect(p1,p2,grid)
%p:交点集合
%p1,p2:已知直线的两个端点
%grid:网格参数，默认为规整网格，均由水平线和竖直线构成
%     包括如下信息:grid.xg,grid.yg，分别规定了网格的x,y坐标点
%                 grid.N,规定了网格的大小

p=[];
%1 判断是否特殊情况
if (abs(p1.x-p2.x)>1e-4 && abs(p1.y-p2.y)>1e-4)
    %2.直线方程为(y-y1)/(x-x1)=(y2-y1)/(x2-x1)
    %计算其与所有水平线y=d以及所有竖直线x=c的交点
    sy=(grid.xg-p1.x)*(p2.y-p1.y)/(p2.x-p1.x)+p1.y;
    sx=(grid.yg-p1.y)/(p2.y-p1.y)*(p2.x-p1.x)+p1.x;
    p.x=[grid.xg sx];
    p.y=[sy grid.yg];
    
    %3剔除不在网格范围内的交点
    p.y(p.x<min(grid.xg) | p.x>max(grid.xg))=[];
    p.x(p.x<min(grid.xg) | p.x>max(grid.xg))=[];
    p.x(p.y<min(grid.yg) | p.y>max(grid.yg))=[];
    p.y(p.y<min(grid.yg) | p.y>max(grid.yg))=[];
    
    %剔除重复计算的交点 这点如何能方便做到
    [p.x,i]=unique(p.x);
    p.y=p.y(i);
end

%特殊情况直接返回交点
if (abs(p1.x-p2.x)<1e-4)
    if (p1.x>min(grid.xg) && p1.x<max(grid.xg))
        p.y=grid.yg;
        p.x=ones(size(p.y))*p1.x;
    end;
end
if (abs(p1.y-p2.y)<1e-4)
    if (p1.y>min(grid.yg) && p1.y<max(grid.yg))
        p.x=grid.xg;
        p.y=ones(size(p.x))*p1.y;
    end;
end
end
function [lsect,isect]=sectinfo(p,grid)
%由交点信息提取出交线长度和交线所在的像素值
%p:交点列表
%image:像素值 
%grid:网格参数，与linesect的相同
%lsect:交线长度
%isect:交线所在的像素对应的一维id

np=length(p.x);
%1.将交点排序
[ps,ix]=sort(p.x);
p.x=ps;
p.y=p.y(ix);

%2.计算交线长度
lsect=sqrt(diff(p.x).^2+diff(p.y).^2);

%3.计算中点坐标
cp.x=(p.x(1:np-1)+p.x(2:np))/2;
cp.y=(p.y(1:np-1)+p.y(2:np))/2;

%4.由中点坐标计算交线所属的像素id
isect=xy2id(cp,grid);

end
function Idx=xy2id(p,grid)
%根据点p的xy坐标转换成图像的像素id
%注意y坐标和像素id排列方向相反

x0=min(grid.xg);y0=min(grid.yg);
dx=grid.xg(2)-grid.xg(1);
dy=grid.yg(2)-grid.yg(1);
n=1+floor((p.x-x0)/dx);%列
n(n>grid.N)=grid.N;
m=grid.N-floor((p.y-y0)/dy);%行
m(m<1)=1;
Idx=sub2ind([grid.N,grid.N],m,n);
end