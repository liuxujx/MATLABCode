function [outp]  =  Gravity(inp,auxdata)
% 计算引力加速度
% 输入状态量、参数
% 输出引力加速度在位置矢量和纬度方向上的分量
% 参考文献：远程火箭飞行动力学与制导，P40-48
gE = auxdata.gE; %地球表面引力加速度 [m/s^2]
GM = auxdata.GM; %引力常数 [m^3/s^2]  
R0 = auxdata.R0; %球形引力场平均半径 [m]
g0 = auxdata.g0; %表面引力加速度 [m/s^2]
a0 = auxdata.a0; %赤道平均半径 [m]
J2 = auxdata.J2; %2阶带谐系数 无单位

Rscale = auxdata.Rscale;
Vscale = auxdata.Vscale;
ascale = auxdata.ascale;
tscale = auxdata.tscale;

scale = auxdata.scale; %是否无量纲化标志
if scale == 0
    r = inp(1);
    theta = inp(2);
    phi = inp(3);
    V = inp(4);
    gamma = inp(5);
    psi = inp(6);
elseif scale == 1 %还原为有量纲变量
    r = inp(1)*Rscale;
    theta = inp(2);
    phi = inp(3);
    V = inp(4)*Vscale;
    gamma = inp(5);
    psi = inp(6);
else
    error('无量纲化标志 scale=0 / 1');
end

%理想球体引力在位置矢量方向上的大小，没有方向
gr_i = GM/r^2;
%J2项引力在位置矢量方向上的大小，没有方向
gr_J2 = GM/r^2*(3/2*J2*(a0/r)^2*(1-3*sin(phi)^2));
%J2项引力在纬度方向上的大小，没有方向
gphi_J2 = GM/r^2*(3*J2*(a0/r)^2*sin(phi)*cos(phi));

gr = gr_i+gr_J2;
gphi = gphi_J2;

outp = {gr,gphi};
end

