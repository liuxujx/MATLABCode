function [rho, rho_a] = Density(altitude,auxdata,planetary,model)
% 地球：基于 1976 US 标准大气模型拟合大气密度
% 火星：参考文献：Advances in sapcecraft atmospheric entry gudiance, P50-51
% altitude = 高度 m (0 < alt < 130e3)
% rho = 大气密度 kg/m^3
% rho_a = d(rho)/d(alt) (kg/m^3/m)
% planetary 1=地球 2=火星
% model 1=一阶指数模型 2=高阶指数模型

alt = altitude;
rho0 = auxdata.rho0;
hs = auxdata.hs;
if planetary==1 && model==1
    rho = rho0*exp(-alt/hs);
    rho_a = -1/hs*rho;
elseif planetary==1 && model==2
    a= -0.009334500409893852;
    b= -0.07713480853494097;
    c= -0.00382948508292597;
    d= 5.242803244759748E-05;
    e= 6.454862591920205E-07;
    f= -2.031133609734722E-08;
    g= 1.568378909033718E-10;
    h= -3.928350728483702E-13;
    alt= altitude/1000;
    if alt >130.0
        alt = 130.;
    end
    y = a+alt.*(b+alt.*(c+alt.*(d+alt.*(e+alt.*(f+alt.*(g+alt.*h))))));
    rho = rho0*exp(y);
    dy = (b+alt.*(2.*c+alt.*(3*d+alt.*(4*e+alt.*(5*f+alt.*(6*g+7*alt.*h))))));
    rho_a= rho.*dy;
elseif planetary==2 && model==1
    rho = rho0*exp(-alt/hs);
    rho_a = -1/hs*rho;
elseif planetary==2 && model==2
    if alt<=10e3 || alt>=140e3
        warning('火星大气密度模型在[10~140]km之外不准确');
    end
    b0 = -4.324;
    b1 = -9.204e-5;
    b2 = -1.936e-11;
    b3 = -7.507e-15;
    b4 = -4.195e-20;
    rho = b0+b1*alt+b2*alt.^2+b3*alt.^3+b4*alt.^4;
    rho_a = b1+2*b2*alt+3*b3*alt.^2+4*b4*alt.^3;
end

end