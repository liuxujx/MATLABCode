function [outp]=Atmosphere(altitude,planetary)
%根据地势高度计算大气
% 输入高度 [m]
% 输出大气温度[K]、声速[m/s]、气压[Pa]、密度[kg/m^3]
%参考文献：[1] Sub_Height2AirQuality.m [2]一种火星大气密度三维解析模型,表1
% [3] Advances in sapcecraft atmospheric entry gudiance, P50-51
if planetary==1 % 1=地球 2=火星
    H = altitude;
    R0=6356766;
    Z=R0*H/(R0+H); %几何高度
    if (Z<174.9)
        warning('几何高度最低为176米');
        AtmosTemp=0;
        AtmosPress=0;
        AtmosDensity=0;
        SoundSpeed=0;
    elseif (174.9<=Z)&&(Z<11000)
        AtmosTemp=288.15-0.0065*Z;
        AtmosPress=1.01325e+5*power(1-0.225577e-4*Z,5.25588);
        AtmosDensity=1.225*power(1-0.225577e-4*Z,4.25588);
        SoundSpeed=20.0468*sqrt(AtmosTemp);
    elseif (11000<=Z)&&(Z<20000)
        AtmosTemp=216.65;
        AtmosPress=2.263204e+4*exp(-1.576885e-4*(Z-11000));
        AtmosDensity=0.3639176*exp(-1.576885e-4*(Z-11000));
        SoundSpeed=20.0468*sqrt(AtmosTemp);
    elseif (20000<=Z)&&(Z<32000)
        AtmosTemp=216.65+0.001*(Z-20000);
        AtmosPress=5.474849e+3*power(1+4.615740e-6*(Z-20000),-34.16322);
        AtmosDensity=8.803471e-2*power(1+4.615740e-6*(Z-20000),-35.16322);
        SoundSpeed=20.0468*sqrt(AtmosTemp);
    elseif (32000<=Z)&&(Z<47000)
        AtmosTemp=228.65+0.0028*(Z-32000);
        AtmosPress=8.680160e+2*power(1+1.224579e-5*(Z-32000),-12.20115);
        AtmosDensity=1.322497e-2*power(1+1.224579e-5*(Z-32000),-13.20115);
        SoundSpeed=20.0468*sqrt(AtmosTemp);
    elseif (47000<=Z)&&(Z<51000)
        AtmosTemp=270.65;
        AtmosPress=1.109058e+2*exp(-1.262266e-4*(Z-47000));
        AtmosDensity=1.427527e-3*exp(-1.262266e-4*(Z-47000));
        SoundSpeed=20.0468*sqrt(AtmosTemp);
    elseif (51000<=Z)&&(Z<71000)
        AtmosTemp=270.65-0.0028*(Z-51000);
        AtmosPress=66.93853*power(1-1.034546e-5*(Z-51000),12.20115);
        AtmosDensity=8.616011e-4*power(1-1.034546e-5*(Z-51000),11.20115);
        SoundSpeed=20.0468*sqrt(AtmosTemp);
    elseif (71000<=Z)&&(Z<84852)
        AtmosTemp=214.65-0.002*(Z-71000);
        AtmosPress=3.956392*power(1-9.317494e-6*(Z-71000),17.08161);
        AtmosDensity=6.421057e-5*power(1-9.317494e-6*(Z-71000),16.08161);
        SoundSpeed=20.0468*sqrt(AtmosTemp);
    else
        if (84852<=H)&&(H<91000)
            AtmosTemp=186.87;
            AtmosPress=1.01325e-1*(2.2730+1.042e-6*H)*exp((87284.8-H)/5470);
            AtmosDensity=4.4603475e-6*exp((87284.8-H)/5470);
            SoundSpeed=275.7302;
        elseif (91000<=H)&&(H<110000)
            AtmosTemp=263.1905-76.3232*sqrt(1-power((H-91000)/(-19942.9),2));
            AtmosPress=1.01325e-1*(2.2730+1.042e-6*H)*exp((87284.8-H)/5470);
            AtmosDensity=4.4603475e-6*exp((87284.8-H)/5470);
            SoundSpeed=275.7302;
        elseif (110000<=H)&&(H<120000)
            AtmosTemp=240.0+0.012*(H-110000);
            AtmosPress=1.01325e-1*(2.2730+1.042e-6*H)*exp((87284.8-H)/5470);
            AtmosDensity=4.4603475e-6*exp((87284.8-H)/5470);
            SoundSpeed=275.7302;
        elseif (120000<=H)&&(H<1000000)
            AtmosTemp=1000-640.0*exp(-0.000018758*(H-120000)*(R0+120000)/(R0+H));
            AtmosPress=1.01325e-1*(2.2730+1.042e-6*H)*exp((87284.8-H)/5470);
            AtmosDensity=4.4603475e-6*exp((87284.8-H)/5470);
            SoundSpeed=275.7302;
        else
            error('高度必须小于1000km');
        end
    end
elseif planetary==2 %文献[2]表1日本大学拟合数据
    H = altitude;
    if H<7000
        AtmosTemp=241.0-0.999*(H/1000);
    else
        AtmosTemp=241.0-2.22*(H/1000);
    end
    if H<=10e3 || H>=120e3
        warning('火星大气声速模型在[10~120]km之外不准确');
    end
    %计算声速参考 Advances in sapcecraft atmospheric entry gudiance, P50-51
    SoundSpeed=223.8-0.0002004*H-1.5888e-8*H^2+1.404e-13*H^3;
    AtmosPress=700*exp(-0.09*(H/1000));
    AtmosDensity=AtmosPress./(188.95110711075*AtmosTemp);
else
    error('planetary=1 地球，planetary=2 火星');
end
outp=[AtmosTemp, SoundSpeed, AtmosPress, AtmosDensity];
end