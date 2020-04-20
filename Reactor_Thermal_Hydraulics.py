#调用数学计算模块
import math
#物性参数计算函数
def tgetVolum1(t):
    waterVolum = 2e-12 * t**3 + 3e-9 * t**2 - 2e-7 * t + 0.001
    return waterVolum

def tgetVolum2(t):
    gasVolumn = -1e-5 * t**3 + 0.0041 * t**2 - 0.573 * t + 27.922
    return gasolumn

def tgetHan1(t):
    waterHan = 1e-5 * t**3 - 0.0029 * t**2 + 4.47 * t - 9.8432
    return waterHan

def tgetHan2(t):
    gasHan = -3e-6 * t**3 - 0.0023 * t**2 + 2.1172 * t + 2489.3
    return gasHan

def tgetPressure(t):
    satPres = 5e-7 * t**3 - 9e-5 * t**2 + 0.007 * t - 0.1992
    return satPres

def tgetDensity(t):
    waterDensity = 2e-6 * t**3 - 0.003 * t**2 - 0.1932 * t + 1005.1
    return waterDensity

def tgetViscous1(t):
    waterViscous = -1e-5 * t**3 + 0.0065 * t**2 - 1.1997 * t + 93.887
    return waterViscous*1e-8

def tgetViscous2(t):
    waterDynavis = -2e-6 * t**3 + 0.0008 * t**2 - 0.1359 * t + 10.48002
    return waterDynavis*1e-4

def tgetViscous3(t):
    gasViscous = -2e-8 * t**3 + 5e-6 * t**2 + 0.0029 * t + 0.9047
    return gasViscous*1e-5

def tgetPrant(t):
    waterPrant = -2e-6 * t**3 + 0.0007 * t**2 - 0.1058 * t + 7.37089
    return waterPrant

def tgetConduct(t):
    waterConduct = 4e-8 * t**3 - 2e-5 * t**2 + 0.0028 * t + 0.55912
    return waterConduct

def pgetVolum1(p):
    waterVolum = 0.0011 * p**3 - 0.0008 * p**2 + 0.0003 * p + 0.00102
    return waterVolum

def pgetVolum2(p):
    gasVolum = -351.36 * p**3 + 210.31 * p**2 - 46.578 * p + 4.5985
    return gasVolum

def pgetHan1(p):
    waterHan = 10031 * p**3 - 6920.9 * p**2 + 2246.7 * p + 252.07
    return waterHan

def pgetHan2(p):
    gasHan = 4246.3 * p**3 - 2919.1 * p**2 + 892.38 * p + 2610.9
    return gasHan

def pgetViscous1(p):
    waterViscous = -126.59 * p**3 - 80.335 * p**2 + 20.354 * p + 4.1855
    return waterViscous*1e-4

def pgetViscous2(p):
    gasViscous = 8.3936 * p**3 - 5.7858 * p**2 + 1.8634 * p + 1.085
    return gasViscous*1e-5

#定义计算函数，（传热系数、速度、雷诺数、摩擦系数）
def convertCoeffi(t):
    Tb = 128
    P = 0.43
    Wnc = 0.46
    beta = 0.25
    Hout = (Tb - t)**-0.6 * (10189.3 + 90416.4 * P - (4314.4 + 46537 * P) * math.log10(100 * Wnc)) * (1.102-1.165 * beta)
    return Hout

def velocity(d,t,n):
    velocity = 4 * (M / n) / (math.pi * d**2 * tgetDensity(t))
    return velocity

def renuoNum(d,t,n):
    Re = velocity(d, t, n) * d / tgetViscous1(t)
    return Re

def mochaXishu(a):
    if a>4000 and a<1e5:
        L = 0.3164 * a**-0.25
    elif a>1e5 and a<3e6:
        L = 0.0032 + 0.221 * a** - 0.237
    else:
        L=a/64
    return L

#输入初始参数
M = eval(input('请输入初始流量(kg/s):'))
Tin = 100
Cp = 4.2211e3
Te = 128
N = 46
D0 = 0.309
D1 = 0.1
D2 = 0.108
s = 0.03
Hw = 4
Hd2 = 6
pi = math.pi
A = pi * D2 * Hd2
conductEffi = 22.2 #W/(m·k)

#请输入初始流量(kg/s):60
#Heat Exchanger Caculate 换热器计算
Tout = 106
Ta = (Tout + Tin) / 2
Tw = 119
Q1 = M * Cp * (Tout - Tin)
while True:
    Hout = convertCoeffi(Tw)
    Tw0 = Te - Q1 / (A * N * Hout)
    if abs(Tw0 - Tw) <= 0.1:
        break
    elif Tw0 > Tw:
        Tw  += 0.1
        continue
    else:
        Tw -= 0.1
        continue
Hin = tgetConduct(Ta) * 0.023 * pow(renuoNum(D1, Ta, N), 0.8) * pow(tgetPrant(Ta), 0.4) / D1
K = 1 / (D2 / (Hin * D1) + D2 * math.log(D2/D1) * 0.5 / conductEffi + 1 / Hout)
Q2 = N * A * K * (Te - Ta)
if Q1 > 1.1 and Q1 < 1.7 and abs(Q1 - Q2) <= 0.5:
    print(1)
    #fluidCalculate()
#else:continue
print(Hin, Hout)
print(Q1, Q2, K, Tw)

#1310.003536202410（Hin）1949.3496302860267（Hout）
#1519596.000000000(Q1) 1535540.5111285856(Q2) 655.9014027481747(K) 119.59999999999997(Tw)

#下降段
Rou1 = tgetDensity(Tin)
u1 = velocity(D0,Tin,1)
P0 = 101325
g = 9.8
Hd1 = 19
Ld1 = 25.8
Ew = 0.294
Ein = 0.5
Pg1 = Rou1 * g * Hd1
Pf1 = mochaXishu(renuoNum(D0, Tin, 1)) * Ld1 * Rou1 * u1**2 / (D0 * 2)
Pw1 = (6 * Ew + Ein) * Rou1 * u1**2 / 2

#换热器段
Rou2 = tgetDensity(Ta)
u2 = velocity(D1,Ta,N)
Pg2 = Rou2*g*Hd2
Pf2 = N * mochaXishu(renuoNum(D1, Ta, N)) * Hd2 * Rou2 * u2**2 / (D1 * 2)
A1 = pi * D0**2 / 4
A2 = D0 * (D2 * N + s * (N - 1)) + 0.1
a = A1 / A2
Pin2 = Rou2 * u1**2 * 0.5 * (1 - A1 / A2)**2

#上升段
Rou3 = tgetDensity(Tout)
u3 = u1
E = -2.7196 * a**5 + 4.3123 * a**4 - 2.0665 * a**3 + 0.1614 * a**2 - 0.3431 * a + 0.5034
Po = E * Rou3 * u3**2 / 2
G = M * 4 / (pi * D0**2)
Dp = mochaXishu(renuoNum(D0, Tout, 1)) * G**2 / (D0 * Rou3 * 2) + Rou3 * g
Pout = 1e5 + Pg1 - Pf1 - Pw1 - Pg2 - Pf2 - Pin2
Ps = tgetPressure(Tout) * 1e6
Zs0 = (Pout - Ps) / Dp
if Zs0 <= 5:
    Pw3 = 0
    print(1)
elif Zs0 < 8:
    Pw3 = Ew * Rou3 * u3**2
    print(2)
elif Zs0 < 11:
    Pw3=2*Ew*Rou3*u3**2
    print(3)
Dpf = mochaXishu(renuoNum(D0, Tout, 1)) * Rou3 * u3**2 / (D0 * 2)
Zs1 = (Pout - Dp * 8 - Pw3 - Dpf * 6 - Ps) / Dp + 8
print(Pg2, Pf2, Pin2, Dp, Pout)
print(Zs0, Zs1, mochaXishu(renuoNum(D0, Tin, 1)))

#56186.8606152(Pg2) 808.2491768238833(Pf2) 309.4936182778883(Pin) 9355.084871887515(Dp) 219946.0402118641(Pout)
#9.928080983099077(Zs0) 9.87813782983087(Zs1) 0.011876296865665317(mochaXishu)

#计算两相段A 段
H0 = pgetHan1(Ps / 1e6)
Hda = 11 - Zs1
Pae = 0.117039049
V3 = pgetVolum1(Ps / 1e6)
V4 = pgetVolum1(Pae)
v3 = pgetVolum2(Ps / 1e6)
v4 = pgetVolum2(Pae)
H4 = pgetHan1(Pae)
h4 = pgetHan2(Pae)
x1 = (H0 - H4) / (h4 - H4)
Rou4 = (Rou3 + 1 / (x1 * v4 + (1 - x1) * V4)) / 2
Pg4 = Rou4 * g * Hda
Fio = (1 + 0.5 * x1 * ((v3 + v4) / (V3 + V4) - 1))* (1 + 0.5 * x1 * ((pgetViscous1(Ps) + pgetViscous1(Pae)) / (pgetViscous2(Ps) + pgetViscous2(Pae)) - 1))**-0.25
Dpf4 = 0.3164 * (0.5 * G * D0 / (pgetViscous1(Ps / 1e6) + pgetViscous1(Pae)))**-0.25 * G**2 * (V3 + V4) * 0.5
Pf4 = Dpf4 * Fio * Hda
Pa4 = G**2 * x1 * ((v3 + v4) / 2 - (V3 + V4) / 2)
print(Hda, H0, H4, h4, x1, Rou4)
print(Pg4, Fio, Pf4, Pa4)
print(Pg4 + Pf4 + Pa4 - (Ps - 0.117e6))

#1.1218621701691305(Hda) 446.3872624209623(H0) 436.30005258423967(H4) 2682.164810047373(h4) 0.004491459159863608(x1) 542.2691232837783(Rou4)
#5961.841911535943(Pg4) 4.464227090969963(Fio) 57.29196103577298（Pf4） 4050.2049669631506（Pa4） 
#1.3388395350702922

Dpw4 = Ew * V4 * G**2 / 2
Ds = 1.1 / 3
Pw4 = Dpw4 * (1 + (v4 / V4 - 1) * (2 * x1 * (1 - x1) * Ds / Ew + x1))
Pc = Ps - Pg4 - Pf4 - Pa4 - 2 * Pw4
print(Dpw4, Pw4, Ps, Pc)

#98.42451196979871 2253.056546477635 127067.9999999998 112492.54806750965

V5 = pgetVolum1(Pc / 1e6)
V6 = pgetVolum1(P0 / 1e6)
v5 = pgetVolum2(Pc / 1e6)
v6 = pgetVolum2(P0 / 1e6)
H5 = pgetHan1(Pc / 1e6)
h5 = pgetHan2(Pc / 1e6)
H6 = pgetHan1(P0 / 1e6)
h6 = pgetHan2(P0 / 1e6)
x2 = (H0 - H5) / (h5 - H5)
x3 = (H0 - H6) / (h6 - H6)
xa = (x2 + x3) / 2
Rou5 = 1 / (x2 * v5 + (1 - x2) * V5)
Rou6 = 1 / (x3 * v6 + (1 - x3) * V6)
Pg5 = (Rou5 + Rou6) / 2 * g * 2
Fio = (1 + 0.5 * xa * ((v3 + v4) / (V3 + V4) - 1)) * (1 + 0.5 * xa * ((pgetViscous1(Pc / 1e6) + pgetViscous1(P0 / 1e6)) / (pgetViscous2(Pc / 1e6) + pgetViscous2(P0 / 1e6)) - 1))**-0.25
Dpf5 = 0.3164 * (0.5 * G * D0 / (pgetViscous1(Pc / 1e6) + pgetViscous1(P0 / 1e6)))**-0.25 * G**2 * (V3 + V4) * 0.5
Pf5 = Dpf5 * Fio * 2
print(x2, x3)
print(Pg5, Fio, Dpf5, Pf5)

#0.006617371300735808(x2) 0.012093269771848688(x3)
#1344.1177515873962(Pg5) 6.975137068659573(Fio) 11.464122492968384(Dpf5) 159.92765152071556(Pf5)

beta5 = 1 / (1 + (1 - x2) * V5 / (x2 * v5))
beta6 = 1 / (1 + (1 - x3) * V6 / (x3 * v6))
Pa5 = G**2 * (((1 - x3)**2 * V6 / (1 - beta6) + x3**2 * v6 / beta6) - ((1 - x2)**2 * V5 / (1 - beta5) + x2**2 * v5 / beta5))
print(beta5, beta6, Pa5)
# 0.9064354302344891 0.9515160195530594 6505.266032436612
print('Dpr=', Pg5 + Pf5 + Pa5)
print('Dpi=', (Pc - P0))
print(Pg5 + Pf5 + Pa5 - (Pc - P0))
#Dpr= 8009.311435544723 (Dpr)
#Dpi= 11167.548067509648 (Dpi)
#-3158.2366319649245(Dpr-Dpi)