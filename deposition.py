import numpy as np


# variables
## particle variables
rho = 2000 # Particle density 颗粒密度, kg m-3 ()

## environment variables
Tw = 273.15 # Sea surface temperature, K (Meguro, 2004)

## Gao给出的参量
Dd = 2.5e-6 # 干颗粒直径, m
z = 20 # Height Z, m (Gao, 2020)
rh = 0.57 # Relative humidiy, <1 
Ta = 273 # Air temperature, K
wspd = 14 # Wind speed, 在水面1Om以上的风速. m/s
# constant
nu = 1.338e-5 # m2/s
phi = 1.729e-5 # air viscosity (kg/m-s): μ 
Kc = 0.4 # Von Karman's constant, 卡门常数
k = 1.380649e-23 # J/K, Boltzmann’s constant, 玻尔兹曼常数
g = 9.8 # gravitational acceleration constant, 重力加速度常数 (m s-2)   
u_starIni = 0.02 # Initialization: Fiction velocity    
err = 0.0000001 # recursive stop condition in calc_ustar()

c_ss = {# sea salt
    "c1":0.7674, 
    "c2":3.079,
    "c3":2.573e-11,
    "c4":-1.424
}
c_ur = {# urban
    "c1":0.3926, 
    "c2":3.101,
    "c3":4.190e-11,
    "c4":-1.404
}
c_ru = {# rural
    "c1":0.2789, 
    "c2":3.115,
    "c3":5.415e-11,
    "c4":-1.399
}
c_ns = {# (NH4)2SO4
    "c1":0.4809, 
    "c2":3.082,
    "c3":3.110e-11,
    "c4":-1.428
}
pm_type_dict = {
    'ss':c_ss,
    'ur':c_ur,
    'ru':c_ru,
    'ns':c_ns
}

def calc_ustar(u_starIni,err, nu, Ta, Tw, wspd, z, g, Kc): 
    '''Calculate z0, psi_h & psi_m by using an interactive procedure (Williams, 1982)
    z0: Roughtness length
    psi_h: Stability function for momentum
    '''
    curuu = u_starIni
    preuu = 0
    while abs(curuu - preuu) > err:
        z0 = nu/9.1*curuu + 0.016*curuu**2/g
        stability = g*(Ta-Tw)*z*np.log(z/z0)/Ta*wspd**2# 'z/L'
        if stability>=0:
            psi_h = -5.2*stability
            psi_m = psi_h
        else:
            psi_h = np.exp(0.598 + 0.39*np.log(-1*stability) - 0.09*(np.log(-1*stability))**2)
            psi_m = np.exp(0.032 + 0.448*np.log(-1*stability) - 0.132*(np.log(-1*stability))**2)
        preuu = curuu
        curuu = Kc*wspd/(np.log(z/z0)-psi_m)
    return [curuu, z0, psi_h, psi_m]
    
def calc_Dw_Quinn1998(rh, rho, Dd, Ta):
    ''''wet' deposition velocity (Quinn & Ondov, 1998)
    rh:相对湿度
    '''
    rho_phys = rho * 1e-3 # kg m-3 -> g cm-3
    rho_aero = 1 # g cm-3, aerodynamic densities
    ''' physical densities
    1. only fine particles (< 1.8 um) grow, and 1.9 g cm-3 is the initial density (like sulfuric acid and ammonium sulfate particles), but the density approaches 1 g cm-3 as particle size increases.
    2. larger particles (those on stages 0, 1, and 2) do not grow, and their densities are 2.4 g cm-3
    '''
    if rh < 0.75:
        Dw_aero = 0.42
    else:
        # Dw_aero = (-0.1157/np.log(rh)+0.0203)**(1/3) # (Koutrakis et al., 1989)
        Dw_aero = Dd*((-0.1157/np.log(rh)+0.0203)/0.42)**(1/3) # (Quinn & Ondov, 1998)
    # aerodynamic to physical (Cphys, Caero: slip-correction factors)
    Cphys = cc(dp=Dw_aero, P=1000e2, Ta=Ta)
    Caero = cc(dp=Dw_aero, P=1000e2, Ta=Ta)
    Dw = Dw_aero * ((rho_phys*Cphys)/(rho_aero*Caero))**(0.5)
    return Dw
def cc(dp, P, Ta, k):
    '''Cunningham slip correction factors 
    dp: particle diameter;
    P: the system pressure;
    Ta: the temperature
    '''
    lbd = k*Ta/(2**0.5*np.pi*P*dp**2) # lambda
    '''
    lbd: molecular mean free path;
    dm: the collision diameter of the molecules;
    k: Boltzmann constant;
    '''
    Kn = 2*lbd/dp # Knudsen number, 克努森数=平均自由程/半径
    alpha = 1.257
    beta = 0.4
    gamma = 1.1
    Cc = 1 + Kn * (alpha + beta * np.exp(-1*gamma/Kn)) 
    # alpha, beta, gamma are empirical constants which depend on the gas types and particle material (Allen and Raabe, 1985)
    return Cc

def calc_Dw_Gerber1985(rh,Dd,pm_type):
    '''计算湿颗粒直径 (Gerber, 1985)
    rh:相对湿度
    Dd:干颗粒直径
    pm_type:颗粒类型，分为：
        ss: sea salt
        ur: urban
        ru: rural
        ns: (NH4)2SO4
    '''    
    var = pm_type_dict[pm_type]
    c1, c2, c3, c4 = var['c1'], var['c2'], var['c3'], var['c4']
    Dw = 2 * (c1*(Dd/2)**c2/(c3*(Dd/2)**c4-np.log10(rh)) + (Dd/2)**3) **(1/3)
    return Dw

def calc_vd():
    [u_star, z0, psi_h, psi_m] = calc_ustar()
    Dw = calc_Dw_Gerber1985(rh=rh, Dd=Dd, pm_type='ss')
    # Dw = calc_Dw_Quinn1998(rh=rh)
    # Dw = 4.5239701e-6

    a = 1.7e-6*wspd**3.75 # (Wu, 1979)
    Dc = (2.38e-7/Dw)*(1+0.163/Dw+0.0548*np.exp(-6.66*Dw)/Dw)# particle's diffusivity, 粒子扩散系数 (Davies, 1966)
    vgd = rho*g*(Dd**2)/(18*phi)# Stokes Law, m/s
    vgw = rho*g*(Dw**2)/(18*phi)# Stokes Law, m/s
    st = (vgw*u_star**2)/(g*nu) # Stokes number
    sc = nu/Dc # Schmidt number
    x = Kc*u_star/(np.log(z/z0)-psi_h) # km, k_ab, k_as. Slinn (1976) and Slinn and Slinn (1980)
    y = (u_star**2/(Kc*wspd))*(10**(-(3/st))+sc**(-0.5))# k_ss, k_bs. (Hess and Hicks, 1975)
    A = x*((1-a)*x+a*x+vgd)+(1-a)*(x+vgd)*a*(x+y+vgw) # (Williams, 1982)
    B = x*(x+y+vgw)+(1-a)*a*(x+y+vgw)**2 # (Williams, 1982)
    vd = (A/B)*((1-a)*(y+vgw)+(x*a*(y+vgw))/(x+a*(x+y+vgw)))+a**2*(y+vgw)*(x+vgd)/(x+a*(x+y+vgw)) # (Williams, 1982)
    return vd

def main():
    
    [u_star, z0, psi_h, psi_m] = calc_ustar()
    Dw = calc_Dw_Gerber1985(rh=rh, Dd=Dd, pm_type='ss')
    # Dw = calc_Dw_Quinn1998(rh=rh)
    # Dw = 4.5239701e-6

    a = 1.7e-6*wspd**3.75 # (Wu, 1979)
    Dc = (2.38e-7/Dw)*(1+0.163/Dw+0.0548*np.exp(-6.66*Dw)/Dw)# particle's diffusivity, 粒子扩散系数 (Davies, 1966)
    vgd = rho*g*(Dd**2)/(18*phi)# Stokes Law, m/s
    vgw = rho*g*(Dw**2)/(18*phi)# Stokes Law, m/s
    st = (vgw*u_star**2)/(g*nu) # Stokes number
    sc = nu/Dc # Schmidt number
    x = Kc*u_star/(np.log(z/z0)-psi_h) # km, k_ab, k_as. Slinn (1976) and Slinn and Slinn (1980)
    y = (u_star**2/(Kc*wspd))*(10**(-(3/st))+sc**(-0.5))# k_ss, k_bs. (Hess and Hicks, 1975)
    A = x*((1-a)*x+a*x+vgd)+(1-a)*(x+vgd)*a*(x+y+vgw) # (Williams, 1982)
    B = x*(x+y+vgw)+(1-a)*a*(x+y+vgw)**2 # (Williams, 1982)
    vd = (A/B)*((1-a)*(y+vgw)+(x*a*(y+vgw))/(x+a*(x+y+vgw)))+a**2*(y+vgw)*(x+vgd)/(x+a*(x+y+vgw)) # (Williams, 1982)
    
    print('Dd: {:.2e} m, Dw: {:.2e} m, Dc: {:.7f}\n'.format(Dd,Dw,Dc))
    print('vgd: {:.2e} m/s, vgw: {:.2e} m/s\n'.format(vgd,vgw))
    print('st:{:.7f}, sc:{}\n'.format(st,sc)) 
    print(('a:{:.7f}, km=k_ab=k_as:{:.7f}, k_ss=k_bs:{:.7f}\n').format(a,x,y))
    print('A:{:.7f}, B:{:.7f}\n'.format(A,B))
    print('vd:{:.7f}'.format(vd))
    
def calc_vd(u_starIni, err, nu, Ta, Tw, wspd, z, g, Kc, rh, Dd, rho, phi):
    [u_star, z0, psi_h, psi_m] = calc_ustar(u_starIni,err, nu, Ta, Tw, wspd, z, g, Kc)
    Dw = calc_Dw_Gerber1985(rh=rh, Dd=Dd, pm_type='ss')
    a = 1.7e-6*wspd**3.75 # (Wu, 1979)
    Dc = (2.38e-7/Dw)*(1+0.163/Dw+0.0548*np.exp(-6.66*Dw)/Dw)# particle's diffusivity, 粒子扩散系数 (Davies, 1966)
    vgd = rho*g*(Dd**2)/(18*phi)# Stokes Law, m/s
    vgw = rho*g*(Dw**2)/(18*phi)# Stokes Law, m/s
    st = (vgw*u_star**2)/(g*nu) # Stokes number
    sc = nu/Dc # Schmidt number
    x = Kc*u_star/(np.log(z/z0)-psi_h) # km, k_ab, k_as. Slinn (1976) and Slinn and Slinn (1980)
    y = (u_star**2/(Kc*wspd))*(10**(-(3/st))+sc**(-0.5))# k_ss, k_bs. (Hess and Hicks, 1975)
    A = x*((1-a)*x+a*x+vgd)+(1-a)*(x+vgd)*a*(x+y+vgw) # (Williams, 1982)
    B = x*(x+y+vgw)+(1-a)*a*(x+y+vgw)**2 # (Williams, 1982)
    vd = (A/B)*((1-a)*(y+vgw)+(x*a*(y+vgw))/(x+a*(x+y+vgw)))+a**2*(y+vgw)*(x+vgd)/(x+a*(x+y+vgw)) # (Williams, 1982)

if __name__ == '__main__':
    # main()
    calc_vd()