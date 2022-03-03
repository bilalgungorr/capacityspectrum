import sys
import numpy as np
import matplotlib.pyplot as plt

def interpolate(x2, x0, x1, y0, y1):
    y2 = y1 + (x2 - x1)*(y1 - y0)/(x1 - x0)
    return y2


def find_value(x, xarray, yarray):
    for i in range(len(xarray) - 1):
        x1, x2 = xarray[i: i+2]
        if x1 <= x <= x2 or x1 >= x >= x2:
            s1 = slice(i, i+2)
            return interpolate(x, *xarray[s1], *yarray[s1])

def find_value2(x, xarray, yarray):

    if x <= min(xarray):
        return interpolate(x, *xarray[0: 2], *yarray[0: 2])
    elif x >= max(xarray):
        return interpolate(x, *xarray[-3: -1], *yarray[-3: -1])
    return find_value(x, xarray, yarray)


class Structure:
    def __init__(self, name, BKS):
        """
        Parameters:
            name (str): raporda kullanilmak uzere yapi adi.
            BKS (float): Bina kullanim sinifi.
            R (float): Tasiyici sistem davranis katsayisi.
            D (float): Dayanim fazlaligi katsayisi.
            soil_class (str): {'ZA', 'ZB', 'ZC', 'ZD', 'ZE'}, zemin sinifi.
        """
        self.name = name 
        self.BKS = BKS
        self.output = True
        self.SDS = None
    

    def find_SDS_SD1(self, SS, S1, soil_class):
        """
        SDS ve SD1 degerleri geri donderir.

        Parameters:
            SS (float): Kisa periyot harita spectral ivme katsayisi
            S1 (float): 1.0 saniye periyot icin harita spectral ivme katsayisi
            soil_class (str): {'ZA', 'ZB', 'ZC', 'ZD', 'ZE'}, zemin sinifi

        Referans:
            TBDY-2.3.4.1, 4A.3.2
        """
        
        self.SS = SS
        self.S1 = S1
        self.soil_class = soil_class
        

        SS_range = [0.25 , 0.50 , 0.75, 1.00 , 1.25 , 1.50 ]
        FS_table = {"ZA": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8],
                    "ZB": [0.9 , 0.9 , 0.9 , 0.9 , 0.9 , 0.9],
                    "ZC": [1.3 , 1.3 , 1.2 , 1.2 , 1.2 , 1.2],
                    "ZD": [1.6 , 1.4 , 1.2 , 1.1 , 1.0 , 1.0],
                    "ZE": [2.4 , 1.7 , 1.3 , 1.1 , 0.9 , 0.8]}


        S1_range = [0.10 , 0.20 , 0.30, 0.40 , 0.50 , 0.60 ]
        F1_table = {"ZA": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8],
                    "ZB": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8],
                    "ZC": [1.5 , 1.5 , 1.5 , 1.5 , 1.5 , 1.4],
                    "ZD": [2.4 , 2.2 , 2.0 , 1.9 , 1.8 , 1.7],
                    "ZE": [4.2 , 3.3 , 2.8 , 2.4 , 2.2 , 2.0]}
    
        
        # tasarim spektral ivmeleri
        # Kısa periyod
        FS = find_value2(SS, SS_range, FS_table[soil_class])
        self.SDS = SS*FS

        # 1sn periyod
        F1 = find_value2(S1, S1_range, F1_table[soil_class])
        self.SD1 = S1*F1

        
        if self.output:
            print(f'S_DS = {self.SDS}; S_D1 = {self.SD1}')


    def designSpectrum(self, T_list, R=None, D=None, Ra=None, includeTA_TB=False):
        """
        R ve D degerleri veya Ra degerinin verilmesi bagli olarak azaltilmis 
            yatay tasarim spectrum degerlerini geri doner. Sayet elastik 
            tasarim spectrumu elde edilmek istenirse Ra=1 degeri verilebilir.

        Parameters:
            T_list (list(float), float): spectral degerinin hesaplanmasi istenilen
                periyot degerleri.
            R (float): Tasiyici sistem davranis katsayisi.
            D (float): Dayanim fazlaligi katsayisi
            Ra (float): Deprem yuku azaltma katsayisi
        Referans:
            TBDY-2.3.4.1, 4A.3.2
        """

        BKS, SDS, SD1 = self.BKS, self.SDS, self.SD1
        self.R = R
        self.D = D

        if SDS == None:
            sys.exit('find_SDS_SD1 metodu once cagirilmalidir.')
            

        # kose periyotlar
        TA = 0.2*SD1/SDS
        TB = SD1/SDS
        TL = 6
        
        if self.output:
            print(f'T_A = {round(TA, 5)}; T_B = {round(TB, 5)}')

        # T_list verisi duzenlemesi
        if includeTA_TB:
            # sinir degerleri T_list objesinde bulundurmak icin duzenleme
            if max(T_list) > TA and TA not in T_list:
                T_list.append(TA)
            if max(T_list) > TB and TB not in T_list:
                T_list.append(TB)
            T_list = sorted(T_list)

        
        # azaltilmis veya elastik tasarim spectral ivme hesabi
        SaR = []
        for T in T_list:
            if Ra is None:
                if T <= TB:
                    Ra = D + (R/BKS - D)*T/TB
                else:
                    Ra = R/BKS
            
            if T < TA :
                Sae = (0.4 + 0.6*(T/TA))*SDS
            elif T <= TB:
                Sae = SDS
            elif T <= TL:
                Sae = SD1/T
            else:
                Sae = SD1*TL/(T**2)
        
            SaR.append(round(Sae/Ra, 5))
        
        
        return T_list, SaR


def Secant(x0, x1, y, func, error, Nmax):
    """
    Parameters:
        x0 (float): initial value
        x1 (float): initial value
        y (float): value to be found.
        func (function): a function that takes x paramater.
        Nmax: number of steps.
    """

    step = 1
    condition=True
    x2 = x0
    f0 = func(x0)

    while condition: 
        f1 = func(x1)
        if f1 != f0:
            x2 = x1 + (y - f1)*(x1 - x0)/(f1 - f0)
            x0, x1 = x1, x2
            f0 = f1
            step += 1
            condition = abs(y - f1) > error and step < Nmax
        else:
            x2 = (x0 + x1)/2
            break
    if step < Nmax - 1:
        return x2


def find_yield_disp_force(deformations, forces, keff=None):
    len_data = len(deformations)
    du = deformations[-1]
    fu = forces[-1]
    
    # area under capacity curve
    area_curve = 0.
    for i in range(len_data - 1):
        d1, d2 = deformations[i: i+2]
        f1, f2 = forces[i: i+2]
        area_curve += 0.5*(f1 + f2)*(d2 - d1)

    def area_bil_func(dy, ke):
        """Compute area under bilinear curve"""
        fy = ke*dy
        return 0.5*((fu + fy)*(du - dy) + fy*dy)
    

    # for binding
    dy_try, fy_try = 0, 0

    def function(ke):
        nonlocal dy_try, fy_try
        
        # initial deformation
        d0 = deformations[len_data//5]
        d1 = 1.2*d0
        # initial areas based on initial deformation and stiffness
        a0, a1 = area_bil_func(d0, ke), area_bil_func(d1, ke)
        
        # deformation that makes area under bilinear curve equal to area_curve
        dy_try = interpolate(area_curve, a0, a1, d0, d1)
        print(area_curve, area_bil_func(dy_try, ke))

        fy_try = ke*dy_try
        
        # finding deformation at which corresponding force 
        # equal to 60 % of yield force (fy_try)
        dy06 = find_value(0.6*fy_try, forces, deformations)
        fy06 = ke*dy06
        
        error = (fy06 - fy_try*0.6)/fy06 
        return error

    if keff is None:
        ke0 = forces[1]/deformations[1]
        ke1 = ke0*0.8
        ke = Secant(ke0, ke1, 0, function, 1e-6, 20)
    
    else:
        ke = keff
    
    function(ke)
    return dy_try, fy_try 


def nameIt(deformations, forces, M1, phiN1, uN1, gamma1, omega_eff=None,
        T_array=None, Sae_array=None, TB=None, structure=None, plot=True):
    grav = 9.806
    PI = np.pi
    keff = None
    if omega_eff is not None:
        keff = (omega_eff**2*M1)/(gamma1*phiN1)

    # capacity curve
    Sd_curve = deformations/(gamma1*phiN1)
    Sa_curve = forces/M1 
    Dy, Vy = find_yield_disp_force(deformations, forces, keff=keff)

    Sdy = Dy/(gamma1*phiN1)
    Say = Vy/M1

    weff = (Say/Sdy)**0.5
    Teff = 2*PI/weff
    

    # demand curve
    if structure is not None:
        if T_array is None:
            T_array = [round(i*0.05, 2) for i in range(81)]
        
        TB = structure.SD1/struct.SDS
        res = structure.designSpectrum(list(T_array), Ra=1, includeTA_TB=True)
        T_array, Sae_array = np.array(res)
        Sae_array = Sae_array*grav
        Sae_eff = structure.designSpectrum([Teff], Ra=1)[1][0]*grav
    
    elif T_array is not None and Sae_array is not None:
        Sae_array = np.array(Sae_array)
        T_array = np.array(T_array)
        Sae_eff = find_value(Teff, T_array, Sae_array)
    
    else: 
        sys.exit('Must be given T_array and Sae_array or structure')

    w_array = (2*PI)/T_array
    Sde_array = Sae_array/w_array**2
    
    Sde_eff = Sae_eff/weff**2

    Ry = Sae_eff/Say
    if Teff < TB:
        Cr = (1 + (Ry - 1)*TB/Teff)/Ry
    else:
        Cr = 1
    
    Mu = Cr*Sde_eff/Sdy
    Du = Cr*Sde_eff*phiN1*gamma1
    
    print(f'mu: {Mu} Ry: {Ry}')
    print(f'Sdy {Sdy} Say: {Say} Dy:{Dy}  Du:{Du} Sde_eff: {Sde_eff}\
            Sae_eff: {Sae_eff}Cr: {Cr} Teff: {Teff} weff: {weff}')
    if plot:
        # capacity curve
        plt.figure()
        plt.plot([0, Dy, deformations[-1]], [0, Vy, forces[-1]], label='bilinear')
        plt.plot(deformations, forces, '.--', label='capacity curve', ms=3)
        plt.legend()
        plt.show()

        # capacity spectrum
        plt.figure()
        plt.plot([0, Sdy, Sd_curve[-1]], [0, Say, Sa_curve[-1]], label='bilinear')
        plt.plot(Sd_curve, Sa_curve, '--', label='capacity spectrum')
        plt.plot(Sde_array, Sae_array, '--', label='demand spectrum')
        plt.plot([Sdy, Sde_eff, Sde_eff, Sde_eff*Cr, Du], [Say, Sae_eff, 0, 0, 0], '.--')
        plt.text(Du, 0, f'{Du:.4f}')
        plt.legend()
        plt.show()



def main():
    # if omega_eff is not given, then will be calculated according to ATC-40.
    omege_eff = None
    M1 = 66.14
    phiN1 = 0.156987
    uN1 = 0.1
    gamma1 = 8.1328
    displacement, force = np.loadtxt('data.txt').T # capacity curve


    #########################################################################
    # 1. if elastic spectrum is given by user:
    # -----------------------------------------------------------------------
    periods, accel = np.loadtxt('data1.txt')[:200].T # elastik spectrum
    accel = accel*9.806
    TB = 0.7072
    nameIt(displacement, force, M1, phiN1, uN1, gamma1, omega_eff=None,
            T_array=periods, Sae_array=accel, TB=TB, structure=None, plot=True)

    #########################################################################
    # 2. if elastic spectrum will computed with Structure's method.
    # -----------------------------------------------------------------------
    SS = 0.726
    S1 = 0.213
    BKS = 1
    soil_class = 'ZE'
    struct = Structure('name', BKS)
    struct.find_SDS_SD1(SS, S1, soil_class)

    # Giving T_array is up to user.
    T_array = None 
    T_array = [round(i*0.05, 2) for i in range(81)]

    nameIt(displacement, force, M1, phiN1, uN1, gamma1, omega_eff=None,
            T_array=T_array, Sae_array=None, TB=None, structure=struct, plot=True)

if __name__ == '__main__':
    main()



