import numpy as N

def check_valid(prop):
    def wrapper_check_valid(self, T):
        if hasattr(T,'__len__'):
            Tlow = (T<self.Tmin).any()
            Thigh = (T>self.Tmax).any()
        else:
            Tlow = (T<self.Tmin)
            Thigh = (T>self.Tmax)
        if Tlow or Thigh:
            return N.nan
        else:
            return prop(self, T)
    return wrapper_check_valid


class Pipe_material():
    def __init__(self, Tmin, Tmax):
        self.Tmin = Tmin
        self.Tmax = Tmax

class SS304(Pipe_material):
    '''
    Stainless steel 304 properties. So far only thermal conductivity extracted from the data in: http://www.inductor-jmag.ru/files/content/a129160.pdf
    '''
    def __init__(self, Tmin=273., Tmax=1600.):
        Pipe_material.__init__(self, Tmin, Tmax)

    @check_valid
    def k(self, T):
        return -2e-06*T**2. + 0.0176*T + 9.8445

class SS316(Pipe_material):
    '''
    Stainless steel 316 properties. So far only thermal conductivity extracted from the data in: http://www.inductor-jmag.ru/files/content/a129160.pdf
    '''
    def __init__(self, Tmin=273., Tmax=1600.):
        Pipe_material.__init__(self, Tmin, Tmax)
    @check_valid
    def k(self, T):
        return -2E-06*T**2. + 0.0179*T + 8.3005


class Inconel601(Pipe_material):
    '''
    Inconel 601 properties. All units SI. Source: http://www.specialmetals.com/assets/smc/documents/alloys/inconel/inconel-alloy-601.pdf
    '''
    def __init__(self, Tmin=273., Tmax=1273.15):
        Pipe_material.__init__(self, Tmin, Tmax)

    @check_valid
    def k(self, T):
        return 1.686e-2*T+6.342

    def rho(self):
        return 8.11e3

    @check_valid
    def E(self, T):
        return -4.983e4*T**2.-3.264e6*T+2.105e11

    @check_valid
    def nu(self, T):
        return 7.194e-11*T**3.-8.477e-8*T**2.+1.073e-4*T+2.462e-1

    @check_valid
    def alpha(self, T):
        return 2.011e-12*T**2.+1.023e-9*T+1.325e-5

class Haynes230(Pipe_material):
    '''
    Haynes 230 properties from manufacturer data. All units SI.
    '''
    def __init__(self, Tmin=273., Tmax=1171.15):
        Pipe_material.__init__(self, Tmin, Tmax)

    @check_valid
    def k(self, T):
        return 1.996e-2*T+2.9807

    @check_valid
    def rho(self):
        return 9.05e3

    @check_valid
    def E(self, T):
        return 6.6039*T**3.-2.6694e4*T**2.-3.4176e7*T+2.2319e11

    @check_valid 
    def nu(self, T):
        return 2.0451e-14*T**4.-6.454e-12*T**3.-3.9514e-8*T**2.+5.9632e-5*T+0.30569

    @check_valid
    def alpha(self, T):
        return 9.848e-13*T**2.+2.179e-9*T+1.175e-5

    @check_valid
    def UTS(self, T):
        T0 = N.logical_and(T>=273.15, T<573.)
        T1 = N.logical_and(T>=573., T<850.)
        T2 = N.logical_and(T>=850., T<1173.15)
        if hasattr(T,'__len__'):
            UTS = N.zeros(len(T))
            UTS[T0] = 1e6*(-2.822e-7*T[T0]**3.+3.390e-4*T[T0]**2.-1.339e-1*T[T0]+2.244e2)
            UTS[T1] = 1e6*(-3.994e-8*T[T1]**4.+1.088e-4*T[T1]**3.-1.102e-1*T[T1]**2.+4.913e1*T[T1]-7.927e3)
            UTS[T2] = 1e6*(-3.381e-6*T[T2]**3.+1.177e-2*T[T2]**2.-1.388e1*T[T2]+5.549e3)

        else:
            if T0:
                UTS = 1e6*(-2.822e-7*T**3.+3.390e-4*T**2.-1.339e-1*T+2.244e2)
            elif T1:
                UTS = 1e6*(-4.183e-12*T**6.+1.755e-8*T**5.-3.058e-5*T**4.+2.834e-2*T**3.-1.473e1*T**2.+4.069e3*T-4.667e5)
            elif T2:
                UTS = 1e6*(-3.381e-6*T**3.+1.177e-2*T**2.-1.388e1*T+5.547e3)
        return UTS

class Inconel617(Pipe_material):
    '''
    Inconel 617 properties from manufacturer data (Haynes). All units SI.
    '''
    def __init__(self, Tmin=273., Tmax=1293.15):
        Pipe_material.__init__(self, Tmin, Tmax)

    @check_valid
    def k(self, T):
        return 1.553e-2*T+8.873

    def rho(self):
        return 8.36e3

    @check_valid
    def E(self, T):
        return -5.3770*T**3.-8.6056e3*T**2.-4.8497e7*T+2.2597e11

    @check_valid 
    def nu(self, T):
        return 1.3022e-18*T**6 - 5.5851e-15*T**5 + 9.6060e-12*T**4 - 8.4475e-09*T**3 + 3.9905e-06*T**2 - 9.5677e-04*T + 3.9074e-01

    @check_valid
    def alpha(self, T):
        return 2.431e-22*T**6 - 1.207e-18*T**5 + 2.394e-15*T**4 - 2.407e-12*T**3 + 1.281e-09*T**2 - 3.335e-07*T + 4.438e-05


class Inconel625(Pipe_material):
    '''
    Inconel 625 properties from manufacturer data (Haynes). All units SI.
    '''
    def __init__(self, Tmin=273., Tmax=1293.15):
        Pipe_material.__init__(self, Tmin, Tmax)

    @check_valid
    def k(self, T):
        return 6.576E-09*T**3 - 1.244E-05*T**2 + 2.196E-02*T + 4.195E+00

    def rho(self):
        return 8.44e3

    @check_valid
    def E(self, T):
        return -4.7974E+01*T**3 + 6.4592E+04*T**2 - 8.3741E+07*T + 2.2763E+11

    @check_valid 
    def nu(self, T):
        return 1.7388E-15*T**5 - 7.0343E-12*T**4 + 1.0599E-08*T**3 - 7.3294E-06*T**2 + 2.3685E-03*T - 6.6376E-03

    @check_valid
    def alpha(self, T):
        return 2.955E-12*T**2 + 2.874E-10*T + 1.228E-05

class Inconel740H(Pipe_material):
    '''
    Inconel 740H properties from https://inldigitallibrary.inl.gov/sites/sti/sti/Sort_19892.pdf: INL report: Creep-Fatigue Behavior and Damage Accumulation of a Candidate Structural Material for Concentrating Solar Power Thermal Receiver Quarter 6 Report.
    '''
    def __init__(self, Tmin=273., Tmax=1193.15):
        Pipe_material.__init__(self, Tmin, Tmax)

    @check_valid
    def k(self, T):
        return 9.487E-09*T**3 - 1.810E-05*T**2 + 2.490E-02*T + 4.307E+00

    def rho(self):
        return 8.05e3

    @check_valid
    def E(self, T):
        return -5.1309*T**3 - 1.6160E+04*T**2 - 3.6773E+07*T + 2.3363E+11

    @check_valid  
    def nu(self, T):
        return N.ones(len(T))*0.31

    @check_valid
    def alpha(self, T):
        return 1.194e-14*T**3-2.359E-11*T**2 + 1.888E-8*T + 8.603E-06


if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    plt.figure()
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    metals = [Haynes230(), Inconel617(), Inconel625(), Inconel740H()]
    plt.subplot(221)
    for mat in metals:
        Ts = N.linspace(mat.Tmin, mat.Tmax, 100)
        plt.plot(Ts-273.15, mat.k(Ts), label=mat.__class__.__name__)
    plt.ylabel(r'k (W.m${^{-1}}$.K${^{-1}}$)')
    plt.xlabel('T (C)')
    plt.subplot(222)
    for mat in metals:
        Ts = N.linspace(mat.Tmin, mat.Tmax, 100)
        plt.plot(Ts-273.15, mat.E(Ts)/1e9, label=mat.__class__.__name__)
    plt.legend()
    plt.ylabel('E (GPa)')
    plt.xlabel('T (C)')
    plt.subplot(223)
    for mat in metals:
        Ts = N.linspace(mat.Tmin, mat.Tmax, 100)
        plt.plot(Ts-273.15, mat.alpha(Ts)*1e6)
    plt.ylabel(r'${\alpha}$(${\mu}$m.m${^{-1}}$)')
    plt.xlabel('T (C)')
    plt.subplot(224)
    for mat in metals:
        Ts = N.linspace(mat.Tmin, mat.Tmax, 100)
        plt.plot(Ts-273.15, mat.nu(Ts))
    plt.xlabel('T (C)')
    plt.ylabel('Poisson')
    
    plt.savefig(sys.path[0]+'/Alloy_props.png')
    
