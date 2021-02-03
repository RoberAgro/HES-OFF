''' Module with classes and methods to evaluate the electrical grid '''

import numpy as np


class ES:
    '''
    The ES object defines and stores parameters of an energy storage system

    Parameters
    ----------
    powerCharge : float
        Rated power for charging in Watts.
    powerDisch : float
        Rated power for discharging in Watts.        
    rocPCharge : float
        Maximum rate of change of power for charging in Watts per second.
    rocPDisch : float
        Maximum rate of change of power for discharging in Watts per second.        
            
    '''
    def __init__(self, powerCharge, powerDisch, rocPCharge, rocPDisch):
        self.__initComplete = False
        self.powerCharge = powerCharge
        self.powerDisch = powerDisch
        self.rocPCharge = rocPCharge
        self.rocPDisch = rocPDisch

        self.__initComplete = True

    @property
    def powerCharge(self):
        return self.__powerCharge
    @powerCharge.setter
    def powerCharge(self, powerCharge):
        if powerCharge < 0:
            raise ValueError("ES powerCharge must be positive.")
        self.__powerCharge = powerCharge

    @property
    def powerDisch(self):
        return self.__powerDisch
    @powerDisch.setter
    def powerDisch(self, powerDisch):
        if powerDisch < 0:
            raise ValueError("ES powerDisch must be positive.")
        self.__powerDisch = powerDisch

    @property
    def rocPCharge(self):
        return self.__rocPCharge
    @rocPCharge.setter
    def rocPCharge(self, rocPCharge):
        if rocPCharge < 0:
            raise ValueError("ES rocPCharge must be positive.")
        self.__rocPCharge = rocPCharge

    @property
    def rocPDisch(self):
        return self.__rocPDisch
    @rocPDisch.setter
    def rocPDisch(self, rocPDisch):
        if rocPDisch < 0:
            raise ValueError("ES rocPDisch must be positive.")
        self.__rocPDisch = rocPDisch


class LD:
    '''
    The LD object defines and stores parameters of electric loads

    Parameters
    ----------
    continuous : float
        Installed power of continuous loads in Watts.
    flexible : float
        Installed power of flexible loads in Watts.        
    largest : float
        Installed power of the largest load in Watts.        
    powerIn : float
        Current active power demand in Watts.                    
    '''
    def __init__(self, continuous, flexible, largest, powerIn=0.0):
        self.__initComplete = False
        self.continuous = continuous
        self.flexible = flexible
        self.largest = largest
        self.powerIn = powerIn

        self.__initComplete = True

    @property
    def continuous(self):
        return self.__continuous
    @continuous.setter
    def continuous(self, continuous):
        if continuous < 0:
            raise ValueError("LD continuous cannot be negative.")
        self.__continuous = continuous

    @property
    def flexible(self):
        return self.__flexible
    @flexible.setter
    def flexible(self, flexible):
        if flexible < 0:
            raise ValueError("LD flexible cannot be negative.")
        self.__flexible = flexible

    @property
    def largest(self):
        return self.__largest
    @largest.setter
    def largest(self, largest):
        if largest < 0:
            raise ValueError("LD largest cannot be negative.")
        elif largest > max(self.continuous, self.flexible):
            raise ValueError("LD largest is invalid.")
        self.__largest = largest

    @property
    def powerIn(self):
        return self.__powerIn
    @powerIn.setter
    def powerIn(self, powerIn):
        if powerIn < 0:
            raise ValueError("LD powerIn cannot be negative.")
        elif powerIn > (self.continuous + self.flexible):
            raise ValueError("LD powerIn is invalid.")
        self.__powerIn = powerIn


class GT:
    '''
    The GT object defines and stores parameters of a gas turbine

    Parameters
    ----------
    model : str
        Gas turbine model.
        'LM2500': General Electric LM2500 aeroderivative
        'LM6000': General Electric LM6000 aeroderivative
    outMin : float
        Minimum power ouput in Watts.
    outMax : float
        Maximum power ouput in Watts.
    powerOut : float
        Current power ouput in Watts.
        If powerOut = 0, switches off GT. Otherwise, outMin <= powerOut <=outMax

    Attributes
    ----------
    powerRated : float
        Rated power of the turbogenerator in Watts. 
    inertia : float
        Inertia constant of the turbogenerator in seconds.
    rocP : float
        Maximum rate of change of power in Watts per second.
    
    '''

    def __init__(self, model, outMin, outMax, powerOut=0):
        self.__initComplete = False
        self.model = model
        self.outMin = outMin
        self.outMax = outMax

        self.__initComplete = True
        if not self.minmaxValid(self.outMin, self.outMax):
            raise ValueError("GT min or max output are invalid.") 
        
        self.powerOut = powerOut

    def minmaxValid(self, outMin, outMax):
        '''
        Check if the min and max outputs values are valid.

        Parameters
        ----------
        outMin : float
            Minimum power ouput in Watts.
            Must be non-negative, <= outMax and <= powerRated. 
        outMax : float
            Maximum power ouput in Watts.
            Must be non-negative, >= outMin and <= powerRated. 

        Returns
        ----------
        True or False

        '''
        return (outMin <= outMax) and (outMax <= self.powerRated)

    @property
    def model(self):
        return self.__model
    @model.setter
    def model(self, model):
        self.__model = model
        if model == 'LM2500':
            # Based on the datasheet values from Edvard Grieg
            #  Document 23380E-KVEST-251-E-DS-00002
            self.__powerRated = 30.8e6
            self.__inertia = 1.85
            self.__rocP = 0.15 * self.powerRated 
        elif model == 'LM6000':
            # Based on the document "LM6000-60 HZ Gas Turbine Generator Set
            #  Product Specification" and Cigré Report 238 "Modeling of Gas
            #  Turbines and Steam Turbines in Combined-Cycle Power Plants"
            self.__powerRated = 44.7e6
            self.__inertia = 1.8
            self.__rocP = 0.15 * self.powerRated 
        else:
            raise ValueError("GT model is invalid.")
    
    @property
    def outMin(self):
        return self.__outMin
    @outMin.setter
    def outMin(self, outMin):
        if outMin < 0:
            raise ValueError("GT min ouput cannot be negative.")
        elif self.__initComplete:
            if not self.minmaxValid(outMin, self.outMax):
                raise ValueError("GT min ouput is invalid.")
        self.__outMin = outMin
        self.powerOut = self.powerOut

    @property
    def outMax(self):
        return self.__outMax
    @outMax.setter
    def outMax(self, outMax):
        if outMax < 0:
            raise ValueError("GT max output cannot be negative.")
        elif self.__initComplete:
            if not self.minmaxValid(self.outMin, outMax):
                raise ValueError("GT max ouput is invalid.")
        self.__outMax = outMax
        self.powerOut = self.powerOut

    @property
    def powerRated(self):
        return self.__powerRated

    @property
    def inertia(self):
        return self.__inertia

    @property
    def rocP(self):
        return self.__rocP
    
    @property
    def powerOut(self):
        return self.__powerOut
    @powerOut.setter
    def powerOut(self, powerOut):
        if (powerOut == 0) or (self.outMin <= powerOut <= self.outMax):
            self.__powerOut = powerOut 
        else:
            self.__powerOut = max(self.outMin, min(powerOut, self.outMax))


class WT:
    '''
    The WT object defines and stores parameters of a wind turbine

    Parameters
    ----------
    model : str
        Wind turbine model.
        'NREL': NREL 5MW
        'Hywind': Hywind Tampen 6MW

    Attributes
    ----------
    powerRated : float
        Rated power of the wind turbine in Watts. 
    
    '''
    def __init__(self, model):
        self.__initComplete = False
        self.model = model

        self.__initComplete = True

    def powerOut(self, windSpeed):
        '''
        Given the wind speed at hub height, provides the output power.

        Parameters
        ----------
        windSpeed : float
            The wind speed at hub height in meters per second.

        Returns
        ----------
        powerOut : float
            The ouput power in Watts.

        '''
        
        return np.interp(windSpeed, self.__inWindSpeed, self.__outPower, left=0, right=0) * self.powerRated

    @property
    def model(self):
        return self.__model
    @model.setter
    def model(self, model):
        self.__model = model
        if model == 'NREL':
            # Based on report NREL/TP-500-38060
            self.__powerRated = 5e6
            self.__inWindSpeed = list(range(2,12)) + [25,26]
            self.__outPower = [0, 0.034, 0.0782, 0.1462, 0.2346, 0.3504, 0.5068, 0.6904, 0.9116, 1, 1, 0]
        elif model == 'Hywind':
            # Based on based on https://www.uib.no/sites/w3.uib.no/files/attachments/hywind_energy_lab.pdf (page 25)
            self.__powerRated = 6e6
            self.__inWindSpeed = list(range(3,16)) + [25,26]
            self.__outPower = [0, 0.029, 0.0725, 0.1304, 0.2101, 0.3261, 0.4638, 0.6232, 0.7754, 0.8913, 0.9565, 0.9855, 1, 1, 0]
        else:
            raise ValueError("WT model is invalid.")

    @property
    def powerRated(self):
        return self.__powerRated

        


def elgridEval(freqMax, freqMin, GT, WT, ES, LD):
    '''
    For a given configuration (defined by the array of objects GT, WT, ES, LD), checks the frequency stability of the electrical grid.  

    Parameters
    ----------
    freqMax: float
        Maximum frequency allowed in Hertz.
    freqMin: float
        Minimum frequency allowed in Hertz.        
    GT : GT
        Array of gas turbine objects.
    WT : WT
        Array of wind turbine objects.
    ES : ES
        Array of energy storage objects.
    LD: LD
        Array of energy storage objects.

    Returns
    ----------
    True or False

    '''
    
    ## Initialize local variables    
    # Base power in MVA
    SBASE = 100e6

    # GTs aggregated power in per unit of SBASE
    powerGTpu = 0
    # GTs aggregated inertia in per unit of SBASE
    inertiaGTpu = 0
    # GTs minimum rocP
    rocPGT = 1e20

    # WTs aggregated power in per unit of SBASE
    powerGTpu = 0

    ## Calculation of aggregated values     
    for _ in GT:
        powerGTpu += _.powerRated / SBASE
        inertiaGTpu += ( _.inertia * _.powerRated / SBASE)
        rocPGT = min(rocPGT, _.rocP)
        
    for _ in WT:
        powerWTpu += _.powerRated / SBASE
         

    return True
        
        
        
        

        

