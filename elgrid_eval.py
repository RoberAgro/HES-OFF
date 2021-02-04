''' Module with classes and methods to evaluate the electrical grid '''

import numpy as np


class ES:
    '''
    The ES object defines and stores parameters of an energy storage system

    Parameters
    ----------
    outMax : float
        Rated power for charging in Watts.
    outMin : float
        Rated power for discharging in Watts.        
    powerIn : float
        Current power input in Watts.
        If = 0, switches off ES.
        If < 0, ES is discharging.
        If > 0, ES is charging.
    
    '''
    def __init__(self, outMax, outMin, powerIn=0):
        self.__initComplete = False
        self.outMax = outMax
        self.outMin = outMin
        self.powerIn = powerIn

        self.__initComplete = True

    @property
    def outMax(self):
        return self.__outMax
    @outMax.setter
    def outMax(self, outMax):
        if outMax < 0:
            raise ValueError("ES outMax must be positive.")
        self.__outMax = outMax

    @property
    def outMin(self):
        return self.__outMin
    @outMin.setter
    def outMin(self, outMin):
        if outMin < 0:
            raise ValueError("ES outMin must be positive.")
        self.__outMin = outMin

    @property
    def powerIn(self):
        return self.__powerIn
    @powerIn.setter
    def powerIn(self, powerIn):
        if (powerIn == 0) or (-self.outMin <= powerIn <= self.outMax):
            self.__powerIn = powerIn 
        else:
            self.__powerIn = max(-self.outMin, min(powerIn, self.outMax))


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
        If = 0, switches off GT. Otherwise, outMin <= powerOut <=outMax

    Attributes
    ----------
    powerRated : float
        Rated power of the turbogenerator in Watts. 
    inertia : float
        Inertia constant of the turbogenerator in seconds.
    rocP : float
        Maximum rate of change of power in pu of powerRated per second.
    
    '''

    def __init__(self, model, outMin, outMax, powerOut=0):
        self.__initComplete = False
        self.__powerOut = 0
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
            self.__rocP = 0.15 
        elif model == 'LM6000':
            # Based on the document "LM6000-60 HZ Gas Turbine Generator Set
            #  Product Specification" and Cigré Report 238 "Modeling of Gas
            #  Turbines and Steam Turbines in Combined-Cycle Power Plants"
            self.__powerRated = 44.7e6
            self.__inertia = 1.8
            self.__rocP = 0.15 
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


class LD:
    '''
    The LD object defines and stores parameters of electric loads

    Parameters
    ----------
    continuous : float
        Installed power of continuous loads in Watts.
    largest : float
        Installed power of the largest continuous load in Watts.        
    flexible : float
        Installed power of flexible loads in Watts.        
    powerIn : float
        Current active power demand in Watts.                    
    '''
    def __init__(self, continuous, largest, flexible, powerIn=0.0):
        self.__initComplete = False
        self.continuous = continuous
        self.largest = largest
        self.flexible = flexible
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
    def largest(self):
        return self.__largest
    @largest.setter
    def largest(self, largest):
        if largest < 0:
            raise ValueError("LD largest cannot be negative.")
        elif largest > self.continuous:
            raise ValueError("LD largest is invalid.")
        self.__largest = largest

    @property
    def flexible(self):
        return self.__flexible
    @flexible.setter
    def flexible(self, flexible):
        if flexible < 0:
            raise ValueError("LD flexible cannot be negative.")
        self.__flexible = flexible

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


class WT:
    '''
    The WT object defines and stores parameters of a wind turbine

    Parameters
    ----------
    model : str
        Wind turbine model.
        'NREL': NREL 5MW
        'Hywind': Hywind Tampen 6MW
    windSpeed : float
        Current wind speed at hub height in meters per second. 

    Attributes
    ----------
    powerRated : float
        Rated power of the wind turbine in Watts. 
    
    '''
    def __init__(self, model, windSpeed = 0):
        self.__initComplete = False
        self.model = model
        self.windSpeed = windSpeed

        self.__initComplete = True

    def powerOut(self, *args):
        '''
        Given the wind speed at hub height, provides the output power.
        If only one input is given, stores it for future use.
        If no argument is given, outputs the power of the stored wind speed.
        If several arguments are given, return a list of output powers.   

        Parameters
        ----------
        *args : list of float
            The wind speed at hub height in meters per second.

        Returns
        ----------
        powerOut : list of float
            The ouput power in Watts.

        '''

        powerOut = []
        for _ in args:
            powerOut.append(np.interp( _, self.__pcurveWindSpeed, self.__pcurvePower, left=0, right=0) * self.powerRated)         

        if len(powerOut) == 0:
            powerOut = np.interp(self.windSpeed, self.__pcurveWindSpeed, self.__pcurvePower, left=0, right=0) * self.powerRated
        elif len(powerOut) == 1:
            powerOut = powerOut[0]
            self.__windSpeed = args[0]

        return powerOut

    @property
    def model(self):
        return self.__model
    @model.setter
    def model(self, model):
        self.__model = model
        if model == 'NREL':
            # Based on report NREL/TP-500-38060
            self.__powerRated = 5e6
            self.__pcurveWindSpeed = list(range(2,12)) + [25,26]
            self.__pcurvePower = [0, 0.034, 0.0782, 0.1462, 0.2346, 0.3504, 0.5068, 0.6904, 0.9116, 1, 1, 0]
        elif model == 'Hywind':
            # Based on based on https://www.uib.no/sites/w3.uib.no/files/attachments/hywind_energy_lab.pdf (page 25)
            self.__powerRated = 6e6
            self.__pcurveWindSpeed = list(range(3,16)) + [25,26]
            self.__pcurvePower = [0, 0.029, 0.0725, 0.1304, 0.2101, 0.3261, 0.4638, 0.6232, 0.7754, 0.8913, 0.9565, 0.9855, 1, 1, 0]
        else:
            raise ValueError("WT model is invalid.")

    @property
    def windSpeed(self):
        return self.__windSpeed
    @windSpeed.setter
    def windSpeed(self, windSpeed):
        if 0 <= windSpeed <= 50:
            self.__windSpeed = windSpeed 
        else:
            self.__windSpeed = 0

    @property
    def powerRated(self):
        return self.__powerRated

        


def elgridEval(freqMax, freqMin, GT, WT, ES, LD):
    '''
    For a given configuration (defined by the array of objects GT, WT, ES, LD), checks the frequency stability of the electrical grid.  

    Parameters
    ----------
    freqMax: float
        Maximum frequency allowed in per unit.
    freqMin: float
        Minimum frequency allowed in per unit.        
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
    
    # Base power in MVA for per unit (pu) / normalization
    SBASE = 100e6

    # Calculation of GTs aggregated values     
    powerGTpu = 0
    inertiaGTpu = 0
    rocPGTpu = 1e20
    for _ in GT:
        powerGTpu += _.powerOut
        inertiaGTpu += _.inertia
        rocPGTpu = min(rocPGTpu, _.rocP * _.powerRated)
    powerGTpu /= SBASE
    inertiaGTpu *= _.powerRated / SBASE
    rocPGTpu /= SBASE
        
    # Calculation of WTs aggregated values     
    powerWTpu = 0
    maxpowerOutStepWT = 0
    for _ in WT:
        powerWTpu += _.powerOut()
        maxpowerOutStepWT = max(maxpowerOutStepWT, _.powerOut())
    powerWTpu /= SBASE
    maxpowerOutStepWT /= SBASE 
             
    # Calculation of ESs aggregated values     
    powerESpu = 0
    maxpowerOutStepES = 0
    maxpowerInStepES = 0
    for _ in ES:
        powerESpu += _.powerIn
        if _.powerIn < 0:
            maxpowerOutStepES = max(maxpowerOutStepES, -_.powerIn)
        else:
            maxpowerInStepES = max(maxpowerInStepES, _.powerIn)
    powerESpu /= SBASE
    maxpowerOutStepES /= SBASE 
    maxpowerInStepES /= SBASE
    
    # Calculation of ESs aggregated values     
    powerLDpu = 0
    for _ in LD:
        powerLDpu += _.powerIn
    powerLDpu /= SBASE

    # Calculate the worst possible load increase   
        
        
        
        

        

