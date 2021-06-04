## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
##                         ___   ___  ________    ________         ______    ________ ________                        ##
##                        |  |  |  | |   ____|   /       |        /  __  \  |   ____||   ____|                        ##
##                        |  |__|  | |  |__     |   (----` ______|  |  |  | |  |__   |  |__                           ##
##                        |   __   | |   __|     \   \    |______|  |  |  | |   __|  |   __|                          ##
##                        |  |  |  | |  |____.----)   |          |  `--'  | |  |     |  |                             ##
##                        |__|  |__| |_______|_______/            \______/  |__|     |__|                             ##
##                                                                                                                    ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##

# Import external packages
import numpy as np
from collections import namedtuple


## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------------------------------------------ ##


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
    powerMin : float
        Minimum power ouput in Watts.
    powerMax : float
        Maximum power ouput in Watts.

    Attributes
    ----------
    powerRated : float
        Rated power of the turbogenerator in Watts. 
    inertia : float
        Inertia constant of the turbogenerator in seconds.
    rocP : float
        Maximum rate of change of power in per unit of powerRated per second.
    damping : float
        Current damping / P(f) droop in per unit of powerRated per per unit of rated frequency
    powerOut : float
        Current power ouput in Watts.
        If = 0, switches off GT. Otherwise, outMin <= powerOut <=outMax
    
    '''

    def __init__(self, model, powerMin=None, powerMax=None, powerOut=0):
        self.__initComplete = False
        self.__powerOut = 0
        self.model = model
        self.powerMin = powerMin
        self.powerMax = powerMax

        self.__initComplete = True
        if not self.minmaxValid(self.powerMin, self.powerMax):
            raise ValueError("GT min or max output are invalid.") 
        
        self.powerOut = powerOut
        self.set_damping(0, 0)

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
    def powerMin(self):
        return self.__outMin
    @powerMin.setter
    def powerMin(self, outMin):
        if outMin < 0:
            raise ValueError("GT min ouput cannot be negative.")
        elif self.__initComplete:
            if not self.minmaxValid(outMin, self.powerMax):
                raise ValueError("GT min ouput is invalid.")
        self.__outMin = outMin
        self.powerOut = self.powerOut

    @property
    def powerMax(self):
        return self.__outMax
    @powerMax.setter
    def powerMax(self, outMax):
        if outMax < 0:
            raise ValueError("GT max output cannot be negative.")
        elif self.__initComplete:
            if not self.minmaxValid(self.powerMin, outMax):
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
    def damping(self):
        return self.__damping

    def set_damping(self, freqMin, freqMax):
        '''
        Set damping based on the current values for outMin and outMax.
        If freqMin = freqMax, then damping = 0.

        Parameters
        ----------
        freqMin : float
            Minimum frequency deviation in per unit.
            Must be <= 0.
        freqMax : float
            Maximum frequency deviation in per unit.
            Must be >= 0.

        Returns
        ----------
        damping : float

        '''

        if (freqMin > 0) or (freqMax < 0):
            raise ValueError("GT.setdamping() invalid use.")

        if (freqMin == freqMax):
            self.__damping = 0
        else:
            self.__damping = (self.powerMax - self.powerMin) / self.powerRated / (freqMax - freqMin)

        return self.damping


    @property
    def powerOut(self):
        return self.__powerOut
    @powerOut.setter
    def powerOut(self, powerOut):
        if (powerOut == 0) or (self.powerMin <= powerOut <= self.powerMax):
            self.__powerOut = powerOut 
        else:
            self.__powerOut = max(self.powerMin, min(powerOut, self.powerMax))


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
    damping : float
        Current damping / P(f) droop in per unit of flexible per per unit of rated frequency
        
    '''
    def __init__(self, continuous, largest, flexible, powerIn=0.0):
        self.__initComplete = False
        self.continuous = continuous
        self.largest = largest
        self.flexible = flexible
        self.powerIn = powerIn

        self.__initComplete = True

    def setdamping(self, freqMin, freqMax):
        '''
        Set damping based on the current values for flexible.
        If freqMin = freqMax, then damping = 0. 

        Parameters
        ----------
        freqMin : float
            Minimum frequency deviation in per unit.
            Must be <= 0. 
        freqMax : float
            Maximum frequency deviation in per unit.
            Must be >= 0. 

        Returns
        ----------
        damping : float

        '''

        if (freqMin > 0) or (freqMax < 0):
            raise ValueError("LD.setdamping() invalid use.")

        if (freqMin == freqMax):
            self.__damping = 0 
        else:
            self.__damping = 1 / (freqMax - freqMin)
            
        return self.damping


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

    @property
    def damping(self):
        return self.__damping


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

        


def elgridEval(GT, WT, ES, LD, freqSSMax=0.02, freqSSMin=-0.02, freqTrMax=0.05, freqTrMin=-0.05):
    '''
    For a given configuration (defined by the array of objects GT, WT, ES, LD), evaluates the frequency stability of the electrical grid.  

    Parameters
    ----------
    freqSSMax: float
        Maximum allowed steady-state frequency deviation in per unit.
    freqSSMin: float
        Minimum allowed steady-state frequency deviation in per unit.        
    freqTrMax: float
        Maximum allowed transient frequency deviation in per unit.
    freqTrMin: float
        Minimum allowed transient frequency deviation in per unit.        
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
    freqOk: bool
        if True, the grid is frequency stable (freqSSMin <= frequency <= freqSSMax).
    freqMinOk: bool
        if True, the grid is frequency stable for the largest load increase or generation reduction (freqSSMin <= frequency).
    freqMinCalc: float
        Minimum calculated steady-state frequency deviation in per unit.
    freqMaxOk: bool
        if True, the grid is frequency stable for the largest load reduction or generation increse (frequency <= freqSSMax).
    freqMaxCalc: float
        Maximum calculated steady-state frequency deviation in per unit.
    rocPOk: bool (not implemented)
        if True, the largest rate of change of active power (rocP) of all gas turbines is below their allowed limits.
    rocPGT: float (not implemented)
        list with the calculated rocP of each GT object in pu per second. 
   
 
    '''
    
    # Base power in MVA for per unit (pu) / normalization
    SBASE = 100e6

    # Calculation of GTs aggregated values     
    powerGTpu = 0
    inertiaGTpu = 0
    dampGTpu = 0
    for _ in GT:
        powerGTpu += _.powerOut
        inertiaGTpu += _.inertia * _.powerRated
        dampGTpu += _.set_damping(freqSSMin, freqSSMax) * _.powerRated
    powerGTpu /= SBASE
    inertiaGTpu /= SBASE
    dampGTpu /= SBASE
        
    # Calculation of WTs aggregated values     
    powerWTpu = 0
    PbNegWT = 0
    for _ in WT:
        powerWTpu += _.powerOut()
        PbNegWT = max(PbNegWT, _.powerOut())
    powerWTpu /= SBASE
    PbNegWT /= SBASE 
             
    # Calculation of ESs aggregated values     
    powerESpu = 0
    PbNegES = 0
    PbPosES = 0
    for _ in ES:
        powerESpu += _.powerIn
        if _.powerIn < 0:
            PbNegES = max(PbNegES, -_.powerIn)
        else:
            PbPosES = max(PbPosES, _.powerIn)
    powerESpu /= SBASE
    PbNegES /= SBASE 
    PbPosES /= SBASE
    
    # Calculation of LDs aggregated values     
    powerLDpu = 0
    PbNegLD = 0
    PbPosLD = 0
    dampLDpu = 0
    for _ in LD:
        powerLDpu += _.powerIn
        if _.powerIn > (_.continuous - _.largest):
            PbPosLD = max(PbPosLD, _.largest)
            PbNegLD = max(PbNegLD, _.continuous - _.powerIn)
        else:
            PbPosLD = max(PbPosLD, _.powerIn)
            PbNegLD = max(PbNegLD, _.largest)
        dampLDpu += _.set_damping(freqSSMin, freqSSMax) * _.flexible

    powerLDpu /= SBASE
    PbNegLD /= SBASE 
    PbPosLD /= SBASE
    dampLDpu /= SBASE

    # Check worst negative power imbalance (load > generation)
    PbNeg = max(PbNegWT, PbNegES, PbNegLD)
    PbNegMax = freqEval(freqSS=-freqSSMin, freqTr=-freqTrMin, Dmin=(dampGTpu + dampLDpu))
    freqMinOk = PbNeg < PbNegMax
    freqMinCalc = -freqEval(freqTr=-freqTrMin, Dmin=(dampGTpu + dampLDpu), Pb=PbNeg)

    # Check worst positive power imbalance (load < generation)
    PbPos = max(PbPosES, PbPosLD)
    PbPosMax = freqEval(freqSS=freqSSMax, freqTr=freqTrMax, Dmin=dampGTpu)
    freqMaxOk = PbPos < PbPosMax
    freqMaxCalc = freqEval(freqTr=freqTrMax, Dmin=dampGTpu, Pb=PbPos)

    # TO-DO: Check if GTs rocP are within the allowed limits
    rocFreq = (dampGTpu + dampLDpu) / inertiaGTpu
    rocPOk = True
    rocPMax = 0
    rocPGT = []
    '''
    for _ in GT:
        rocPGT.append(_.damping * rocFreq) 
        rocPOk = rocPOk and (rocPGT[-1] < _.rocP)
        rocPMax = max(rocPMax, rocPGT[-1])
    '''
    
    # Return results
    Desc = namedtuple("Desc", ["freqOk", "freqMinOk", "freqMinCalc", "freqMaxOk", "freqMaxCalc", "rocPOk", "rocPGT"])
    
    return Desc(
        freqMinOk and freqMaxOk,
        freqMinOk,
        freqMinCalc,
        freqMaxOk,
        freqMaxCalc,
        rocPOk,
        rocPGT
    )
        
        
def freqEval(freqSS=0, freqTr=0, Dmin=0, Pb=0):
    '''
    Evaluates the frequency stability of an electrical grid as described in https://doi.org/10.1109/TPWRS.2020.3039832    

    Parameters
    ----------
    freqSS: float
        Maximum steady-state frequency deviation in per unit.
    freqTr: float
        Maximum transient frequency deviation in per unit.
    Dmin: float
        Minimum system damping in per unit.
    Dmin: float
        Maximum active power disturbance in per unit.

    Returns
    ----------
    The parameter not defined.
    Raises error if all parameters are defined or more than one is not defined. 

    '''
    
    if (freqSS, freqTr, Dmin, Pb).count(0) > 1:
        raise ValueError("freqEval() invalid use.") 
    
    if freqSS == 0:
        return Pb / (Dmin * (1 - freqTr))
    elif freqTr == 0:
        return 1 - (Pb / (Dmin * freqSS) )
    elif Dmin == 0:
        return Pb / (freqSS * (1 - freqTr))
    elif Pb == 0:
        return Dmin * freqSS * (1 - freqTr)
    else:
        raise ValueError("freqEval() invalid use.")
    
    
    
