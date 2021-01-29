''' Module with classes and methods to evaluate the electrical grid '''

import numpy as np


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

    Attributes
    ----------
    powerRated : float
        Rated power of the turbogenerator in Watts. 
    inertia : float
        Inertia constant of the turbogenerator in seconds.
    
    '''

    def __init__(self, model, outMin, outMax):
        self.__initComplete = False
        self.model = model
        self.outMin = outMin
        self.outMax = outMax

        self.__initComplete = True
        if not self.powerValidate(self.outMin, self.outMax):
            raise ValueError("GT min or max output are invalid.") 


    def powerValidate(self, outMin, outMax):
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
        elif model == 'LM6000':
            # Based on the document "LM6000-60 HZ Gas Turbine Generator Set
            #  Product Specification" and Cigré Report 238 "Modeling of Gas
            #  Turbines and Steam Turbines in Combined-Cycle Power Plants"
            self.__powerRated = 44.7e6
            self.__inertia = 1.8
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
            if not self.powerValidate(outMin, self.outMax):
                raise ValueError("GT min ouput is invalid.")
        self.__outMin = outMin


    @property
    def outMax(self):
        return self.__outMax

    @outMax.setter
    def outMax(self, outMax):
        if outMax < 0:
            raise ValueError("GT max output cannot be negative.")
        elif self.__initComplete:
            if not self.powerValidate(self.outMin, outMax):
                raise ValueError("GT max ouput is invalid.")
        self.__outMax = outMax


    @property
    def powerRated(self):
        return self.__powerRated


    @property
    def inertia(self):
        return self.__inertia


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
    inertia : float
        Inertia constant of the turbogenerator in seconds.
    
    '''
    def __init__(self, model):
        self.__initComplete = False
        self.model = model

        self.__initComplete = True

    def powerCurve(self, windSpeed):
        '''
        Given the wind speed, provides the output power.

        Parameters
        ----------
        windSpeed : float
            The wind speed in meters per second.

        Returns
        ----------
        powerOut : float
            The ouput power in Watts.

        '''
        
        return np.interp(windSpeed, self.__windSpeed, self.__outPower, left=0, right=0) * self.__powerRated

    @property
    def model(self):
        return self.__model

    @model.setter
    def model(self, model):
        self.__model = model
        if model == 'NREL':
            # Based on report NREL/TP-500-38060
            self.__powerRated = 5e6
            self.__windSpeed = list(range(2,12)) + [25,26]
            self.__outPower = [0, 0.034, 0.0782, 0.1462, 0.2346, 0.3504, 0.5068, 0.6904, 0.9116, 1, 1, 0]
        elif model == 'Hywind':
            # Based on based on https://www.uib.no/sites/w3.uib.no/files/attachments/hywind_energy_lab.pdf (page 25)
            self.__powerRated = 6e6
            self.__windSpeed = list(range(3,16)) + [25,26]
            self.__outPower = [0, 0.029, 0.0725, 0.1304, 0.2101, 0.3261, 0.4638, 0.6232, 0.7754, 0.8913, 0.9565, 0.9855, 1, 1, 0]
        else:
            raise ValueError("WT model is invalid.")


testGT = GT('LM2500', 10e6, 30e6)
print(testGT.inertia)

testGT.outMin = 20e6

print(testGT.outMin)

testWT = WT('Hywind')
print(testWT.powerCurve(35))
