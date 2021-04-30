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

# Import packages
import hes_off_object_oriented.elgrid_models as elgev


def test_GT_initialization():

    # Empty initialization
    testGT = elgev.GT(model='LM2500', powerMin=10e6, powerMax=30e6)
    print()
    print('{:<20} \t {}'.format("GT model:", testGT.model))
    print('{:<20} \t {}'.format("GT rated power:", testGT.powerRated))
    print('{:<20} \t {}'.format("GT min power:", testGT.powerMin))
    print('{:<20} \t {}'.format("GT max power:", testGT.powerMax))
    print('{:<20} \t {}'.format("GT moment of inertial:", testGT.inertia))
    print('{:<20} \t {}'.format("GT max power RoC:", testGT.rocP))
    print('{:<20} \t {}'.format("GT frequency damping:", testGT.damping))

def test_GT_set_damping():

    return


# # Test area for GT class
# print('## GT class ##')
# testGT = elgev.GT('LM2500', 10e6, 30e6)
# print(testGT.inertia)
# print(testGT.rocP)
# testGT.outMin = 20e6
# print(testGT.outMin)



test_GT_initialization()

# # Test area for WT class
# print('## WT class ##')
# testWT = elgev.WT('Hywind')
# print(testWT.windSpeed)
# print(testWT.powerOut())
# print(testWT.powerOut([3,4]))
# print(testWT.powerOut(12))
# print(testWT.windSpeed)
# print(testWT.powerOut())
#
#
# # Test area for ES class
# print('## ES class ##')
# testES = elgev.ES(6e6, 2.5e6, 1.3478e6)
# testES.powerOut = -3e6
# print(testES.powerIn)
#
# # Test area for LD class
# print('## LD class ##')
# testLD = elgev.LD(36.8e6, 7e6, 6.8e6)
# print(testLD.continuous)
# testLD.powerIn = 10e6
# print(testLD.powerIn)
#
# # Test area for elgridEval function
# print('## elgridEval function ##')
# myTestGTs = [('LM2500', 10e6, 30e6, 15e6), ('LM6000', 10e6, 30e6, 15e6)]
# myGTList = [elgev.GT(_model, _outMin, _outMax, _powerOut) for (_model, _outMin, _outMax, _powerOut) in myTestGTs]
#
# myTestWTs = [('Hywind', 12), ('NREL', 12)]
# myWTList = [elgev.WT(_model, _windSpeed) for (_model, _windSpeed) in myTestWTs]
#
# myTestESs = [(6e6, 0, 1.3478e6), (0, 2.5e6, 0)]
# myESList = [elgev.ES(_outMax, _outMin, _powerIn) for (_outMax, _outMin, _powerIn) in myTestESs]
#
# myTestLDs = [(36.8e6, 16.8e6, 7e6, 21e6)]
# myLDList = [elgev.LD(_continuous, _largest, _flexible, _powerIn) for (_continuous, _largest, _flexible, _powerIn) in myTestLDs]
#
# testelgridEval = elgev.elgridEval(myGTList, myWTList, myESList, myLDList)
#
# print(testelgridEval)
# print(testelgridEval.freqOk)