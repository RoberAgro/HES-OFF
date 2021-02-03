import elgrid_eval as elgev

# Test area for GT class
testGT = elgev.GT('LM2500', 10e6, 30e6)
print(testGT.inertia)
print(testGT.rocP)

testGT.outMin = 20e6

print(testGT.outMin)

# Test area for WT class
testWT = elgev.WT('Hywind')
print(testWT.powerOut(12))

# Test area for ES class
testES = elgev.ES(6e6, 2.5e6, 6e6/300, 2.5e6/300)
print(testES.powerCharge)
print(testES.rocPDisch)

# Test area for LD class
testLD = elgev.LD(36.8e6, 6.8e6, 7e6)
print(testLD.continuous)
testLD.powerIn = 10e6
print(testLD.powerIn)

# Test area for elgridEval method
myTestGTs = [('LM2500', 10e6, 30e6), ('LM2500', 10e6, 30e6)]
myGTList = [elgev.GT(_model,_outMin,_outMax) for (_model,_outMin,_outMax) in myTestGTs]

myTestLDs = [(36.8e6, 6.8e6, 7e6)]
myLDList = [elgev.LD(_continuous,_flexible,_largest) for (_continuous,_flexible,_largest) in myTestLDs]