import elgrid_eval as elgev

testGT = elgev.GT('LM2500', 10e6, 30e6)
print(testGT.inertia)

testGT.outMin = 20e6

print(testGT.outMin)

testWT = elgev.WT('Hywind')
print(testWT.powerCurve(35))

testES = elgev.ES(6e6, 2.5e6)
print(testES.powerCharge)
print(testES.powerDisch)


testLD = elgev.LD(6e6, 2.5e6)
print(testES.powerCharge)
print(testES.powerDisch)

