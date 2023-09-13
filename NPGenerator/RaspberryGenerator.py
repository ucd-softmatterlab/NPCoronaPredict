import numpy as np

def buildRaspberry( innerRadius,outerRadius, beadRadius):
    beadList = []
    gridPoints =2* (np.round(3.0*outerRadius/(2*np.sqrt(6)*beadRadius))+2).astype(int)
    for i in range(gridPoints):
        for j in range(gridPoints):
            for k in range(gridPoints):
                beadList.append([  beadRadius*(  2*i + ( (j+k)%2   )    ), 
                                 beadRadius*( np.sqrt(3) * (j + 1.0/3.0*(k%2))   ), 
                                 beadRadius*( 2*np.sqrt(6)/3.0 * k)   ])
    beadArray = np.array(beadList)
    beadArray = beadArray - outerRadius
    return beadArray[
        np.logical_and(
            beadArray[:,0]**2 + beadArray[:,1]**2 + beadArray[:,2]**2 < outerRadius**2,
            beadArray[:,0]**2 + beadArray[:,1]**2 + beadArray[:,2]**2 > innerRadius**2 
            )
        
        ]

