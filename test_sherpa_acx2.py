

import numpy as np
import sherpa

from sherpa_acx2 import ACX2


engs = np.linspace(0,1,11)
# print(engs)
flux = 100*engs


obj = ACX2()

print(obj)

pars = [1.0,1.0,1.0,8.0,1.0,0.09,1.0]

obj.calc(pars,engs=engs,flux=flux)