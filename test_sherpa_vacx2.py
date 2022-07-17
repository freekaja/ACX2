

import numpy as np
import sherpa

from sherpa_vacx2 import VACX2


engs = np.linspace(0,1,11)
# print(engs)
flux = 100*engs


obj = VACX2()

print(obj)

pars = [1.0,1.0,1.0,8.0,1.0,0.09,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]

print(len(pars))
obj.calc(pars,engs=engs,flux=flux)