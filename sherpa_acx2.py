import numpy
import os
from sherpa.models import model
# Default spot is at {HOME}/atomdb
# export ATOMDB=/Users/andy/atomdb
# print( str(os.path.expanduser("~")) + "/atomdb")
# os.environ["ATOMDB"] = str(os.path.expanduser("~"))  + "/atomdb"

import acx2 as acx2model

# CHANGE THESE FILE PATHS TO REFLECT YOUR SYSTEM

Hsigmafile  = 'acx2_H_v1_sigma.fits'
Hlinefile   = 'acx2_H_v1_line.fits'
Hcontfile   = 'acx2_H_v1_cont.fits'
Hesigmafile = 'acx2_He_v1_sigma.fits'
Helinefile  = 'acx2_He_v1_line.fits'
Hecontfile  = 'acx2_He_v1_cont.fits'


acx2_acxmodelobject = acx2model.ACXModel()

# para_name, units, default_val, hard_min, soft_min, soft_max, hard_max, fit_delta
# 
acx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
            "collnpar   \"kev/u,km/s\" 1.0 0.01 0.2 100. 1000. 0.01",
            "collntype     \"\"      1.0 1.0 1.0 4.0 4.0 -0.01",
            "acxmodel      \"\"      8.0 1.0 1.0 16.0 16.0 -0.01",
            "recombtype    \"\"      1.0 1.0 1.0 2.0 2.0 -0.01",
            "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
            "abund         \"\"      1.0 0.0 0.0 10.0 10.0 0.01")


def _acx2(params, engs, flux):

  """
  ACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See acx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(acx2, acx2Info, 'add')
    # make a model
    m = xspec.Model('acx2')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)

  # vacx model has the 14 main elements
  elements = [1,2,6,7,8,10,12,13,14,16,18,20,26,28]

  if len(acx2_acxmodelobject.DonorList)==0:
    # Add a donor - hydrogen in this case (default will have H, He)
    acx2_acxmodelobject.add_donor('H', \
                                  Hlinefile, \
                                  Hcontfile, \
                                  Hsigmafile,\
                                  elements = numpy.array(elements))

    acx2_acxmodelobject.add_donor('He', \
                                  Helinefile, \
                                  Hecontfile, \
                                  Hesigmafile,\
                                  elements = numpy.array(elements))

  # check energies
  acx2_acxmodelobject.set_ebins(ebins)

  # check temperature
  acx2_acxmodelobject.set_temperature(params[0])

  acx2_acxmodelobject.set_acxmodel(int(params[3]))


  # set recombination type (1 = single ['solar wind'], 2= full ['comet'])
  acx2_acxmodelobject.set_recombtype(int(params[4]))


  # set abundance vector
  abund = numpy.array(params[6])
  acx2_acxmodelobject.set_abund(abund, elements=elements)

  # set the H & He fraction from HeFrac
  acx2_acxmodelobject.set_donorabund(['H','He'], [1-params[5], params[5]])

  # set the collision type (1,2,3 or 4)
  cp = int(params[2])
  if cp== 1:
    cpunits = 'kev/amu'
  elif cp in [2,3,4]:
    cpunits = 'km/s'
  else:
    raise(ValueError, "Unknown value %i for collntype (should be 1,2,3 or 4)"%(cp))

  acx2_acxmodelobject.set_collisiontype(cp, cpunits)


  # get the spectrum
  spec = acx2_acxmodelobject.calc_spectrum(params[1])

  # return the flux.
  flux = spec*1e10
  return flux

class ACX2(model.RegriddableModel1D):
    def __init__(self, name ='ACX2'):
        self.temperature = model.Parameter(name,'temperature',units='keV',val=1.0,hard_min=0.00862,min=0.00862,max=86,hard_max=86)
        self.collnpar =  model.Parameter(name,'collnpar',units='kev/u,km/s',val=1.0,hard_min=0.01,min=0.2,max=100,hard_max=1000)
        self.collntype =  model.Parameter(name,'collntype',val=1.0,hard_min=1.0,min=1.0,max=4.0,hard_max=4.0)
        self.acxmodel =  model.Parameter(name,'acxmodel',val=8.0,hard_min=1.0,min=1.0,max=16.0,hard_max=16.0)
        self.recombtype =  model.Parameter(name,'recombtype',val=1.0,hard_min=1.0,min=1.0,max=2.0,hard_max=2.0)
        self.Hefrac =  model.Parameter(name,'Hefrac',val=0.09,hard_min=0.0,min=0.0,max=1.0,hard_max=1.0)
        self.abund =  model.Parameter(name,'abund',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)

        model.RegriddableModel1D.__init__(self, name,
                                          (self.temperature, self.collnpar, self.collntype,
                                           self.acxmodel,self.recombtype,self.Hefrac,self.abund))

    def calc(self,pars,engs,flux,*args,**kwargs):
        """
        Evaluate the model
        """
        return _acx2(pars,engs=engs,flux=flux)

