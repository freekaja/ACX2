import numpy
import os
from sherpa.models import model
# Default spot is at {HOME}/atomdb
# export ATOMDB=/Users/andy/atomdb
# print( str(os.path.expanduser("~")) + "/atomdb")
os.environ["ATOMDB"] = str(os.path.expanduser("~"))  + "/atomdb"

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

vvacx2Info = ("temperature   \"keV\"   1.0 0.00862 0.00862 86. 86. 0.01",
              "collnpar    \"kev/u,km/s\" 1.0 0.01 0.2 100. 1000. 0.01",
              "collntype     \"\"      1.0 1.0 1.0 4.0 4.0 -0.01",
              "acxmodel      \"\"      8.0 1.0 1.0 16.0 16.0 -0.01",
              "recombtype    \"\"      1.0 1.0 1.0 2.0 2.0 -0.01",
              "Hefrac        \"\"      0.09 0.0 0.0 1.0 1.0 -0.01",
              "H             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "He            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Li            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Be            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "B             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "C             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "N             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "O             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "F             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ne            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Na            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Mg            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Al            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Si            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "P             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "S             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Cl            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ar            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "K             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ca            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Sc            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ti            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "V             \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Cr            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Mn            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Fe            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01",
              "Ni            \"\"      1.0 0.0 0.0 10.0 10.0 -0.01")



def _vvacx2(params,engs ,flux):

  """
  VVACX2 model for data

  PARAMETERS
  ----------
  engs : list[float]
    The energy bin edges (from xspec)
  params : list[float]
    The parameter list. See vvacx2Info for definition
  flux : list[float]
    The array to fill with return values

  RETURNS
  -------
  None
    Fills out the flux array with photon cm3 s-1 bin-1, x1e10

  USAGE
  -----
    # load the model into XSPEC
    xspec.AllModels.addPyMod(vvacx2, vvacx2Info, 'add')
    # make a model
    m = xspec.Model('vvacx2')
  """

  # This is the call that will return everything. So set everything!
  ebins = numpy.array(engs)

  # vacx model has the all elements up to nickel except cobalt
  elements =[ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,\
             11, 12, 13, 14, 15, 16, 17, 18, 19, 20,\
             21, 22, 23, 24, 25, 26, 28]

  if len(acx2_acxmodelobject.DonorList)==0:
    # Add a donor - hydrogen in this case (default will have H, He)
    acx2_acxmodelobject.add_donor('H', \
                                  Hlinefile, \
                                  Hcontfile, \
                                  Hsigmafile,\
                        elements = elements)
    acx2_acxmodelobject.add_donor('He', \
                                  Helinefile, \
                                  Hecontfile, \
                                  Hesigmafile,\
                        elements = elements)


  # check energies
  acx2_acxmodelobject.set_ebins(ebins)

  # check temperature
  acx2_acxmodelobject.set_temperature(params[0])

  acx2_acxmodelobject.set_acxmodel(int(params[3]))

  # set recombination type (1 = single ['solar wind'], 2= full ['comet'])
  acx2_acxmodelobject.set_recombtype(int(params[4]))

  # set abundance vector
  abund = numpy.array(params[6:33])
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
  flux= spec*1e10
  return flux




class VVACX2(model.RegriddableModel1D):
    def __init__(self, name ='VVACX2'):
        self.temperature = model.Parameter(name,'temperature',units='keV',val=1.0,hard_min=0.00862,min=0.00862,max=86,hard_max=86)
        self.collnpar =  model.Parameter(name,'collnpar',units='kev/u,km/s',val=1.0,hard_min=0.01,min=0.2,max=100,hard_max=1000)
        self.collntype =  model.Parameter(name,'collntype',val=1.0,hard_min=1.0,min=1.0,max=4.0,hard_max=4.0)
        self.acxmodel =  model.Parameter(name,'acxmodel',val=8.0,hard_min=1.0,min=1.0,max=16.0,hard_max=16.0)
        self.recombtype =  model.Parameter(name,'recombtype',val=1.0,hard_min=1.0,min=1.0,max=2.0,hard_max=2.0)
        self.Hefrac =  model.Parameter(name,'Hefrac',val=0.09,hard_min=0.0,min=0.0,max=1.0,hard_max=1.0)
        self.H =  model.Parameter(name,'H',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.He =  model.Parameter(name,'He',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Li =  model.Parameter(name,'Li',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Be =  model.Parameter(name,'Be',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.B =  model.Parameter(name,'B',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.C =  model.Parameter(name,'C',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.N =  model.Parameter(name,'N',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.O =  model.Parameter(name,'O',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.F =  model.Parameter(name,'F',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Ne =  model.Parameter(name,'Ne',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Na =  model.Parameter(name,'Na',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Mg =  model.Parameter(name,'Mg',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Al =  model.Parameter(name,'Al',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Si =  model.Parameter(name,'Si',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.P =  model.Parameter(name,'P',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.S =  model.Parameter(name,'S',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Cl =  model.Parameter(name,'Cl',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Ar =  model.Parameter(name,'Ar',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.K =  model.Parameter(name,'K',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Ca =  model.Parameter(name,'Ca',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Sc =  model.Parameter(name,'Sc',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Ti =  model.Parameter(name,'Ti',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.V =  model.Parameter(name,'V',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Cr =  model.Parameter(name,'Cr',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Mn =  model.Parameter(name,'Mn',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Fe =  model.Parameter(name,'Fe',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)
        self.Ni =  model.Parameter(name,'Ni',val=1.0,hard_min=0.0,min=0.0,max=10.0,hard_max=10.0)

        model.RegriddableModel1D.__init__(self, name,
                                          (self.temperature, self.collnpar, self.collntype,
                                           self.acxmodel,self.recombtype,self.Hefrac,self.H,
                                           self.He,self.Li,self.Be,self.B,self.C,self.N,self.O,
                                           self.F,self.Ne,self.Na,self.Mg,self.Al,self.Si,self.P,
                                           self.Si,self.P,self.S,self.Ar,self.K,self.Ca,self.Sc,
                                           self.Ti,self.V,self.Cr,self.Mn,self.Fe,self.Ni))

    def calc(self,pars,engs,flux,*args,**kwargs):
        """
        Evaluate the model
        """
        return _vvacx2(pars,engs=engs,flux=flux)

