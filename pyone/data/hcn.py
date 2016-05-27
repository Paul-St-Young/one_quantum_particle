#! /usr/bin/env python

import numpy as np
# ======= setup machine and jobs =======
from nexus import settings,Job,obj
machine_settings = obj(
    runs          = 'grid', # default to 'runs'
    generate_only = 0,
    status_only   = 0,
    sleep         = 1,
    machine       = 'ws1',
    ericfmt       = '/opt/intel14/gamess-2014R1/auxdata/ericfmt.dat'
)
settings(**machine_settings)
gms_job = Job(serial=True,app="rungms",cores=1)

# ======= structure =======
from nexus import generate_physical_system, Structure
from nexus import generate_gamess

CH = 1.0655 # ang
CN = 1.1538 # ang
basis = 'ccd'

# ======= simulations =======
'''
MODERN BASIS SET FAMILIES (CCN, PCN, MCP_NZP, IMCP) ARE INTENDED
 FOR USE ONLY AS SPHERICAL HARMONIC BASIS SETS. PLEASE SET ISPHER=1
  IN THE $CONTRL GROUP IN ORDER TO USE GBASIS=ACCT
'''

rhf_runs = []
x = np.linspace(-1.5,1.5,32)
z = np.linspace(-1.0,2.0,32)
for i in range(len(x)):
  for j in range(len(z)):
    myid = "x%2.5f_z%2.5f" % (x[i],z[j])
    struct = Structure(elem=['H','C','N'],
            pos=[[x[i],0,z[j]],
                 [0,0,-CH],
                 [0,0,-(CH+CN)]]
            ,units='A')
    system = generate_physical_system(
      structure  = struct,
      net_charge = 0,
      net_spin   = 0
    )
    rhf_inputs = obj(
        identifier = myid + 'rhf',
        scftyp     = 'rhf', 
        gbasis     = basis,
        ispher     = 1,
        runtyp     = 'energy', 
        coord      = 'unique',
        system     = system,
        dirscf     = True, # use ram for speed
        job        = gms_job,
    )
    rhf = generate_gamess(
        path = myid+'/rhf',
        **rhf_inputs
    )
    rhf_runs.append(rhf)
  # end for j
# end for i

# ======= main ======
if __name__=='__main__':
    from nexus import ProjectManager
    pm = ProjectManager()
    pm.add_simulations(rhf_runs)
    pm.run_project()

    data = []
    for rhf in rhf_runs:
      ga = rhf.load_analyzer_image()
      hline = rhf.input.data.text.split("\n")[2]
      myx,myy,myz = map(float,hline.split()[-3:])
      try:
        Etot = ga.energy["total"]
      except:
        Etot = 0 # if Hartree-Fock fails, then H is too close to C
      # end try
      data.append({
        "x":myx,
        "z":myz,
        "E":Etot
      })
    # end for

    import json
    json.dump(data,open("hf-ccd-grid.json","w"))

# end __main__
