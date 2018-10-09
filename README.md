# VTI_AVO
Modeling impact of VTI on AVO  
##  CONCEPT  
Using the AVO equations the program generates a model of the amplitude variation with offset of the given parameters. The required paramters are the Vp, Vs, and Rhob of both the top layer and bottom layer. If you input Thomsen parameters, i.e. epsilon and delta the model computes the resulting AVO response using Ruger equations. 
##  COMMAND LINE INTERFACE  
```
>python vtiavo.py  -h
usage: vtiavo.py [-h] [--vp VP VP] [--vs VS VS] [--rb RB RB]
                 [--angles ANGLES ANGLES] [--epsilonp EPSILONP EPSILONP]
                 [--epsilons EPSILONS EPSILONS] [--delta DELTA DELTA]
                 [--model {iso,vti,vva_delta,vva_epes}] [--hideplot]

AVO w/ VTI modeling and inversion Mar 6, 2015

optional arguments:
  -h, --help            show this help message and exit
  --vp VP VP            Vp top bottom. dfv= 9500 10000
  --vs VS VS            Vs top bottom. dfv= 4700 5200
  --rb RB RB            Rhob top bottom. dfv= 2.70 2.65
  --angles ANGLES ANGLES
                        Minimum Maximum angles.dfv= 0 60 deg
  --epsilonp EPSILONP EPSILONP
                        Epsilon P of upper and lower layer. dfv= 0.10 0.20
  --epsilons EPSILONS EPSILONS
                        Epsilon S of upper and lower layer. dfv= 0.05 0.08
  --delta DELTA DELTA   Delta Delta of upper and lower layer. dfv= 0.10 0.20
  --model {iso,vti,vva_delta,vva_epes}
                        modeling options:isotropic,anisotropic VTI, Phase
                        velocities with angle epsilon, delta, Phase velocities
                        with angle epsilonP epsilonS -> delta is calculated
                        .dfv= iso
  --hideplot            Do not display plots, only save to pdf  
  ```  
  No file is needed to run the program. The default is to run an isotropic model and generate a plot of an AVO curve. If you want to run VTI, then simply us ``--model vti``  should use the default paramters and generate plots. 
  
  >  Will update more details later.
