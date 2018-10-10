# VTI_AVO
Modeling impact of VTI on AVO  
##  CONCEPT  
Using the AVO equations the program generates a model of the amplitude variation with offset of the given parameters. The required paramters are the Vp, Vs, and Rhob of both the top layer and bottom layer. If you input Thomsen parameters, i.e. epsilon and delta the model computes the resulting AVO response using Banik equations. 
##  COMMAND LINE INTERFACE  
```
python vtiavo.py -h
usage: vtiavo.py [-h] [--vp VP VP] [--vs VS VS] [--rb RB RB]
                 [--angles ANGLES ANGLES] [--epsilonp EPSILONP EPSILONP]
                 [--epsilons EPSILONS EPSILONS] [--delta DELTA DELTA]
                 [--model {iso,vti,vva_delta,vva_epsilon}] [--hideplot]

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
  --model {iso,vti,vva_delta,vva_epsilon}
                        modeling options:isotropic,anisotropic VTI, Phase
                        velocities with angle epsilon, delta, Phase velocities
                        with angle epsilonP epsilonS -> delta is calculated
                        .dfv= iso
  --hideplot            Do not display plots, only save to pdf
  ```  
  No file is needed to run the program. The default is to run an isotropic model and generate a plot of an AVO curve. If you want to run VTI, then simply use ``--model vti``  should use the default paramters and generate plots. 
  
 Most of the equations are taken from the book *__Seismic Reflection Processing With Special Reference to Anisotropy__*  Springer Verlag, 2004.  
 Simple AVO equations supply the variation of reflectivity with angle of incidence in the isotropic case. If Thomsen parameters are available then you can use the epsilon and delta of Thomsen to compute the VTI reflectivity change with angle of incidence. 

 Using Banik definitions, epsilon p and epsilon s are required. The former is equivalent to that of Thomsen, whereas the latter is similarly defined but for the SV (at an angle of 45 degrees, as opposed to the horizontal, i.e. 90 degrees for the epsilon p).  
 
 The epsilon p and epsilon s have and impact on the velocities variation with angle. Banik's delta, has an impact on the amplitudes of both the P and SV. The equations indicate that the magnitude of the differences between layer 2 vs layer 1 are the main contributors to the anisotropic behaviour as seen on both the velocities and amplitudes.  
 
 Run the program with the ``--model vti`` option to generate the amplitude variation with angle and plot the proposed intercept and gradient line fitted to the squared of the angle of incidence.  
 
 Run the program with the ``--model vva_delta`` to model the impact of the Banik delta on the amplitudes of the PP and the PSV.  
 
 Run the program with the  ``--model vva_epsilon`` to model the impact of the Banik epsilon p and epsilon s on the phase velocities with incidence angle.  
 
 All the generated plots are saved as pdf files.  
 
 
 
 