# xsec
is a module that contains all the functions necessary to simulate the scattering cross section & efficiency for an elastic sphere embedded in another solid

The code is based on
[1] Ying and Truell, Journal of Applied Physics 27, 1086 (1956) - Longitudinal Scattering
[2] Iwashimazu, Journal of Sound and Vibration, 40(2), 267-271 (1975) - Transverse Scattering

The main function is:
sigma_sphere(k,p,MatParams):
    """
    Returns:
        sigma: scattering cross section
        scat_eff:  
    Inputs:
        k: a vector of wavenumbers
        p: polarization index (p=1 (longitudinal), p=2,3 (transverse))
        MatParams: Object containing the materials properties
    """
    
 There is an important helper object "MatParams" which packs in all the materials properties for easy interfacing.  
 Call example:
  anp = 3e-9; # particle radius
  C11m = 166e9; C44m = 80e9;  rhom = 2300; # matrix properties
  C11np = 126e9; C44np = 68e9; rhonp = 5323; # particle properties
  MP = xsec.create_matparams(anp,C11m,C44m,rhom,C11np,C44np,rhonp); # <--- MP object can be passed to functions in xsec module
  
  # Call_Truell 
  contains an example code that uses the xsec module
 
