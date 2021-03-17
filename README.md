# Square-well-range-R-a12-reff-
Accompanying code to:  https://journals.aps.org/pra/abstract/10.1103/PhysRevA.99.023618
preprint in: https://arxiv.org/abs/1811.04436

Contains a python script "find_R_from_a12_reff.py" where one can determine the diameter of attractive square well from the input s-wave scattering length and effective range reff. Works only for reff>0, a12<0, for a potential which does not support a two-body bound state. In a script "quantum_bose-bose_mixture_functional.py", one can also calculate energy per particle for a given density, s-wave scattering length a12<0 and effective range reff>0 of the attractive interactions.  
In both scripts, density is in units $a_{11}^{-3}$, lengths in units $a_{11}$, and energy in $\hbar^2 / (2 m a_{11}^2)$.
