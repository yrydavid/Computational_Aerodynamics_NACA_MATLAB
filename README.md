# Computational_Aerodynamics_NACA_MATLAB
# Computational Aerodynamics: NACA Airfoil Generator & Hess-Smith Panel Method

A complete MATLAB framework for low-speed aerodynamics analysis, integrating parametric airfoil geometry generation, classical Thin Airfoil Theory (Glauert), and the Hess-Smith Panel Method. Developed as part of the Aerodynamics academic coursework at the National University of Science and Technology POLITEHNICA Bucharest (UPB).

## Project Structure & Modules
1. **Geometry Generation (`generate_naca_profile.m`):** Parametric generation for NACA 4-digit and 5-digit series (e.g., NACA 44015, NACA 23015) using cosine spacing discretization to refine density near leading and trailing edges.
2. **Thin Airfoil Theory (`thinAirfoilTheory.m`):** Implementation of Glauert's integral equations using Fourier series expansion ($A_0$, $A_1$, $A_2$) to compute theoretical lift ($C_L$) and pitching moment ($C_{m,c/4}$).
3. **Hess-Smith Panel Method (`hessSmith_academic.m`):** Source-vortex singularity distribution, global linear system assembly, Kutta condition enforcement at the trailing edge, and pressure coefficient ($C_p$) distribution extraction.

## Validation & Results
All numerical models are cross-validated against experimental wind tunnel data extracted from standard aerospace references.

## References & Bibliography
- **Abbott, I. H., & von Doenhoff, A. E.** (1959). *Theory of Wing Sections, including a summary of airfoil data*. Dover Publications.
- **Katz, J., & Plotkin, A.** (2001). *Low-Speed Aerodynamics: From Vortex Theory to Numerical Methods*. Cambridge University Press.
- **UPB Aerospace Engineering Faculty** – *Aerodynamics Laboratory Course Notes & Guidelines*.
