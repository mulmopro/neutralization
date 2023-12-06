# Neutralization in square/circular cross-sectional Tmixer
3D solver that models dilution and neutralization reactions (in this case is tailored for NaOH and HCl)

twoLiquidMixingFoamNeutralization: Please compile first "Cphi" class, then "chemReac" class and, at the end, the executable solver "twoLiquidMixingFoamNeutralization".

testCase: please run ./Allrun (this command generates the circular cross-sectional Tmixer mesh and it execute the solver "twoLiquidMixingFoamNeutralization")

# Experimental setup
The experimental setup is reported:
![image](https://github.com/mulmopro/neutralization/assets/102947817/dece5b56-0b70-4d40-bc49-6329281196e9)
Two geometries (D<sub>h</sub> = 3 mm) are investigated both experimentally and numerically: the former is a square cross-sectional T-mixer (T<sub>sq</sub>), the latter is a circular cross-sectional T-mixer (T<sub>cr</sub>)

# Reference
Should you use these codes for your own research, we would greatly appreciate if you could cite our work:

G. Battaglia, S. Romano, A. Raponi, F. Volpe, L. Bellanca, M. Ciofalo, D. Marchisio, A. Cipollina, G. Micale, A. Tamburini,
Mixing phenomena in circular and rectangular cross-sectional T-mixers: Experimental and numerical assessment,
Chemical Engineering Research and Design,
Volume 201,
2024,
Pages 228-241,
ISSN 0263-8762,
https://doi.org/10.1016/j.cherd.2023.11.056.
(https://www.sciencedirect.com/science/article/pii/S0263876223007669)

## BibTeX entry
@article{BATTAGLIA2024228,
title = {Mixing phenomena in circular and rectangular cross-sectional T-mixers: Experimental and numerical assessment},
journal = {Chemical Engineering Research and Design},
volume = {201},
pages = {228-241},
year = {2024},
issn = {0263-8762},
doi = {https://doi.org/10.1016/j.cherd.2023.11.056},
url = {https://www.sciencedirect.com/science/article/pii/S0263876223007669},
author = {G. Battaglia and S. Romano and A. Raponi and F. Volpe and L. Bellanca and M. Ciofalo and D. Marchisio and A. Cipollina and G. Micale and A. Tamburini},
keywords = {T-mixer, Digital image analysis, RANS k- ε tuning, Mixing, Asymmetric flow rates}
}
