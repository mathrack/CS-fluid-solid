The present data and the associated source code are freely available under the GNU GPL v3 licence.

# Direct access to our LES results (Microsoft Excel files)

[Reynolds number based on the friction velocity = 150, Prandtl = 0.71](/../raw/master/LES_database/wale_150_pr071/wale_150_pr071.xlsm)

[Reynolds number based on the friction velocity = 395, Prandtl = 0.71](/../raw/master/LES_database/wale_395_pr071/wale_395_pr071.xlsm)

[Reynolds number based on the friction velocity = 395, Prandtl = 1](/../raw/master/LES_database/wale_395_pr1/wale_395_pr1.xlsm)

[Reynolds number based on the friction velocity = 1020, Prandtl = 0.71](/../raw/master/LES_database/wale_1020_pr071/wale_1020_pr071.xlsm)

[Reynolds number based on the friction velocity = 1020, Prandtl = 1](/../raw/master/LES_database/wale_1020_pr1/wale_1020_pr1.xlsm)

# Description of the repository

The folder [LES_database](/LES_database/) contains the data and source code used to perform LES of the turbulent channel flow with conjugate heat transfer using *Code_Saturne*.

The folder [LES_VnV](/LES_VnV/) contains the data and source code used to verify and validate the internal coupling module implemented in *Code_Saturne*.

Some subsets of the present work were presented at the international conferences NURETH-17 (Xi'An, China) and NENE2017 (Bled, Slovenia) during September 2017.
The corresponding proceedings are available on [HAL](https://hal.archives-ouvertes.fr/hal-01631515) and at [https://repo.ijs.si/CFLAG/NENE2017](https://repo.ijs.si/CFLAG/NENE2017).
The authors have detected and corrected some errors in the scilab post-processing scripts after publication of the proceedings.
To be more specific, the error estimator proposed is incorrect: the reconstructed statistics at the wall automatically satisfy the compatibility condition.
The associated error is thus zero.

# Acknowledgements

The author and coworkers thank the Slovenian Research Agency and EDF R&D for funding the study (research projects P2-0026 and PR-07184).
We also thank framasoft and the Institut Jozef Stefan for providing the gitlab services that host the present project at [https://repo.ijs.si/CFLAG/CS-Fluid-Solid](https://repo.ijs.si/CFLAG/CS-Fluid-Solid) and [https://framagit.org/CFLAG/CS-fluid-solid](https://framagit.org/CFLAG/CS-fluid-solid).
