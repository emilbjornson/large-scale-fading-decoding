Large-Scale-Fading Decoding in Cellular Massive MIMO Systems with Spatially Correlated Channels
==================

This is a code package is related to the follow scientific article:

Trinh Van Chien, Christopher Mollén and Emil Björnson, “[Large-Scale-Fading Decoding in Cellular Massive MIMO Systems with Spatially Correlated Channels](https://arxiv.org/abs/1807.08071),” IEEE Transactions on Communications, To appear.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

Massive multiple-input–multiple-output (MIMO) systems can suffer from coherent intercell interference due to the phenomenon of pilot contamination. This paper investigates a two-layer decoding method that mitigates both coherent and non-coherent interference in multi-cell Massive MIMO. To this end, each base station (BS) first estimates the channels to intra-cell users using either minimum mean-squared error (MMSE) or element-wise MMSE (EW-MMSE) estimation based on uplink pilots. The estimates are used for local decoding on each BS followed by a second decoding layer where the BSs cooperate to mitigate inter-cell interference. An uplink achievable spectral efficiency (SE) expression is computed for arbitrary two-layer decoding schemes. A closed-form expression is then obtained for correlated Rayleigh fading, maximum-ratio combining, and the proposed large-scale fading decoding (LSFD) in the second layer. We also formulate a sum SE maximization problem with both the data power and LSFD vectors as optimization variables. Since this is an NP-hard problem, we develop a low-complexity algorithm based on the weighted MMSE approach to obtain a local optimum. The numerical results show that both data power control and LSFD improve the sum SE performance over single-layer decoding multi-cell Massive MIMO systems.


## Content of Code Package

The article contains 8 simulation figures, numbered 3-10. Simulation_Fig3_4_5.m generates Figures 3-5, Simulation_Fig6_7.m generates Figures 6-7, Simulation_Fig8_9.m generates Figures 8-9, and Simulation_Fig10.m generates Figure 10. The package also contains 13 Matlab functions that are used by some of the scripts.

See each file for further documentation.


## Acknowledgements

This paper was supported by the European Union’s Horizon 2020 research and innovation programme under grant agreement No 641985 (5Gwireless). It was also supported by ELLIIT and CENIIT.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
