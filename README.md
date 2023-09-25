# Optimization of RIS-aided MIMO Systems via the Cutoff Rate

This repo is cloned from [this CodeOcean capsule](https://codeocean.com/capsule/1300715/tree/v1), which contains the code for the following scientific paper:

Nemanja Stefan Perović, Le-Nam Tran, Marco Di Renzo, and Mark F. Flanagan, "Optimization of RIS-aided MIMO Systems via the Cutoff Rate," IEEE Wireless Commun. Lett. ,  vol. 10, no.8, pp. 1692-1696, Aug. 2021.

## Instructions
"**Algorithm1PGM.m**" and "**Algorithm2SCA.m**" implement the projected PGM and SCA method proposed in the above paper.

Note that, in this implementation of the SCA method, we use [Yalmip](https://yalmip.github.io/) as a parser (instead of CVX) as it has been found to achieve better efficiency. We still leave the CVX code in the implementation for reference.
