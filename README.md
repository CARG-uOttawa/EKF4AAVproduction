
# Monitoring the recombinant Adeno-Associated Virus Production by Extended Kalman Filter

Adeno-associated virus (AAV) is a flexible viral vector technology for gene therapy applications 
that holds promise as the safest and most effective way to repair single-gene abnormalities in nondividing cells. It represents one of the most straightforward vehicles for direct translation of genomic revolution in medicine. However, improving the viral titer productivity in the AAV production remains challenging.
The fist step in this direction is monitoring the state variables of AAV production efficiently, since the incorporation of such variables in process operation can provide improved control performance with enhanced productivity. But the traditional approaches to monitoring are expensive, laborious and time-consuming.
This paper presents Extended Kalman Filter (EKF) approach to reduce the errors in measurements of viable cell density (Xv) and to estimate the others state variables (GLC, GLN, LAC, AMM and AAV Viral titer) of AAV production that are measured at low sampling frequency. 
The proposed EKF uses an Unstructured Mechanistic Kinetic Model (UMKM) that depend of Xv and it is applicable in Cell Expansion Phase (CEP) and Viral Vector Production Phase (VVPP) phases of upstream process. Three datasets were used to parameters estimation, calibration and test of the proposed EKF, and the data was collected from production of rAAV by triple plasmid transfection of HEK293SF-3F6 cells.
The parameters used in the UMKM were estimated with Neural Ordinary Differential Equation and Bayesian inference approaches and were the used as initial parameters in EKF.
The evaluation showed that the proposed approach was able to estimate the GLC, LAC and AAV viral titer efficiently beside reduce the noise of Xv.
It has a high potential to be extended to a online soft-sensor or be a sub-model of a model predictive control and be classified as a low cost and fast approach to monitoring AAV production. 



# Contents

* Matlab code of EKF implemantation
* Julia code NODE and Bayesian inference implementation
* dataset
