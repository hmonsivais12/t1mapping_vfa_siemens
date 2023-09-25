# t1mapping_vfa_siemens
Calculates a T1 map with B1 correction using .nii images from two flip angles. Tested on Siemens Prisma vendor sequence. 

Usage:
         T1_B1cor=T1cal_nii;

The output will be T1_B1cor.nii
B1 correction is done by smoothing FA3 20 times.
--- written by Humberto Monsivais Jan. 12, 2022 -------------
--- based on Dr. Chien-Lin's code for T1 mapping --------

Separated the individual functions for easier reading and changed the script to  accept .nii files from the Siemens scanner.

Please have SPM installed in MATLAB====

User selects .nii files for the two flip angles used in increasing order and also the B1 map file that has already been coregistered to the T1w images
