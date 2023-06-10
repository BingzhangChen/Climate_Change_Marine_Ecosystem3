2023 Croucher foundation Climate Change and Marine Ecosystems summer course README

## Brief summary

This repository provides the metadata for the several R codes and source data files provided by Bingzhang Chen as part of the summer course.

## Authors

Bingzhang Chen

## Contact email

[bingzhang.chen\@strath.ac.uk](mailto:bingzhang.chen@strath.ac.uk)

## Contributor

Bingzhang Chen is responsible for collecting the data and writing the R code.

## LICENSE

All the codes are covered by the MIT license. Please see the **LICENSE** file for details.

# Metadata

## Software and packages

All the codes are written in R 4.3.0.

## How to run the code

1.  Download the source code and data. The most updated code and data are available at <https://github.com/BingzhangChen/Croucher2023.git> which can be obtained either by using git clone or directly downloaded from Github. It is recommended that the student should go through each line to understand the code.

2.  To run the NPZ model, run **NPZ.R** in R.

3.  To run the linear model of MTE on the data of metabolism in Hatton et al. (2019), run **MTE_LM.R** in R.

4.  To run the nonlinear model to estimate the Temperature Performance Curve (and thermal niche), run **Topt_nls.R** in R.

5.  To run the artificial neural network on picophytoplankton abundances in the South China Sea, run **neuralnet.R** in R.

6.  To run the random forest model on picophytoplankton abundances in the South China Sea, run **RF.R** in R.

## Data files

1.  Metabolism.csv: Raw data of metabolic rate \~ temperature of different groups of organisms in Hatton et al. (2019).

2.  Corrected_Metabolism.Rdata: Metabolic rate data corrected by temperature by Hatton et al. (2019). Once this .Rdata file is loaded in R, a dataframe called **w** can be found with the column headings.

3.  Merged_PHY.Rdata: merged data of phytoplankton growth rate \~ temperature from Chen and Laws (2017) and Kremer et al. (2017) with some additions. Once this .Rdata file is loaded in R, a dataframe called **newdat** can be found with the column headings.

4.  SCS_Pico.Rdata: Chlorophyll and picophytoplankton abundance data in the South China Sea from Chen et al. (2020). Once this .Rdata file is loaded in R, a dataframe called **np** can be found with the column headings.

## Funding

This work has been supported by both the Croucher Foundation and a Leverhulme Trust Research Project Grant (RPG-2020-389).

## References

1.  Chen, B., and E. A. Laws. 2017. Is there a difference of temperature sensitivity between marine phytoplankton and heterotrophs? Limnology and Oceanography 62:806--817.

2.  Chen, B., H. Liu, W. Xiao, L. Wang, and B. Huang. 2020. A machine-learning approach to modeling picophytoplankton abundances in the South China Sea. Progress in Oceanography, 189:102456.

3.  Hatton, I. A., A. P. Dobson, D. Storch, E. D. Galbraith, and M. Loreau. 2019. Linking scaling laws across eukaryotes. Proceedings of the National Academy of Sciences 116: 21616--21622. <doi:10.1073/pnas.1900492116>

4.  Kremer, C. T., M. K. Thomas, and E. Litchman. 2017. Temperature- and size-scaling of phytoplankton population growth rates: Reconciling the Eppley curve and the metabolic theory of ecology. Limnology and Oceanography 62:1658--1670.
