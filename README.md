# SynchronyMeasures-Robustness

This repository contains the scripts and data used in the study:
**"Comparison of Different Spike Train Synchrony Measures Regarding Their Robustness to Erroneous Data from Bicuculline-Induced Epileptiform Activity"**

## About the Study

In this study, various spike train synchrony measures were systematically compared to evaluate their robustness to erroneous data. The data was derived from experiments involving bicuculline-induced epileptiform activity. The results provide insights into the performance and reliability of different synchrony measures under challenging conditions.

## Publication

For full details of the study, refer to the original publication:
**Ciba, M., Bestel, R., Nick, C., De Arruda, G. F., Peron, T., Henrique, C. C., ... & Thielemann, C. (2020).**
*Comparison of Different Spike Train Synchrony Measures Regarding Their Robustness to Erroneous Data from Bicuculline-Induced Epileptiform Activity.*
Neural Computation, 32(5), 887–911.
DOI: [10.1162/neco\_a\_01277](https://doi.org/10.1162/neco_a_01277)
[Link to publication](https://direct.mit.edu/neco/article-abstract/32/5/887/95595/Comparison-of-Different-Spike-Train-Synchrony?redirectedFrom=fulltext)

## Repository Overview

This repository includes:

1. MATLAB scripts for data processing and synchrony measure calculations.
2. Data files for reproducing the analysis.

### Data Structure

* **`MEA-TS-Data.zip`**: Contains spike train data used in the study.
* **`DrCell.zip`**: Includes additional functions required to run the scripts.

## Setup Instructions

To run the scripts (MATLAB 2019a or later recommended):

1. **Unzip Data**:
   Unzip the spike train data file `MEA-TS-Data.zip` into the folder `01_CalcSync`.

2. **Prepare Functions**:

   * Unzip `DrCell.zip` to access all required functions.
   * Add these functions to your MATLAB path:

     * Manually: Use the MATLAB path editor to include the unzipped folder.
     * Programmatically: Run the script `DrCell.m` to add the functions to your path.

3. **Run Scripts**:
   Navigate to `01_CalcSync` and execute the MATLAB scripts as required.

## Citation

If you use the scripts or data in your research, please cite the original publication:
**Ciba, M., Bestel, R., Nick, C., De Arruda, G. F., Peron, T., Henrique, C. C., ... & Thielemann, C. (2020).**
*Comparison of Different Spike Train Synchrony Measures Regarding Their Robustness to Erroneous Data from Bicuculline-Induced Epileptiform Activity.*
Neural Computation, 32(5), 887–911.
DOI: [10.1162/neco\_a\_01277](https://doi.org/10.1162/neco_a_01277)

## License

This project is distributed under the MIT License. See the `LICENSE` file for details.

## Contact

For inquiries or issues related to this repository, please contact:

Manuel Ciba

manuel.ciba@th-ab.de

or 

BioMEMS lab (University of Applied Sciences Aschaffenburg)

biomems.lab.thab@gmail.com
