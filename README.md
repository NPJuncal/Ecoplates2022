# Ecoplates2022

R code and raw and processed data of the manuscript " Coupling microbial incubations and FTIR-ATR to assess organic matter compounds persistence in salt marsh soils"

## Repository structure

`data/`: folder containing raw and processed data from ecoplates incubations, FTIR spectrums and geochemical variables

-   **`T0.csv`**, **`T1.csv`**, **`T2.csv`**, **`T3.csv`** and **`T4.csv`**: .csv files containing absorbance data from the ecoplates at times 0, 1, 2, 3 and 4.
-   **`FT1.csv`**, **`FT2.csv`**, **`FT3.csv`**, and **`FT4.csv`**: .csv files containing processed data from the ecoplates at times 1, 2, 3 and 4.
-   **`Cores.csv`**: .csv containing geochemical data from the cores (OM, CaCO3 and granulometry)
-   **`Spectra.Ecop.csv`**: .csv containing the FTIR-ATR spectrums of the samples and the Ecoplates compounds

`scripts/`: folder containing the R scripts

-   **`Ecoplates2022_4.0.R`**: R script containing the code to process Ecoplates data and produce statistical analysis and figures
-   **`FTIR.Ecoplates.R`**: R script containing the code to process FTIR-ATR spectrums and the transpose matrix PCA and global PCA, statistical analysis and figures
-   **`FTIR.Ecoplates.R`**: R script containing the code to process geochemical and produce statistical analysis and figures

`results/`: folder containing .csv files with PCA results and plots and figures

## Licence

[GNU General Public License GPL (\>= 3)](https://www.gnu.org/licenses/gpl-3.0.html)

## Contact

For any queries, please contact Nerea Piñeiro Juncal (np.juncal [at] gmail.com).
