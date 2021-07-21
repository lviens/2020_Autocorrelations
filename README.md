# Python codes to compute noise autocorrelation functions (ACFs) and image the Kanto Basin, Japan

The following manuscript has been submitted for publication, the initial pre-print version is available [here](https://www.essoar.org/doi/abs/10.1002/essoar.10504543.1). The paper will be resubmitted soon and an updated version will be made available
- Viens L., C. Jiang, and M. Denolle (to be re-submitted), Imaging the Kanto Basin bedrock with earthquake and noise autocorrelation functions

## Description:
* The **Codes** folder contains:
  - **Compute_AC_PWS.py** to compute noise ACFs from the raw MeSO-net data and stack them with the Phase Weighted Stack (PWS, [Schimmel and Paulssen, 1997](https://academic.oup.com/gji/article/130/2/497/760640)) method (2 days of data at the AYMH station are in the E.AYHM.zip file in the **Data** folder).
  - **Reproduce_Fig_1.py** to reproduce Figure 1 (Map of the Kanto Basin).
  - **Reproduce_Fig_2.py** to reproduce Figure 2 of the paper (Noise and earthquake ACF plots).
  - **Reproduce_Fig_3.py** to reproduce Figure 3 of the paper (Maps of the Kanto Basin basin depth).
  - **Reproduce_Fig_4_Noise.py** and **Reproduce_Fig_4_Earthquake.py** to reproduce the Figure 4 of the paper (Axitra simulated ACFs vs oberved ACFs).

* The **Data** folder contains:
  - Five zip files with all the data needed to run the codes and reproduce the figures of the paper. The codes should unzip the zip files automatically.

* The **Figures** folder contains 10 figures that can be plotted with the 4 codes. 


## Codes and their outputs:

* The **Compute_AC_PWS.py** code processes two days of continuous records at the AYMH station and computes the PWS of the 30-min ACFs. For the two days of data, the 20-min ACFs are relatively consistent. Note that we use 1 month of data in the paper. 
<p align="center">
<img src="https://github.com/lviens/2020_Autocorrelations/blob/master/Figures/ACF_example.png" width=50%/>
</p>


* The **Reproduce_Fig_1.py** code outputs a map of the Kanto Basin with the MeSO-net station locations, the JIVSM model, and the four different lines used in this study. The code also plots an azimuthal equidistant projection map centered on the Kanto Basin including all the Mw 6+ earthquakes which occurred within 30 and 95 degrees from the Kanto Basin between May 2017 and 2020.
<p align="center">
<img src="https://github.com/lviens/2020_Autocorrelations/blob/master/Figures/Fig_1.png" width=60%>
</p>

* The **Reproduce_Fig_2.py** code is used to plot the noise and earthquake ACFs along the 3 lines. The Figure below shows the ACFs along Line 1. 
<p align="center">
<img src="https://github.com/lviens/2020_Autocorrelations/blob/master/Figures/Fig_2_Line_1.png" width=60%>
</p>

* The **Reproduce_Fig_3.py** code compares the bedrock depths from the JIVSM, noise ACF, and earthquake ACF measurements. 
<p align="center">
<img src="https://github.com/lviens/2020_Autocorrelations/blob/master/Figures/Fig_3.png" width=60%>
</p>

* The **Reproduce_Fig_4_Noise.py** and **Reproduce_Fig_4_Earthquake.py** codes show a comparison of the Axitra ACF simulations with the observed earthquake and noise ACFs. 
<p align="center">
<img src="https://github.com/lviens/2020_Autocorrelations/blob/master/Figures/Fig_4_Line_1_Earthquake.png" width=60%>
</p>
