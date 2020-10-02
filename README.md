# Python codes to compute noise autocorrelation functions (ACFs) and image the Kanto Basin, Japan

The manuscript has been submitted for publication, a pre-print version is available here.

## Description:
* The **Codes** folder contains:
  - **Compute_AC_PWS.py** to compute noise ACFs from the raw MeSO-net data and stack them with the Phase Weighted Stack (PWS, [Schimmel and Paulssen, 1997](https://academic.oup.com/gji/article/130/2/497/760640)) method (2 days of data at the AYMH station are in the E.AYHM.zip file in the **Data** folder).
  - **Reproduce_Fig_1.py** to reproduce Figure 1 (Map of the Kanto Basin).
  - **Reproduce_Fig_2_3.py** to reproduce Figures 2 and 3 of the paper (Noise and earthquake ACF plots).
  - **Reproduce_Fig_4.py** to reproduce Figure 4 of the paper (Maps of the Kanto Basin basin depth).

* The **Data** folder contains:
  - Four zip files with all the data needed to run the codes and reproduce the figures of the paper. The codes should unzip the zip files automatically.

* The **Figures** folder contains the 7 figures that can be plotted with the 4 codes. 


## Example:

* The **Compute_AC_PWS.py** code processes two days of continuous records at the AYMH station and computes the PWS of the 30-min ACFs. For the two days of data, the 20-min ACFs are relatively consistent. Note that we use 1 month of data in the paper. 
<p align="center">
<img src="https://github.com/lviens/2020_Autocorrelations/blob/master/Figures/ACF_example.png" width=50%/>
</p>


* The **Reproduce_Fig_1.py** code outputs a map of the Kanto Basin with the MeSO-net station locations, the JIVSM model, and the four different lines used in this study. The code also plots an azimuthal equidistant projection map centered on the Kanto Basin including all the Mw 6+ earthquakes which occurred within 30 and 95 degrees from the Kanto Basin between May 2017 and 2020.
<p align="center">
<img src="https://github.com/lviens/2020_Autocorrelations/blob/master/Figures/Fig_1.png" width=60%>
</p>

* The **Reproduce_Fig_2_3.py** is used to plot the noise and earthquake ACFs along the 4 lines. The Figure below shows the ACFs along Line 1. 
<p align="center">
<img src="https://github.com/lviens/2020_Autocorrelations/blob/master/Figures/Line_1.png" width=60%>
</p>

* The **Reproduce_Fig_4.py** compares the bedrock depths from the JIVSM, noise ACF, and earthquake ACF measurements. 
<p align="center">
<img src="https://github.com/lviens/2020_Autocorrelations/blob/master/Figures/Fig_4.png" width=60%>
</p>
