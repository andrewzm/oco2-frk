# On statistical approaches to generate Level 3 products from remote sensing retrievals

Reproducible code for generating the results in Sections 2 and 3 of the article `On statistical approaches to generate Level 3 products from remote sensing retrievals' by Zammit-Mangion, A., Cressie, N., and Shumack, C.'

![alt text](https://raw.githubusercontent.com/andrewzm/oco2-frk/master/img/FRK_OCO2E_no_sub.png)

*Figure 8 (e): The difference between the FRK Version 8r predictions and those of FRK Version 7r for 13 May 2016.*

## Instructions

For reproducing the results in Section 2 please run the `R` script `Section2_Running_example.R`. 

For reproducing the Level 3 products and the results in Section 3 you will need to run the scripts in the three folders. 

**1\_Generate\_FRK\_L3**: Scripts in this folder download the V7 and V8 data, generate the Level 3 products, and filter the product at the TCCON locations for comparison later.

**2_Filter\_TCCON**: Scripts in this folder download the TCCON data and subset them according the overpass time of OCO-2.

**3_Analyse\_Results**: Scripts in this folder are used to generate the comparison tables and figures in Section 3.


## Other resources

The FRK products can be downloaded directly from [here](https://niasra.uow.edu.au/cei/oco2level3/index.html).

A pre-print version of the paper will be available shortly.

If you have any queries please do not hesitate to [e-mail](mailto:azm@uow.edu.au) me.