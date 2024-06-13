# Greenness Of Trees - FDA

A functional data analysis of the greenness of trees for 200 locations, over 1 year at 8-day intervals. To run the R script, one just has to set the working directory to the csv file provided. The R script falls into 8 parts:

#### 1 - Visualising the raw data

#### 2- Conduct a Generalised Cross Validation (GCV) of the data with a saturated B-spline basis and standard roughness penalty

#### 3 - Conduct a GCV of the data with a saturated B-spline basis and a harmonic acceleration roughness penalty

#### 4 - Obtain the first and second derivatives of the data

#### 5 - Perform an unpenalised functional Principal Component Analysis (fPCA) of the data

#### 6 - Use GCV to find the best smoothing parameter for the principal components

#### 7 - Use Leave One Out Cross Validcation (LOOCV) to find the best smoothing parameter for the principal components

#### 8 - Check orthogonality of unpenalised and penalised principal components (using the inprod function)
