context("lagr")

test_that("str_length is number of characters", {
    data(poverty)
    df = poverty[poverty$year=="1990",]
    
    #Define which variables we'll use as predictors of poverty:
    f = logitindpov ~ pag + pex + pman + pserve + pfire + potprof
    
    #Lasso model
    #tune = lagr.sel(formula=f, data=df, coords=c('x', 'y'), longlat=TRUE, varselect.method='AICc', kernel=epanechnikov, tol.bw=0.01, bw.type='knn', bwselect.method='AICc', resid.type='pearson')
    model = lagr(formula=f, data=df, coords=c('x','y'), longlat=TRUE, varselect.method='AICc', bw=0.2, kernel=epanechnikov, bw.type='knn', resid.type='pearson', verbose=TRUE)
    
})

