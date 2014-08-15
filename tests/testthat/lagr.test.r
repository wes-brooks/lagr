context("lagr")

test_that("str_length is number of characters", {
    data(poverty)
    df = poverty[poverty$year=="1970",]
    
    #Define which variables we'll use as predictors of poverty:
    f = logitindpov ~ pag + pex + pman + pserve + pfire + potprof
    
    #LAGR model
    tune = lagr.tune(formula=f, data=df, coords=c(x,y), longlat=TRUE, varselect.method='AIC', kernel=epanechnikov, tol.bw=0.01, bw.type='knn', bwselect.method='AIC', verbose=TRUE)
    model = lagr(formula=f, data=df, coords=c('x','y'), longlat=TRUE, varselect.method='AIC', bw=0.2, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
})

