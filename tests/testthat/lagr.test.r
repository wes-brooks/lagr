context("lagr")

test_that("str_length is number of characters", {
    data(poverty)
    df = poverty[poverty$year=="1970",]
    
    #Define which variables we'll use as predictors of poverty:
    f = logitindpov ~ pag + pex + pman + pserve + pfire + potprof
    fl = rbind(c(-89.4437, 38.8863), c(-89.5350, 41.4068))
    colnames(fl) = c('x', 'y')
    
    #LAGR model
    tune = lagr.tune(formula=f, data=df, coords=c(x,y), longlat=TRUE, varselect.method='AIC', kernel=epanechnikov, tol.bw=0.01, bw.type='knn', bwselect.method='AIC', verbose=TRUE)
    model = lagr(formula=f, data=df, coords=c('x','y'), fit.loc=fl, longlat=TRUE, varselect.method='AIC', bw=0.2, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
    
    y2 = ifelse(df$logitindpov > mean(df$logitindpov), 1, 0)
    df2 = df
    df2$logitindpov = y2
    model = lagr(formula=f, data=df2, family=binomial(), coords=c('x','y'), fit.loc=fl, longlat=TRUE, varselect.method='AIC', bw=0.5, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
    model = lagr(formula=f, data=df2, family=gaussian, coords=c('x','y'), fit.loc=fl, longlat=TRUE, varselect.method='AIC', bw=0.5, kernel=epanechnikov, bw.type='knn', verbose=TRUE)
    
    sglm = SGL(list(x=as.matrix(df2[,c(3:7,15)]), y=df2$logitindpov), type='logit', alpha=0, index=1:6)
})

