#include <Rcpp11>
#include <iostream>
#include <math.h>
#include <numeric>

using namespace Rcpp;
using namespace std;



double loglik(NumericVector eta, NumericVector y, NumericVector weights, Function linkinv, Function devfun)
{
    try {
        // Compute some values we need for the logLik computation:
        NumericVector mu = linkinv(eta);

        // Calculate the actual log likelihood:
        NumericVector dev_resids = devfun(y, mu, weights);
        double ll = sum(dev_resids);
        return(ll);
    }
    catch (...) {
        cout << "Error in loglik!" << endl;

        cout << "eta: ";
        for( auto c : eta)
            cout << c << ' ';
        cout << endl;

        cout << "y: ";
        for( auto c : y)
            cout << c << ' ';
        cout << endl;

        cout << "weights: ";
        for( auto c : weights)
            cout << c << ' ';
        cout << endl;
    }
    return 0.;
}


void gradCalc(NumericVector eta, NumericVector y, NumericVector weights, Function linkinv, NumericVector ldot)
{
    try {
        // Compute some values we need for the gradient computation:
        NumericVector mu = linkinv(eta);
        double sumw = sum(weights);
    
        // Calculate the actual log likelihood:
        ldot = weights * (mu-y) / sumw;
    }
    catch (...) {
        cout << "Error in gradCalc!" << endl;

        cout << "eta: ";
        for( auto c : eta)
            cout << c << ' ';
        cout << endl;

        cout << "y: ";
        for( auto c : y)
            cout << c << ' ';
        cout << endl;

        cout << "weights: ";
        for( auto c : weights)
            cout << c << ' ';
        cout << endl;

        cout << "ldot: ";
        for( auto c : weights)
            cout << c << ' ';
        cout << endl;
    }

}


void rcppLinSolver(NumericMatrix X, NumericVector y, NumericVector w, NumericVector adaweights, int nrow, int ncol, int numGroup, NumericMatrix beta, Function linkinv, Function devfun, IntegerVector rangeGroupInd, IntegerVector groupLen, NumericVector lambda, int step, int innerIter, double thresh, NumericVector ldot, NumericVector nullBeta, double gamma, NumericVector eta, IntegerVector betaIsZero, int& groupChange, IntegerVector isActive, IntegerVector useGroup, int reset)
{
    NumericVector theta(ncol);
    int startInd = 0;
    double zeroCheck = 0;
    double check = 0;
    int count = 0;
    double diff = 1;
    double norm = 0;
    double uOp = 0;
    double Lnew = 0;
    double Lold = 0;
    double sqNormDelta = 0;
    double iProd = 0;
    NumericVector etaNew(nrow);
    NumericVector etaNull(nrow);
    NumericVector var(nrow);
    
    for(int i = 0; i < numGroup; i++)
    {
        if(useGroup[i] == 1)
        {
            startInd = rangeGroupInd[i];
            
            // Setting up null residuals calc to check if group is 0
            // Null residuals means the coefficients of group i are set to zero before computing residuals.
            for(int k = 0; k < nrow; k++)
            {
                etaNull[k] = eta[k];
                for(int j = startInd; j < rangeGroupInd[i] + groupLen[i]; j++)
                {
                    etaNull[k] = etaNull[k] - X[k + nrow * j] * beta[step*ncol + j]; 
                }
            }
            gradCalc(etaNull, y, w, linkinv, ldot);
            
            // Now compute the correlation of this group with the null residuals.
            NumericVector grad(groupLen[i]);
            for(int j = 0; j < groupLen[i]; j++)
            {
                grad[j] = 0;
                for(int k = 0; k < nrow; k++)
                {
                    grad[j] = grad[j] + X[k + nrow * (j + rangeGroupInd[i])] * ldot[k];
                }
            }
            zeroCheck = sum(pow(grad,2));
            
            // Is this group's correlation with the null residuals smaller than the threshold?
            // If so, set the coefficients of the group to zero.
            if(zeroCheck <= pow(adaweights[i],2)*pow(lambda[step],2)*groupLen[i])
            {
                if(betaIsZero[i] == 0)
                {
                    for(int k = 0; k < nrow; k++)
                    {
                        for(int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++)
                        {
                            eta[k] = eta[k] - X[k + nrow * j] * beta[step*ncol + j];
                        }
                    }
                }
                betaIsZero[i] = 1;
                for(int j = 0; j < groupLen[i]; j++)
                {
                    beta[step*ncol + j + rangeGroupInd[i]] = 0;
                }
            }
            else // The group's coefficients are nonzero
            {
                if(isActive[i] == 0)
                {
                    groupChange = 1;
                }
                isActive[i] = 1;
                
                // Hold the previous values of the coefficients
                for(int k = 0; k < ncol; k++)
                {
                    theta[k] = beta[step*ncol + k];
                }
                
                betaIsZero[i] = 0;
                NumericVector z(groupLen[i]);
                NumericVector U(groupLen[i]);
                NumericVector delta(groupLen[i]);
                NumericVector betaNew(ncol);
                
                count = 0;
                check = 100000;
                double t = 1;
                
                // Until convergence, iteratively recompute the group's coefficients:
                while(count <= innerIter && check > thresh)
                {
                    count++;
                    
                    // Compute the working residuals (i.e. residuals at the current coefficient values)
                    gradCalc(eta, y, w, linkinv, ldot);
                    
                    // Compute the group's correlation with the working residuals
                    for(int j = 0; j < groupLen[i]; j++)
                    {          
                        grad[j] = 0;
                        for(int k = 0; k < nrow; k++)
                        {
                            grad[j] = grad[j] + X[k + nrow * (j + rangeGroupInd[i])] * ldot[k];
                        }
                    }
                    
                    // Compute the log likelihood for the previous coefficient iteration
                    Lold = loglik(eta, y, w, linkinv, devfun);
                    
                    // Back-tracking to optimize the step size for gradient descent:
                    diff = -1;
                    int optim_steps = 0;
                    while(diff < 0)
                    {
                        optim_steps++;
                        
                        for(int j = 0; j < groupLen[i]; j++)
                        {
                            z[j] = beta[step*ncol + j + rangeGroupInd[i]] - t * grad[j];
                        }
                        
                        norm = sum(pow(z, 2));
                        norm = sqrt(norm);
                        
                        if(norm != 0){
                            uOp = (1 - adaweights[i]*lambda[step]*sqrt(double(groupLen[i]))*t/norm);   //Or not?
                        }
                        else{uOp = 0;}
                        
                        if(uOp < 0)
                        {
                            uOp = 0;
                        }
                     
                        for(int j = 0; j < groupLen[i]; j++)
                        {
                            U[j] = uOp*z[j];
                            delta[j] = U[j] - beta[step*ncol + j + rangeGroupInd[i]];    
                        }
                     
                        // Setting up betaNew and etaNew in direction of Grad for descent momentum
                        for(int k = 0; k < nrow; k++)
                        {
                            etaNew[k] = eta[k];
                            for(int j = 0; j < groupLen[i]; j++)
                            {
                                etaNew[k] = etaNew[k] + delta[j] * X[k + nrow*(rangeGroupInd[i] + j)];
                            }
                        }
                        
                        // Compute log likelihood for the working coefficients
                        Lnew = loglik(etaNew, y, w, linkinv, devfun);
                        
                        sqNormDelta = sum(pow(delta, 2));
                        iProd = sum(grad * delta);
                        
                        // Check whether the working coefficients reduce the log likelihood
                        diff = Lold - Lnew + iProd + 0.5/t * sqNormDelta;
                        
                        // Reduce the step size for the next iteration
                        t = t * gamma;
                    }
                    t = t / gamma;
                    
                    check = 0;
                    
                    for(int j = 0; j < groupLen[i]; j++)
                    {
                        // Check the difference between iterations of the coefficients
                        check = check + fabs(theta[j + rangeGroupInd[i]] - U[j]);

                        // Calculate the null etas (i.e. the etas when group i's coefficients are zero)
                        for(int k = 0; k < nrow; k++)
                        {
                            eta[k] = eta[k] - X[k + nrow * (j + rangeGroupInd[i])]*beta[step*ncol + j + rangeGroupInd[i]];
                        }
                        
                        //Apply Nesterov acceleration to calculate the new coefficient values:
                        beta[step*ncol + j + rangeGroupInd[i]] = U[j] + count/(count+3) * (U[j] - theta[j + rangeGroupInd[i]]);
                        theta[j + rangeGroupInd[i]] = U[j];
                        
                        // Compute the new values of eta after iterating group i's coefficients.
                        for(int k = 0; k < nrow; k++)
                        {
                            eta[k] = eta[k] + X[k + nrow * (j + rangeGroupInd[i])]*beta[step*ncol + j + rangeGroupInd[i]];
                        }
                    }
                }
            }
        }
    }
}

// [[Rcpp::export]]
double killTest(Function ll, double x)
{
    NumericVector y = ll(x);
    return y[0];
}

// [[Rcpp::export]]
int rcppLinNest(NumericMatrix X, NumericVector y, NumericVector w, NumericVector adaweights, Function linkinv, Function devfun, int nrow, int ncol, int numGroup, IntegerVector rangeGroupInd, IntegerVector groupLen, NumericVector lambda, NumericMatrix beta, int innerIter, int outerIter, double thresh, double outerThresh, NumericVector eta, double gamma, IntegerVector betaIsZero, int reset)
{
    NumericVector prob(nrow);
    NumericVector nullBeta(ncol);
    int n = nrow;
    int p = ncol;
    NumericVector ldot(n);
    IntegerVector isActive(numGroup);
    IntegerVector useGroup(numGroup);
    IntegerVector tempIsActive(numGroup);
    int nlam = lambda.size();
    
    for (int step=0; step<nlam; step++)
    {
        for(int i=0; i<numGroup; i++)
        {
            isActive[i] = 0;
            useGroup[i] = 1;
        }
        
        //Copy the most recent betas into position
        if (step>0)
        {
            int l = step - 1;
            
            for (int i=0; i<ncol; i++)
            {
                beta[step*ncol + i] = beta[l*ncol + i];
            }
        }
        
        // outer most loop creating response etc...
        int outermostCounter = 0;
        double outermostCheck = 100000;
        NumericVector outerOldBeta(p);
        int groupChange = 1;
            
        while(groupChange == 1)
        {
            groupChange = 0;
            
            rcppLinSolver(X, y, w, adaweights, nrow, ncol, numGroup, beta, linkinv, devfun, rangeGroupInd, groupLen, lambda, step, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, isActive, useGroup, reset);
            
            while(outermostCounter < outerIter && outermostCheck > outerThresh)
            {
                outermostCounter ++;
                for(int i=0; i<p; i++)
                {
                    outerOldBeta[i] = beta[step*ncol + i];
                }
                
                for(int i=0; i<numGroup; i++)
                {
                    tempIsActive[i] = isActive[i];
                }
                
                rcppLinSolver(X, y, w, adaweights, nrow, ncol, numGroup, beta, linkinv, devfun, rangeGroupInd, groupLen, lambda, step, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, isActive, tempIsActive, reset);
                
                outermostCheck = 0;
                for(int i=0; i<p; i++)
                {
                    outermostCheck = outermostCheck + fabs(outerOldBeta[i] - beta[step*ncol + i]);
                }
                outermostCheck = outermostCheck / abs(sum(outerOldBeta));
            }
        }
    }
    return 1;
}
