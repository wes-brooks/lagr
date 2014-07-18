#include <iostream>
#include <math.h>

using namespace std;

extern "C" {
    
    //////////////////////////////////
    //This function returns the weighted mean of the difference between eta and y
    void linGradCalc(int *nrow, double *eta, double *y, double *w, double *ldot)
    {
        double sumw = 0;
        for(int i = 0; i < nrow[0]; i++)
        {
            sumw = sumw + w[i];
        }
        
        for(int i = 0; i < nrow[0]; i++)
        {
            ldot[i] = w[i]*(eta[i] - y[i])/sumw;  /* OR MAYBE NOT? */
        }
    }
    
    //Returns the weighted sum of squared residuals.
    double linNegLogLikelihoodCalc(int *nrow, double *eta, double *y, double *w)
    {
        double squareSum = 0;
        double sumw = 0;
        
        for(int i = 0; i < nrow[0]; i++)
        {
            squareSum = squareSum + w[i]*pow(eta[i] - y[i], 2)/2; 
            sumw = sumw + w[i];
        }
        
        return squareSum/sumw;   /* OR MAYBE NOT? */
    }
    
    
    void linSolver(double *X, double *y, double *w, int* index, double* adaweights, int *nrow, int *ncol, int *numGroup, double *beta, int *rangeGroupInd, int *groupLen, double *lambda, int *innerIter, double *thresh, double *ldot, double *nullBeta, double *gamma, double *eta, int* betaIsZero, int& groupChange, int* isActive, int* useGroup, double *step, int *reset)
    {
        double *theta = new double[ncol[0]];
        int startInd = 0;
        double zeroCheck = 0;
        double check = 0;
        int count = 0;
        double t = step[0];
        double diff = 1;
        double norm = 0;
        double uOp = 0;
        double Lnew = 0;
        double Lold = 0;
        double sqNormG = 0;
        double iProd = 0;
        double *etaNew = NULL;
        etaNew = new double[nrow[0]];
        double *etaNull = NULL;
        etaNull = new double[nrow[0]];
        
        for(int i = 0; i < numGroup[0]; i++)
        {
            if(useGroup[i] == 1)
            {
                startInd = rangeGroupInd[i];
                
                // Setting up null gradient calc to check if group is 0
                for(int k = 0; k < nrow[0]; k++)
                {
                    etaNull[k] = eta[k];
                    for(int j = startInd; j < rangeGroupInd[i] + groupLen[i]; j++)
                    {
                        etaNull[k] = etaNull[k] - X[k + nrow[0] * j] * beta[j]; 
                    }
                }
                
                // Calculating Null Gradient
                linGradCalc(nrow, etaNull, y, w, ldot);
                
                double *grad = NULL;
                grad = new double[groupLen[i]];
                
                // 
                for(int j = 0; j < groupLen[i]; j++)
                {
                    grad[j] = 0;
                    for(int k = 0; k < nrow[0]; k++)
                    {
                        grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
                    }
                }
                
                
                zeroCheck = 0;
                for(int j = 0; j < groupLen[i]; j++)
                {
                    zeroCheck = zeroCheck + pow(grad[j],2);
                }
                
                if(zeroCheck <= pow(adaweights[i],2)*pow(lambda[0],2)*groupLen[i])  //Or not?
                {
                    if(betaIsZero[i] == 0)
                    {
                        for(int k = 0; k < nrow[0]; k++)
                        {
                            for(int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++)
                            {
                                eta[k] = eta[k] - X[k + nrow[0] * j] * beta[j];
                            }
                        }
                    }
                    betaIsZero[i] = 1;
                    for(int j = 0; j < groupLen[i]; j++)
                    {
                        beta[j + rangeGroupInd[i]] = 0;
                    }
                }
                else
                {
                    if(isActive[i] == 0)
                    {
                        groupChange = 1;
                    }
                    isActive[i] = 1;
                    
                    for(int k = 0; k < ncol[0]; k++)
                    {
                        theta[k] = beta[k];
                    }
                    
                    betaIsZero[i] = 0;
                    double *z = NULL;
                    z = new double[groupLen[i]];
                    double *U = NULL;
                    U = new double[groupLen[i]];
                    double *G = NULL;
                    G = new double[groupLen[i]];
                    double *betaNew = NULL;
                    betaNew = new double[ncol[0]];
                    
                    count = 0;
                    check = 100000;
                    
                    while(count <= innerIter[0] && check > thresh[0])
                    {
                        count++;
                        
                        linGradCalc(nrow, eta, y, w, ldot);
                        
                        for(int j = 0; j < groupLen[i]; j++)
                        {		  
                            grad[j] = 0;
                            for(int k = 0; k < nrow[0]; k++)
                            {
                                grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
                            }
                        }
                        
                        diff = -1;
                        //	      t = 0.5;
                        Lold = linNegLogLikelihoodCalc(nrow, eta, y, w);
                        
                        // Back-tracking
                        while(diff < 0)
                        {
                            for(int j = 0; j < groupLen[i]; j++)
                            {
                                z[j] = beta[j + rangeGroupInd[i]] - t * grad[j];
                            }
                            
                            norm = 0;
                            for(int j = 0; j < groupLen[i]; j++)
                            {
                                norm = norm + pow(z[j],2);
                            }
                            norm = sqrt(norm);
                            
                            if(norm != 0){
                                uOp = (1 - adaweights[i]*lambda[0]*sqrt(double(groupLen[i]))*t/norm);   //Or not?
                            }
                            else{uOp = 0;}
                            
                            if(uOp < 0)
                            {
                                uOp = 0;
                            }
                            
                            for(int j = 0; j < groupLen[i]; j++)
                            {
                                U[j] = uOp*z[j];
                                G[j] = 1/t *(beta[j + rangeGroupInd[i]] - U[j]);
                            }
                            
                            // Setting up betaNew and etaNew in direction of Grad for descent step
                            for(int k = 0; k < nrow[0]; k++)
                            {
                                etaNew[k] = eta[k];
                                for(int j = 0; j < groupLen[i]; j++)
                                {
                                    etaNew[k] = etaNew[k] - t*G[j] * X[k + nrow[0]*(rangeGroupInd[i] + j)];
                                }
                            }
                            
                            Lnew = linNegLogLikelihoodCalc(nrow, etaNew, y, w);
                            
                            sqNormG = 0;
                            iProd = 0;
                            
                            for(int j = 0; j < groupLen[i]; j++)
                            {
                                sqNormG = sqNormG + pow(G[j],2);
                                iProd = iProd + grad[j] * G[j];
                            }
                            
                            diff = Lold - Lnew - t * iProd + t/2 * sqNormG;
                            
                            t = t * gamma[0];
                        }
                        t = t / gamma[0];
                        
                        check = 0;
                        
                        for(int j = 0; j < groupLen[i]; j++)
                        {
                            check = check + fabs(theta[j + rangeGroupInd[i]] - U[j]);
                            for(int k = 0; k < nrow[0]; k++)
                            {
                                eta[k] = eta[k] - X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
                            }
                            beta[j + rangeGroupInd[i]] = U[j] + count%reset[0]/(count%reset[0]+3) * (U[j] - theta[j + rangeGroupInd[i]]);
                            theta[j + rangeGroupInd[i]] = U[j];
                            
                            for(int k = 0; k < nrow[0]; k++)
                            {
                                eta[k] = eta[k] + X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
                            }
                        }
                    }
                    delete [] z;
                    delete [] U;
                    delete [] G;
                    delete [] betaNew;
                }
                delete [] grad;
            }
        }
        delete [] etaNew;
        delete [] etaNull;
        delete [] theta;
    }
    
    //' @export
    int linNest(double *X, double* y, double* w, int* index, double* adaweights, int *nrow, int *ncol, int *numGroup, int *rangeGroupInd, int *groupLen, double *lambda, double *beta, int *innerIter, int *outerIter, double *thresh, double *outerThresh, double *eta, double *gamma, int *betaIsZero, double *step, int *reset)
    {
        double* prob = NULL;
        prob = new double[nrow[0]];
        double* nullBeta = NULL;
        nullBeta = new double[ncol[0]];
        int n = nrow[0];
        int p = ncol[0];
        double *ldot = NULL;
        ldot = new double[n];
        int groupChange = 1;
        int* isActive = NULL;
        isActive = new int[numGroup[0]];
        int* useGroup = NULL;
        useGroup = new int[numGroup[0]];
        int* tempIsActive = NULL;
        tempIsActive = new int[numGroup[0]];
        
        for(int i = 0; i < numGroup[0]; i++)
        {
            isActive[i] = 0;
            useGroup[i] = 1;
        }
        
        // outer most loop creating response etc...
        int outermostCounter = 0;
        double outermostCheck = 100000;
        double* outerOldBeta = NULL;
        outerOldBeta = new double[p];
        
        while(groupChange == 1)
        {
            groupChange = 0;
            
            linSolver(X, y, w, index, adaweights, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, isActive, useGroup, step, reset);
            
            while(outermostCounter < outerIter[0] && outermostCheck > outerThresh[0])
            {
                outermostCounter ++;
                for(int i = 0; i < p; i++)
                {
                    outerOldBeta[i] = beta[i];
                }
                
                for(int i = 0; i < numGroup[0]; i++)
                {
                    tempIsActive[i] = isActive[i];
                }
                
                linSolver(X, y, w, index, adaweights, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, isActive, tempIsActive, step, reset);
                
                outermostCheck = 0;
                for(int i = 0; i < p; i++)
                {
                    outermostCheck = outermostCheck + fabs(outerOldBeta[i] - beta[i]);
                }
            }
        }
        
        delete [] nullBeta;
        delete [] outerOldBeta;
        delete [] ldot;
        delete [] isActive;
        delete [] useGroup;
        delete [] tempIsActive;
        delete [] prob;
        
        return 1;
    }
    
    
    
    //////////////////////////////////
    //Logistic regression:
    
    void logitGradCalc(int *nrow, double *prob, int *y, double *ldot)
    {
        for(int i = 0; i < nrow[0]; i++)
        {
            ldot[i] = (-y[i] + prob[i])/nrow[0]; /* OR MAYBE NOT? */
        }
    }
    
    void pCalc(int *nrow, double *eta, double *prob)
    {
        for(int i = 0; i < nrow[0]; i++)
        {
            prob[i] = exp(eta[i]) / (1 + exp(eta[i]));
        }
    }
    
    
    double logitNegLogLikelihoodCalc(int *nrow, double *prob, int *y)
    {
        double logLik = 0;
        
        for(int i = 0; i < nrow[0]; i++)
        {
            logLik = logLik + y[i] * log(prob[i]) + (1 - y[i]) * log(1 - prob[i]);
        }
        
        return -logLik/nrow[0];  /* OR MAYBE NOT? */
    }
    
    void betaZeroSolve(int *nrow, double *betaZero, double *eta, double *prob, double *thresh, int *innerIter, int *y)
    {
        double diff = 10;
        double num = 0;
        double denom = 0;
        int count = 0;
        
        while(pow(diff,2) > pow(thresh[0],2) && count < innerIter[0])
        {
            pCalc(nrow, eta, prob);
            diff = 0;
            
            for(int i = 0; i < nrow[0]; i++)
            {
                num = num + y[i] - prob[i];
                denom = denom + prob[i] * (1 - prob[i]);
            }
            diff = num/denom;
            betaZero[0] = betaZero[0] - diff;
            
            for(int i = 0; i < nrow[0]; i++)
            {
                eta[i] = eta[i] + diff;
            }
        }
    }
    
    
    void logitSolver(double *X, int *y, int* index, int *nrow, int *ncol, int *numGroup, double *beta, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, int *innerIter, double *thresh, double *ldot, double *nullBeta, double *gamma, double *eta, int* betaIsZero, int& groupChange, int* isActive, int* useGroup, double *prob, double *betaZero, double *step)
    {
        double *theta = new double[ncol[0]];
        double *thetaNew = new double[ncol[0]];
        int startInd = 0;
        double zeroCheck = 0;
        double check = 0;
        int count = 0;
        double t = step[0];
        double diff = 1;
        double norm = 0;
        double uOp = 0;
        double Lnew = 0;
        double Lold = 0;
        double sqNormG = 0;
        double iProd = 0;
        double *etaNew = NULL;
        etaNew = new double[nrow[0]];
        double *etaNull = NULL;
        etaNull = new double[nrow[0]];
        
        for(int i = 0; i < numGroup[0]; i++)
        {
            if(useGroup[i] == 1)
            {
                startInd = rangeGroupInd[i];
                
                // Setting up null gradient calc to check if group is 0
                for(int k = 0; k < nrow[0]; k++)
                {
                    etaNull[k] = eta[k];
                    for(int j = startInd; j < rangeGroupInd[i] + groupLen[i]; j++)
                    {
                        etaNull[k] = etaNull[k] - X[k + nrow[0] * j] * beta[j]; 
                    }
                }
                
                // Calculating Null Gradient
                pCalc(nrow, etaNull, prob);
                logitGradCalc(nrow, prob, y, ldot);
                
                double *grad = NULL;
                grad = new double[groupLen[i]];
                
                for(int j = 0; j < groupLen[i]; j++)
                {
                    grad[j] = 0;
                    for(int k = 0; k < nrow[0]; k++)
                    {
                        grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
                    }
                    if(grad[j] < lambda1[0] && grad[j] > -lambda1[0])
                    {
                        grad[j] = 0;
                    }
                    if(grad[j] > lambda1[0])
                    {
                        grad[j] = grad[j] - lambda1[0];
                    }
                    if(grad[j] < -lambda1[0])
                    {
                        grad[j] = grad[j] + lambda1[0];
                    }
                    if(pow(grad[j],2) == pow(lambda1[0],2))
                    {
                        grad[j] = 0;
                    }
                }
                
                zeroCheck = 0;
                for(int j = 0; j < groupLen[i]; j++)
                {
                    zeroCheck = zeroCheck + pow(grad[j],2);
                }
                
                if(zeroCheck <= pow(lambda2[0],2) * groupLen[i])   //Or not?
                {
                    if(betaIsZero[i] == 0)
                    {
                        for(int k = 0; k < nrow[0]; k++)
                        {
                            for(int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++)
                            {
                                eta[k] = eta[k] - X[k + nrow[0] * j] * beta[j];
                            }
                        }
                    }
                    betaIsZero[i] = 1;
                    for(int j = 0; j < groupLen[i]; j++)
                    {
                        beta[j + rangeGroupInd[i]] = 0;
                    }
                }
                else
                {
                    if(isActive[i] == 0)
                    {
                        groupChange = 1;
                    }
                    isActive[i] = 1;
                    
                    for(int k = 0; k < ncol[0]; k++)
                    {
                        theta[k] = beta[k];
                    }
                    
                    betaIsZero[i] = 0;
                    double *z = NULL;
                    z = new double[groupLen[i]];
                    double *U = NULL;
                    U = new double[groupLen[i]];
                    double *G = NULL;
                    G = new double[groupLen[i]];
                    double *betaNew = NULL;
                    betaNew = new double[ncol[0]];
                    
                    count = 0;
                    check = 1000000;
                    
                    while(count <= innerIter[0] && check > thresh[0])
                    {
                        count++;
                        
                        pCalc(nrow, eta, prob);
                        logitGradCalc(nrow, prob, y ,ldot);
                        
                        for(int j = 0; j < groupLen[i]; j++)
                        {		  
                            grad[j] = 0;
                            for(int k = 0; k < nrow[0]; k++)
                            {
                                grad[j] = grad[j] + X[k + nrow[0] * (j + rangeGroupInd[i])] * ldot[k];
                            }
                        }
                        
                        diff = -1;
                        //	      t = 0.5;
                        pCalc(nrow, eta, prob);
                        Lold = logitNegLogLikelihoodCalc(nrow, prob, y);
                        
                        // Back-tracking
                        while(diff < 0)
                        {
                            for(int j = 0; j < groupLen[i]; j++)
                            {
                                z[j] = beta[j + rangeGroupInd[i]] - t * grad[j];
                                if(z[j] < lambda1[0] * t && z[j] > -lambda1[0] * t)
                                {
                                    z[j] = 0;
                                }
                                if(z[j] > lambda1[0] * t)
                                {
                                    z[j] = z[j] - lambda1[0] * t;
                                }
                                if(z[j] < -lambda1[0] * t)
                                {
                                    z[j] = z[j] + lambda1[0] * t;
                                }
                            }
                            
                            norm = 0;
                            for(int j = 0; j < groupLen[i]; j++)
                            {
                                norm = norm + pow(z[j],2);
                            }
                            norm = sqrt(norm);
                            
                            if(norm != 0)
                            {
                                uOp = (1 - lambda2[0]*sqrt(double(groupLen[i]))*t/norm);   //Or not?
                            }
                            else{uOp = 0;}
                            
                            if(uOp < 0)
                            {
                                uOp = 0;
                            }
                            
                            for(int j = 0; j < groupLen[i]; j++)
                            {
                                U[j] = uOp*z[j];
                                G[j] = 1/t *(beta[j + rangeGroupInd[i]] - U[j]);
                            }
                            
                            // Setting up betaNew and etaNew in direction of Grad for descent step
                            for(int j = 0; j < rangeGroupInd[i]; j++)
                            {
                                thetaNew[j] = beta[j];
                            }
                            for(int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++)
                            {
                                thetaNew[j] = beta[j] - t * G[j];
                            }
                            
                            for(int j = rangeGroupInd[i] + groupLen[i]; j < ncol[0]; j++)
                            {
                                thetaNew[j] = beta[j];
                            }
                            for(int k = 0; k < nrow[0]; k++)
                            {
                                etaNew[k] = eta[k];
                                for(int j = 0; j < groupLen[i]; j++)
                                {
                                    etaNew[k] = etaNew[k] - t*G[j] * X[k + nrow[0]*(rangeGroupInd[i] + j)];
                                }
                            }
                            
                            pCalc(nrow, etaNew, prob);
                            Lnew = logitNegLogLikelihoodCalc(nrow, prob, y);
                            
                            sqNormG = 0;
                            iProd = 0;
                            
                            for(int j = 0; j < groupLen[i]; j++)
                            {
                                sqNormG = sqNormG + pow(G[j],2);
                                iProd = iProd + grad[j] * G[j];
                            }
                            
                            diff = Lold - Lnew - t * iProd + t/2 * sqNormG;
                            
                            t = t * gamma[0];
                        }
                        t = t / gamma[0];
                        
                        check = 0;
                        
                        for(int j = 0; j < groupLen[i]; j++)
                        {
                            check = check + fabs(theta[j + rangeGroupInd[i]] - U[j]);
                            for(int k = 0; k < nrow[0]; k++)
                            {
                                eta[k] = eta[k] - X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
                            }
                            beta[j + rangeGroupInd[i]] = U[j] + count/(count+3) * (U[j] - theta[j + rangeGroupInd[i]]);
                            theta[j + rangeGroupInd[i]] = U[j];
                            
                            for(int k = 0; k < nrow[0]; k++)
                            {
                                eta[k] = eta[k] + X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
                            }
                        }
                    }
                    delete [] z;
                    delete [] U;
                    delete [] G;
                    delete [] betaNew;
                }
                delete [] grad;
            }
        }
        betaZeroSolve(nrow, betaZero, eta, prob, thresh, innerIter, y);
        
        delete [] etaNew;
        delete [] etaNull;
        delete [] theta;
        delete [] thetaNew;
    }
    
    
    
    int logitNest(double *X, int* y, int* index, int *nrow, int *ncol, int *numGroup, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, double *beta, int *innerIter, int *outerIter, double *thresh, double *outerThresh, double *eta, double *gamma, int *betaIsZero, double* betaZero, double *step)
    {
        double oldBetaZero = betaZero[0];
        double* prob = NULL;
        prob = new double[nrow[0]];
        double* nullBeta = NULL;
        nullBeta = new double[ncol[0]];
        int n = nrow[0];
        int p = ncol[0];
        double *ldot = NULL;
        ldot = new double[n];
        int groupChange = 1;
        int* isActive = NULL;
        isActive = new int[numGroup[0]];
        int* useGroup = NULL;
        useGroup = new int[numGroup[0]];
        int* tempIsActive = NULL;
        tempIsActive = new int[numGroup[0]];
        
        for(int i = 0; i < numGroup[0]; i++)
        {
            isActive[i] = 0;
            useGroup[i] = 1;
        }
        
        // outer most loop creating response etc...
        int outermostCounter = 0;
        double outermostCheck = 1000000;
        double* outerOldBeta = NULL;
        outerOldBeta = new double[p];
        
        while(groupChange == 1)
        {
            groupChange = 0;
            
            logitSolver(X, y, index, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, isActive, useGroup, prob, betaZero, step);
            
            while(outermostCounter < outerIter[0] && outermostCheck > outerThresh[0])
            {
                outermostCounter ++;
                for(int i = 0; i < p; i++)
                {
                    outerOldBeta[i] = beta[i];
                }
                oldBetaZero = betaZero[0];
                
                for(int i = 0; i < numGroup[0]; i++)
                {
                    tempIsActive[i] = isActive[i];
                }
                
                
                logitSolver(X, y, index, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, ldot, nullBeta, gamma, eta, betaIsZero, groupChange, tempIsActive, isActive, prob, betaZero, step);
                
                outermostCheck = 0;
                for(int i = 0; i < p; i++)
                {
                    outermostCheck = outermostCheck + fabs(outerOldBeta[i] - beta[i]);
                }
                outermostCheck = outermostCheck + fabs(oldBetaZero - betaZero[0]);
            }
        }
        
        delete [] nullBeta;
        delete [] outerOldBeta;
        delete [] ldot;
        delete [] isActive;
        delete [] useGroup;
        delete [] tempIsActive;
        delete [] prob;
        
        return 1;
    } //loigitNest
    
} //extern "C"
