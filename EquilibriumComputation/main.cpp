//
//  main.cpp
//  UniversityFP
//
//  Created by Rodrigo Azuero on 4/26/16.
//  Copyright (c) 2016 Rodrigo Azuero Melo. All rights reserved.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
//#include <RcppArmadillo.h>
//#include <RcppEigen.h>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <unistd.h>
#include <nlopt.hpp>
using std::vector;
using namespace std;




//Define template for random number generator
template<class T>
double gen_normal_3(T &generator)
{
    return generator();
}



//======================================
//V. Fiscal policy tools
//======================================


//V.a. Lump sum tax for households
// [[Rcpp::export]]
double bfinal(double binitial, double tax){
    return(binitial*(1-tax));
}

//V.b Subsidy for interest rate
// [[Rcpp::export]]
double rsub(double ttheta,
            double tthetathreshold,
            double bb,
            double bbthreshold,
            double subsidy){
    
    if (ttheta<tthetathreshold ||
        bb>bbthreshold){
        return(0);
    }
    else{
        return(subsidy);
    }
}

//V.c Borrowing constraint for high ability students
// [[Rcpp::export]]
double AbarC(double ttheta,
             double tthetathreshold,
             double Abar){
    
    if (ttheta<tthetathreshold){
        return(0);
    }
    else{
        return(Abar);
    }
}


double Phis(double bbeta, double r, double T, double S, double ssigma){
  double term = pow(bbeta / pow(1+r,ssigma-1), 1/ssigma);
  return( (1 - pow(term, T-S+1))/(1 - term));
}

double Phi0(double bbeta, double r, double T, double S, double ssigma){
  double term = pow(bbeta / pow(1+r,ssigma-1), 1/ssigma);
  return( (1 - pow(term, S))/(1 - term));
}

double Phiro(double bbeta, double r, double T, double S, double ssigma){
  double term = 1 / (1+r);
  return( (1 - pow(term, T-S+1))/(1 - term));
}

double Phiry(double bbeta, double r, double T, double S, double ssigma){
  double term = 1 / (1+r);
  return( (1 - pow(term, S))/(1 - term));
}


//-----------------------------------------------
// I. First section of file. Computing utilities
//    of different alternatives.
//=----------------------------------------------


// I. a Utility of study - unconstrained: This will be used to find
//      utility of both studying in good and bad university

// [[Rcpp::export]]
double UNCSTUDY( double b,
                double ttheta,
                double ssigma,
                double bbeta,
                double w,
                double P,
                double ggamma,
                double r,
                double Abar,
                double tthetaadmission,
                double tthetaRthreshold,
                double bRthreshold,
                double Rsubsidy,
                double tax, double T, double S){
    // Relevant Abar
    Abar = AbarC(ttheta, tthetaRthreshold, Abar);
    
    //Computing subsidized interest rate
    
    double subR=rsub(ttheta, tthetaRthreshold, b, bRthreshold, Rsubsidy);
    
    //Computing thefinal income
    b=bfinal(b, tax);
    
    
    
    //If person's ttheta is below the threshold of admission, utility of -10000
    if (ttheta<tthetaadmission){
        return(-10000);
    }
    else{
      
      double phi0 = Phi0(bbeta, r, T, S, ssigma);
      double phis = Phis(bbeta, r, T, S, ssigma);
      double phiro = Phiro(bbeta, r, T, S, ssigma);
      double phiry = Phiry(bbeta, r, T, S, ssigma);
      
      double term = w*ttheta*phiro/phis + (b*pow(1+r, S))/phis - P*phiry*pow(1+r, S)/phis;
      double disc = pow(bbeta*(1+r), (-S/ssigma));
      double deno = 1 + disc*(phi0*pow(1+r, S)/phis);
      
      double cons0 = disc*term/deno;
      double cons1 = pow(bbeta*(1+r), (S/ssigma))*cons0;
      
      double bbetatilda = pow(bbeta, S)*phis/phi0;
      
      double utility=pow(cons0,1-ssigma)/(1-ssigma) + bbetatilda*pow(cons1,1-ssigma)/(1-ssigma);
      
      //Finally, we need to see if the borrowing constraint is not satisfied:
      double deuda=cons0+P-b;
      if (deuda>Abar){
          utility=pow(-10.0,5.0);
      }
      return(utility);
    }
}


//I.b Utility if study constrained

// [[Rcpp::export]]
double CSTUDY(double b,
              double ttheta,
              double ssigma,
              double bbeta,
              double w,
              double P,
              double ggamma,
              double r,
              double Abar,
              double tthetaadmission,
              double tthetaRthreshold,
              double bRthreshold,
              double Rsubsidy,
              double tax, double T, double S){
    
    // Relevant Abar
    Abar = AbarC(ttheta, tthetaRthreshold, Abar);
    
    //Subsity to interest rate
    double subR=rsub(ttheta, tthetaRthreshold, b, bRthreshold, Rsubsidy);
    
    //Final income after taxes
    b=bfinal(b, tax);
    
    
    //If person's ttheta is below the threshold of admission, utility of -10000
    if (ttheta<tthetaadmission){
        return(-10000);
    }
    else{
        
      double phi0 = Phi0(bbeta, r, T, S, ssigma);
      double phis = Phis(bbeta, r, T, S, ssigma);
      double phiro = Phiro(bbeta, r, T, S, ssigma);
      double phiry = Phiry(bbeta, r, T, S, ssigma);
      
      //Consumption in t=0
      double cons0 = (1/phi0)*(b+(Abar/pow(1+r, S)) - P*phiry);
      
      
      //Consumption in t=1
      double cons1=(1/phis)*(w*ttheta*phiro-Abar);
      
      double bbetatilda = pow(bbeta, S)*phis/phi0;
      
      //Computing utility:
      double utility=0;
      
      
      //Consumption negative,bad utility
      if (cons0 <0 || cons1<0){
          utility=pow(-10.0,5.0);
      }
      
      //Consumption positive, compute usual utility
      else{
          utility=pow(cons0,1-ssigma)/(1-ssigma) + bbetatilda*pow(cons1,1-ssigma)/(1-ssigma);
      }
      return(utility);
    }
}

//I.c Utility of unconstrained if not study

// [[Rcpp::export]]
double UNCNOTSTUDY(double b,
                   double ttheta,
                   double ssigma,
                   double bbeta,
                   double w,
                   double ggamma,
                   double r,
                   double Abar,
                   double tax, double T, double S){
    
    // Relevant Abar
    Abar = 0;
    
    //Computing final income
    b=bfinal(b,tax);
    
    //Consumption in t=0
    double phi0 = Phi0(bbeta, r, T, S, ssigma);
    double phis = Phis(bbeta, r, T, S, ssigma);
    double phiro = Phiro(bbeta, r, T, S, ssigma);
    double phiry = Phiry(bbeta, r, T, S, ssigma);
    
    double term = w*ttheta*((phiro + pow(1+r,S)*phiry)/phis) + (b*pow(1+r, S))/phis;
    double disc = pow(bbeta*(1+r), (-S/ssigma));
    double deno = 1 + disc*(phi0*pow(1+r, S)/phis);
    
    double cons0 = disc*term/deno;
    double cons1 = pow(bbeta*(1+r), (S/ssigma))*cons0;

    //Computing utility if studies
    
    double bbetatilda = pow(bbeta, S)*phis/phi0;
    
    double utility=pow(cons0,1-ssigma)/(1-ssigma) + bbetatilda*pow(cons1,1-ssigma)/(1-ssigma);
    
    //Finally, check if borrowing constraint is satisfied
    double deuda=cons0-b-w*ttheta;
    if (deuda>Abar){
        utility=pow(-10.0,5.0);
    }
    
    return(utility);
}


//I.d Utility if constrained and not study

// [[Rcpp::export]]
double CNOTSTUDY(double b,
                 double ttheta,
                 double ssigma,
                 double bbeta,
                 double w,
                 double ggamma,
                 double r,
                 double Abar,
                 double tax, double T, double S){
    
    // Relevant Abar
    Abar = 0;
    
    //Final income
    b=bfinal(b, tax);
    
    double phi0 = Phi0(bbeta, r, T, S, ssigma);
    double phis = Phis(bbeta, r, T, S, ssigma);
    double phiro = Phiro(bbeta, r, T, S, ssigma);
    double phiry = Phiry(bbeta, r, T, S, ssigma);
    
    //Consumption in t=0
    double cons0 = (w*ttheta*phiry + b + Abar/pow(1+r, S)) * (1/phi0);
    
    
    //Consumption in t=1
    double cons1=(w*ttheta*phiro - Abar)*(1/phis);
    
    double bbetatilda = pow(bbeta, S)*phis/phi0;
    
    //Computing utility:
    double utility=0;
    
    
    //Consumption negative,bad utility
    if (cons0 <0 || cons1<0){
        utility=pow(-10.0,5.0);
    }
    
    
    //Consumption positive, compute usual utility
    else{
        utility=pow(cons0,1-ssigma)/(1-ssigma) + bbetatilda*pow(cons1,1-ssigma)/(1-ssigma);
    }
    
    return(utility);
}

//For each possibility of the household I will put the level of debt that will be subsidized by the government. I will
//only have to do it for CSTUDY y UNCSTUDY

//Debt if constrained and study

// [[Rcpp::export]]
double DEBTCSTUDY(double b,
                  double ttheta,
                  double ssigma,
                  double bbeta,
                  double w,
                  double P,
                  double ggamma,
                  double r,
                  double Abar,
                  double tthetaadmission,
                  double tthetaRthreshold,
                  double bRthreshold,
                  double Rsubsidy,
                  double tax, double T, double S){
    
    // Relevant Abar
    Abar = AbarC(ttheta, tthetaRthreshold, Abar);
    
    //Subsity to interest rate
    double subR=rsub(ttheta, tthetaRthreshold, b, bRthreshold, Rsubsidy);
    
    //Final income after taxes
    b=bfinal(b, tax);
    
    
    //Debt level:
    double debt=0;
    
    //If person's ttheta is below the threshold of admission, utility of -10000
    if (ttheta<tthetaadmission){
        return(-10000);
    }
    else{
        
        //Consumption in t=0
        double cons0=b-P+Abar;
        
        
        //Consumption in t=1
        double cons1=w*ttheta-(1+r)*Abar+r*subR*min(Abar,P);
        debt=r*subR*min(Abar,P);
        
        //Computing utility:
        double utility=0;
        
        
        //Consumption negative,bad utility
        if (cons0 <0 || cons1<0){
            utility=pow(-10.0,5.0);
        }
        
        //Consumption positive, compute usual utility
        else{
            utility=pow(cons0,1-ssigma)/(1-ssigma)-
            ggamma*(1/(1+ttheta)-0.5)+bbeta*pow(cons1,1-ssigma)/(1-ssigma);
        }
        return(debt);
    }
}


// [[Rcpp::export]]
double DEBTUNCSTUDY( double b,
                    double ttheta,
                    double ssigma,
                    double bbeta,
                    double w,
                    double P,
                    double ggamma,
                    double r,
                    double Abar,
                    double tthetaadmission,
                    double tthetaRthreshold,
                    double bRthreshold,
                    double Rsubsidy,
                    double tax){
    // Relevant Abar
    Abar = AbarC(ttheta, tthetaRthreshold, Abar);
    
    //Computing subsidized interest rate
    
    double subR=rsub(ttheta, tthetaRthreshold, b, bRthreshold, Rsubsidy);
    
    //Computing thefinal income
    b=bfinal(b, tax);
    
    //Getting debt
    double debt=0;
    //If person's ttheta is below the threshold of admission, utility of -10000
    if (ttheta<tthetaadmission){
        return(-10000);
    }
    else{
        //Consumption in t=0
        double num0=bbeta*(1+r);
        num0=pow(num0,-1/ssigma);
        
        double num1=w*ttheta+b*(1+r)-P*(1+r)
        +min(P,Abar)*r*subR;
        
        //Getting the debt level
        debt=min(P,Abar)*r*subR;
        
        
        double den0=bbeta*(1+r);
        den0=pow(den0,-1/ssigma);
        den0=(1+r)*den0;
        den0=1+den0;
        
        double cons0=num0*num1/den0;
        
        
        //Consumption in t=1
        double cons1=num1/den0;
        
        //Computing utility if studies
        
        double utility=pow(cons0,1-ssigma)/(1-ssigma)-
        ggamma*(1/(1+ttheta)-0.5)+bbeta*pow(cons1,1-ssigma)/(1-ssigma);
        
        
        //Finally, we need to see if the borrowing constraint is not satisfied:
        double deuda=cons0+P-b;
        if (deuda>Abar){
            utility=pow(-10.0,5.0);
        }
        return(debt);
    }
}




//I. e. Function comparing utilities and decision rules of
//      different alternatives. Choose the one with highest utility.
//      decision == 1 if individual studies

// [[Rcpp::export]]
vector<double> decision( double b,
                        double ttheta,
                        double ssigma,
                        double bbeta,
                        double w,
                        double wl,
                        double wh,
                        double Pl,
                        double Ph,
                        double tthetaadmissionl,
                        double tthetaadmissionh,
                        double ggamma,
                        double r,
                        double Abar,
                        double tthetaRthreshold,
                        double bRthreshold,
                        double Rsubsidy,
                        double tax, double T, double S){
    //The vector decision will store the decision and the amount of debt that is subject for govt. subsidy for each person.
    vector<double> decision;
    decision.resize(2);
    
    double uncstudy_h = UNCSTUDY(b, ttheta, ssigma, bbeta, wh, Ph, ggamma, r, Abar, tthetaadmissionh, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);

    double uncstudy_l = UNCSTUDY(b, ttheta, ssigma, bbeta, wl, Pl, ggamma, r, Abar, tthetaadmissionl, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
    
    double cstudy_h = CSTUDY(b, ttheta, ssigma, bbeta, wh, Ph, ggamma, r, Abar, tthetaadmissionh, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
    
    double debtcstudy_h=DEBTCSTUDY(b, ttheta, ssigma, bbeta, wh, Ph, ggamma, r, Abar, tthetaadmissionh, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
    
    double cstudy_l = CSTUDY(b, ttheta, ssigma, bbeta, wl, Pl, ggamma, r, Abar, tthetaadmissionl, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
    
    double debtcstudy_l=DEBTCSTUDY(b, ttheta, ssigma, bbeta, wh, Pl, ggamma, r, Abar, tthetaadmissionl, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
    
    double uncnotstudy = UNCNOTSTUDY(b, ttheta, ssigma, bbeta, w, ggamma, r, Abar, tax, T, S);
    double cnotstudy = CNOTSTUDY(b, ttheta, ssigma, bbeta, w, ggamma, r, Abar, tax, T, S);
    
    if(uncstudy_h > uncstudy_l && uncstudy_h > cstudy_l && uncstudy_h > uncnotstudy && uncstudy_h > cnotstudy){
        decision[0] = 1;
    }
    else if (cstudy_h > uncstudy_l && cstudy_h > cstudy_l && cstudy_h > uncnotstudy && cstudy_h > cnotstudy){
        decision[0] = 1;
    }
    else if(uncstudy_l > uncstudy_h && uncstudy_l > cstudy_h && uncstudy_l > uncnotstudy && uncstudy_l > cnotstudy){
        decision[1] = 1;
    }
    else if (cstudy_l > uncstudy_h && cstudy_l > cstudy_h && cstudy_l > uncnotstudy && cstudy_l > cnotstudy){
        decision[1] = 1;
    }
    return(decision);
}




//==========================================================
//II. b. Computing the mean ttheta of students going to univ
//==========================================================
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
vector<double> aaverages(double zh,
                         double zl,
                         vector<double> bbgrid,
                         vector<double> ttgrid,
                         double ssigma,
                         double bbeta,
                         double w,
                         double Pl,
                         double Ph,
                         double ggamma,
                         double Abar,
                         double r,
                         double tthetaadmissionl,
                         double tthetaadmissionh,
                         double tthetaRthreshold,
                         double bRthreshold,
                         double Rsubsidy,
                         double tax, double T, double S){
    
    // Returns a vector with <thetaaverage1, bbaverage1, numberstudents1, thetaaverage2, bbaverage2, numberstudents2>
    vector<double> result;
    result.resize(6);
    
    //0. bbgrid: grid of bequests
    //   ttgrid: grid of ttheta
    
    //0.0. Getting the size of matrix of distributions
    int bbnumber=bbgrid.size(); //Size along b dim
    int ttnumber=ttgrid.size(); //Size along ttheta dimension.
    
    //0.1 Generating matrix
    vector<vector<double> > MATGRIDh;
    vector<vector<double> > MATGRIDl;
    MATGRIDh.resize(bbnumber);
    MATGRIDl.resize(bbnumber);
    for (int n=0; n<bbnumber; n++){
        MATGRIDh[n].resize(ttnumber);
        MATGRIDl[n].resize(ttnumber);
    }
    
    //0.2. Once we have the size we can loop for every bb and every ttheta, the optimal decission for each person given prices and else;
    //Need to generate two counters.
    //One will be for average and the other for
    //the total number of students
    
    double ttstudentsh = 0;
    double aaverageh = 0;
    double bbaverageh = 0;
    
    double ttstudentsl = 0;
    double aaveragel = 0;
    double bbaveragel = 0;
    
    //For each student, identifying if it is
    //worth or not studying
    
    vector<double> res;
    res.resize(3);
#pragma omp parallel for shared(MATGRIDh, MATGRIDl) private(res)
    for (int bb=0; bb<bbnumber; bb=bb+1){
        for (int tt=0; tt<ttnumber; tt=tt+1){
            res = decision(bbgrid[bb], ttgrid[tt], ssigma, bbeta, w, w*(1+zl), w*(1+zh), Pl, Ph, tthetaadmissionl, tthetaadmissionh, ggamma, r, Abar, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
            MATGRIDh[bb][tt] = res[0];
            MATGRIDl[bb][tt] = res[1];
        }
    }
    
    for (int bb=0; bb<bbnumber; bb=bb+1){
        for (int tt=0; tt<ttnumber; tt=tt+1){
            ttstudentsh = ttstudentsh + MATGRIDh[bb][tt];
            aaverageh   = aaverageh + ttgrid[tt]*MATGRIDh[bb][tt];
            bbaverageh  = bbaverageh + bbgrid[bb]*MATGRIDh[bb][tt];
            
            ttstudentsl = ttstudentsl + MATGRIDl[bb][tt];
            aaveragel   = aaveragel + ttgrid[tt]*MATGRIDl[bb][tt];
            bbaveragel  = bbaveragel + bbgrid[bb]*MATGRIDl[bb][tt];
        }
    }
    
    
    if (ttstudentsh<1){
        result[0] = -1000;
        result[1] = -1000;
        result[2] = -1000;
    }
    else{
        aaverageh  = aaverageh / ttstudentsh;
        bbaverageh = bbaverageh / ttstudentsh;
        
        result[0] = aaverageh;
        result[1] = bbaverageh;
        result[2] = ttstudentsh;
    }
    if (ttstudentsl<1){
        result[3] = -1000;
        result[4] = -1000;
        result[5] = -1000;
    }
    else{
        aaveragel  = aaveragel / ttstudentsl;
        bbaveragel = bbaveragel / ttstudentsl;
        
        result[3] = aaveragel;
        result[4] = bbaveragel;
        result[5] = ttstudentsl;
    }    //Compute the z being offered
    return(result);
}



//Function of variable cost for universities
double VariableCHigh(double Nh){
    double result=Nh*0.003+0.00002*pow(Nh,2);
    return(result);
}

double VariableCLow(double Nl){
    double result=Nl*0.003+0.00002*pow(Nl,2);
    return(result);
}


//Zvector spit out
// [[Rcpp::export]]
vector<double> Zs(double aalpha1, double aalpha2,
                  double RMCH, double RMCL,
                  double zHof,
                  double zLof,
                  vector<double> bbgrid,
                  vector<double> ttgrid,
                  double ssigma,
                  double bbeta,
                  double w,
                  double pL,
                  double pH,
                  double ggamma,
                  double Abar,
                  double r,
                  double tthetaLthres,
                  double tthetaHthres,
                  double tthetaRthreshold,
                  double bRthreshold,
                  double Rsubsidy,
                  double tax,
                  double ZHPAR,
                  double ZLPAR, double T, double S){
    
    //0. Block of explanation of parameters:
    //0.1 pH-> Price of high education
    //0.2 tthetaHtrhes-> threshold for attendance in sskill for h
    //0.3 RMCH-> R-C in h
    //0.4 aalpha1-> exponent of mean ttheta in Z prod f
    //0.5 aalpha2-> exponent of I in Z prod fun.
    //0.6 zHof-> zH initial offer.
    
    vector<double> x;
    x.resize(6);
    
    x = aaverages(zHof, zLof, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
    
    //1. Block of computing intermediate functions
    
    double tthetaHbar = x[0]; //Funciton to compute tthetabarra of h
    double tthetaLbar = x[3]; //Function to compute tthetabarra of l
    
    double Nh = x[2]; //Function to compute number of students in h
    double Nl = x[5]; //Function to compute the number of students in l
    
    double VH = VariableCHigh(Nh); //Function to compute variable cost of h
    double VL = VariableCLow(Nl); //Function to compute variable cost of l
    
    double IH = (pH*Nh+(RMCH)-VH)/Nh;
    double IL = (pL*Nl+(RMCL)-VL)/Nl;
    
    if( Nl <= 0 || IL < 0){
        IL = 0;
    }
    if(Nh <=0 || IH < 0){
        IH = 0;
    }
    if(tthetaLbar<=0){
        tthetaLbar=0;
    }
    if(tthetaHbar<=0){
        tthetaHbar=0;
    }
    //2. Block of computing final outputs
    double zH=ZHPAR*pow(tthetaHbar,aalpha1)*pow(IH,aalpha2);
    double zL=ZLPAR*pow(tthetaLbar,aalpha1)*pow(IL,aalpha2);
    if(1==2){
        cout << " -----zs----" << endl;
        cout << Nh << " Nh " << endl;
        cout << Nl << " Nl " << endl;
        cout << VL << " VL " << endl;
        cout << VH << " VH " << endl;
         cout << RMCH << " RMCH " << endl;
        cout << RMCL << " RMCL " << endl;
        cout << pH << " pH " << endl;
        cout << pL << " pL " << endl;
        cout << pH*Nh+(RMCH)-VH << " pH*Nh+(RMCH)-VH  " << endl;
        cout << pL*Nl+(RMCL)-VL  << " pL*NL+(RMCl)-VL  " << endl;
        cout << (pH*Nh+(RMCH)-VH)/Nh << " pH*Nh+(RMCH)-VH  " << endl;
        cout << (pL*Nl+(RMCL)-VL)/Nl  << " pL*NL+(RMCl)-VL  " << endl;
        cout << tthetaHbar << " tthetaHbar " << endl;
        cout << tthetaLbar << " tthetaLvar " << endl;
        cout << IH << " IH " << endl;
        cout << IL << " IL " << endl;
        cout << zHof << " zHof " << endl;
        cout <<  zLof << " zlof " << endl;
        cout << zH  << " zH " << endl;
        cout << zL  << " zL " << endl;
    }
    //3. Block of loading the final return element
    
    //Return element
    vector<double> ZSPIT(2);
    ZSPIT[0]=zH;
    ZSPIT[1]=zL;
    return(ZSPIT);
}

//Funciton to compute the error between the proposed Z and final Z
// [[Rcpp::export]]
double errZ(double aalpha1, double aalpha2,
            double RMCH, double RMCL,
            double zHof,
            double zLof,
            vector<double> bbgrid,
            vector<double> ttgrid,
            double ssigma,
            double bbeta,
            double w,
            double pL,
            double pH,
            double ggamma,
            double Abar,
            double r,
            double tthetaLthres,
            double tthetaHthres,
            double tthetaRthreshold,
            double bRthreshold,
            double Rsubsidy,
            double tax,
            double ZHPAR,
            double ZLPAR, double T, double S){
    
    
    //0. Block of explanation of parameters-> All defined previously
    //   in function Zs
    
    
    //1. Computing intermediate functions.
    //1.1. Vector of zpitted elements
    vector<double> ZSpit=Zs(aalpha1, aalpha2, RMCH, RMCL, zHof, zLof, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    
    //1.2 Computing each error squared
    double zher1=pow((ZSpit[0]-zHof),2);
    double zher2=pow((ZSpit[1]-zLof),2);
    if(1==2){
        cout << " ----errz---- " << endl;
        cout << ZSpit[0] << " ZSpit[0] " << endl;
        cout << ZSpit[1] << " ZSpit[1] " << endl;
        cout << zHof << " zHof " << endl;
        cout << zLof << " zLof " << endl;
    }
    //1.3 Computing the sum of squared errors
    double ERR=pow(zher1,2)+pow(zher2,2);
    return(ERR);
}


//Defining structure of parameters to use
typedef struct FP_params{
    double PARaalpha1;
    double PARaalpha2;
    double PARRMCH;
    double PARRMCL;
    vector<double> PARbbgrid;
    vector<double> PARttgrid;
    double PARssigma;
    double PARbbeta;
    double PARw;
    double PARpL;
    double PARpH;
    double PARggamma;
    double PARAbar;
    double PARr;
    double PARtthetaLthres;
    double PARtthetaHthres;
    double PARtthetaRthreshold;
    double PARbRthreshold;
    double PARRsubsidy;
    double PARtax;
    double PARZHPAR;
    double PARZLPAR;
    double PART;
    double PARS;
    ;}pricesolver;

double FP_rootminsquarederror(unsigned n, const double *x,
                              double *grad,
                              void *params){
    //0. Loading parameters
    struct FP_params *p=(struct FP_params *)params;
    double aalpha1=p->PARaalpha1;
    double aalpha2=p->PARaalpha2;
    double RMCH=p->PARRMCH;
    double RMCL=p->PARRMCL;
    vector<double> bbgrid=p->PARbbgrid;
    vector<double> ttgrid=p->PARttgrid;
    double ssigma=p->PARssigma;
    double bbeta=p->PARbbeta;
    double w=p->PARw;
    double pL=p->PARpL;
    double pH=p->PARpH;
    double ggamma=p->PARggamma;
    double Abar=p->PARAbar;
    double r=p->PARr;
    double tthetaLthres=p->PARtthetaLthres;
    double tthetaHthres=p->PARtthetaHthres;
    double tthetaRthreshold=p->PARtthetaRthreshold;
    double bRthreshold=p->PARbRthreshold;
    double Rsubsidy=p->PARRsubsidy;
    double tax=p->PARtax;
    double ZHPAR=p->PARZHPAR;
    double ZLPAR=p->PARZLPAR;
    double T=p->PART;
    double S=p->PARS;
    if (1==2){
        cout << " -------------------"<<endl;
        cout << " loading FP_rootminsquarederror"<<endl;
        cout <<  aalpha1<< " aalpha1 " << endl;
        cout << aalpha2<< " aalpha2"<< endl;
        cout << RMCH << " RMCH"<< endl;
        cout << RMCL << " RMCL " << endl;
        cout << bbgrid[10] << " bbgrid[11]"<< endl;
        cout << ttgrid[10] << " ttgrid[11]"<< endl;
        cout << ssigma << " ssigma" << endl;
        cout << bbeta << " bbeta " << endl;
        cout << w << " w  " << endl;
        cout << pL << " pL" << endl;
        cout << pH << " pH " << endl;
        cout << ggamma << " ggamma " << endl;
        cout << Abar << " Abar " << endl;
        cout <<  r << " r " << endl;
        cout <<  tthetaLthres << " tthetaLthres  "  << endl;
        cout <<  tthetaHthres  << " tthetaHthres  "  << endl;
        cout <<  tthetaRthreshold  << "tthetaRthreshold   " << endl;
        cout <<  bRthreshold  << "  bRthreshold " << endl;
        cout <<  Rsubsidy   << " Rsubsidy  "<< endl;
        cout <<  tax << " tax  " << endl;
        cout << x[0] << " ZLOW" << endl;
        cout << x[1]+x[0] << " x[1] (ZHIGHPART)" << endl;
        cout << " -------------------"<<endl;
    }
    double residual=errZ( aalpha1,  aalpha2,
                         RMCH,  RMCL,x[1]+x[0],x[0], bbgrid, ttgrid,ssigma,bbeta,
                         w,pL,pH,ggamma,Abar,r,tthetaLthres,tthetaHthres,
                         tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S);
    double arrengment=residual;
    
    if (isnan(residual)==1){
        //I will try to generate normal deviates as min in case we have NAN so that
        //optimizer doesn't get stuck here.
        
        arrengment=1000;
        //arrengment=100;
    }
    
    return(arrengment);
}


//Function to compute fixed point in zh and zl
vector<double>ZFP(double aalpha1, double aalpha2,
                  double RMCH, double RMCL,
                  vector<double> bbgrid,
                  vector<double> ttgrid,
                  double ssigma,
                  double bbeta,
                  double w,
                  double pL,
                  double pH,
                  double ggamma,
                  double Abar,
                  double r,
                  double tthetaLthres,
                  double tthetaHthres,
                  double tthetaRthreshold,
                  double bRthreshold,
                  double Rsubsidy,
                  double tax,
                  double ZHPAR,
                  double ZLPAR, double T, double S){
    
    
    //Loading the structure
    pricesolver parstructest={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pL,pH,ggamma,Abar,r,tthetaLthres,tthetaHthres,
        tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    
    //Set up the optimization algorrithm
    nlopt_opt opt;
    opt=nlopt_create(NLOPT_LN_NELDERMEAD,2);// Dimension 2. Algoritthm cobyla
    nlopt_set_min_objective(opt,FP_rootminsquarederror,&parstructest);
    nlopt_set_xtol_rel(opt, 1.0e-3); //Tolerance
    //const double tolerance=1.0e-5;
    double LB[2]={0.05,0};
    double UB[2]={100,100};
    nlopt_set_lower_bounds(opt, LB);
    nlopt_set_upper_bounds(opt, UB);
    double xtest[2]={};
    xtest[0]=0.7;
    xtest[1]= 0.1;
    //nlopt_set_xtol_abs(opt, &tolerance);
    //Initial guess
    //double xtest[2]={};
    //xtest[0]=log(2);
    //xtest[1]=log(0.5);
    //Starting the optimization algorithm
    double minf;
    nlopt_optimize(opt, xtest, &minf);
    //if (nlopt_optimize(opt, xtest, &minf) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //    printf("found minimum at f(%g) = %0.10g\n", xtest[0], minf);
    //}
    //Will say that it converges only if the investment is actually positive
    
    
    
    //2.1 Matrix identifying the fixed points found
    vector<double> F;
    F.resize(3);
    F[0]=xtest[0];
    F[1]=xtest[1];
    F[2]=minf;
    if(1==2){
        cout << " in FP finding " << endl;
        cout << F[0] << " F[0]"<< endl;
        cout << F[1] << " F[1]"<< endl;
        cout << F[2] << " error " << endl;
    }
    return(F);
}

//Objective function of university
double zOBJ_HIGH(double aalpha1, double aalpha2,
                 double RMCH, double RMCL,
                 vector<double> bbgrid,
                 vector<double> ttgrid,
                 double ssigma,
                 double bbeta,
                 double w,
                 double pL,
                 double pH,
                 double ggamma,
                 double Abar,
                 double r,
                 double tthetaLthres,
                 double tthetaHthres,
                 double tthetaRthreshold,
                 double bRthreshold,
                 double Rsubsidy,
                 double tax,
                 double ZHPAR,
                 double ZLPAR, double T, double S){
    
    //First compute the fixed point and find ZH and ZL
    vector<double> F;
    F.resize(3);
    
    
    F=ZFP(aalpha1, aalpha2, RMCH, RMCL, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    double zLFound=F[0];
    double zHFound=F[0]+F[1];
    //Finding the elements from the consumers
    vector<double> Consumers;
    Consumers.resize(6);
    Consumers=aaverages(zHFound, zLFound, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
    double tthetaHbar=Consumers[0];
    double bbHaverage=Consumers[1];
    double invbbHaverage=1/bbHaverage;
    double Nh=Consumers[2];
    double VH=VariableCHigh(Nh);
    double IH = (pH*Nh+(RMCH)-VH)/Nh;
    //
    if(Nh<0){
        IH=0;
    }
    
    double Answer=ZHPAR*pow(tthetaHbar,aalpha1)*pow(IH,aalpha2)*pow(invbbHaverage,1-aalpha1-aalpha2);
    if (isnan(Answer)==1){
        
        //I will try to generate normal deviates as min in case we have NAN so that
        //optimizer doesn't get stuck here.
        
        Answer=-100;
    }
    
    if(F[2]>=0.1){
        Answer=-50;
    }
    
    
    if(1==2){
        cout << " ---ZobjHIgh----- " << endl;
        cout << pH << " pH " << endl;
        cout << RMCH << " RMCH " << endl;
        cout << Nh << " Nh " << endl;
        cout << VH<< "VH" << endl;
        cout << IH << " IH " << endl;
        cout <<   tthetaHbar << " tthetaHbar"<< endl;
        cout << bbHaverage << " bbHaverage " << endl;
        cout << invbbHaverage << " invbbHaverage " << endl;
        cout << F[0] << " F[0]  " << endl;
        cout << F[1] << " F[1] " << endl;
        cout << F[2] << " error in objhigh" << endl;
    }
    
    
    
   
    
    //Finally, if error is too large, we will not allow that to happen so set a very low objective function
    
    if (1==2){
        cout << " ---- in objHihghhh  ----- " << endl;
        cout << tthetaHthres << " tthetaHthres " << endl;
        cout << pH << " pH" << endl;
        cout << Answer  <<  " Answer MAXIMIZER  " << endl;
    }
    return(Answer);
    
}



//Objective function of low university
double zOBJ_LOW(double aalpha1, double aalpha2,
                double RMCH, double RMCL,
                vector<double> bbgrid,
                vector<double> ttgrid,
                double ssigma,
                double bbeta,
                double w,
                double pL,
                double pH,
                double ggamma,
                double Abar,
                double r,
                double tthetaLthres,
                double tthetaHthres,
                double tthetaRthreshold,
                double bRthreshold,
                double Rsubsidy,
                double tax,
                double ZHPAR,
                double ZLPAR, double T, double S){
    
    //First compute the fixed point and find ZH and ZL
    vector<double> F;
    F.resize(3);
    F=ZFP(aalpha1, aalpha2, RMCH, RMCL, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    double zLFound=F[0];
    double zHFound=F[0]+F[1];
    //Finding the elements from the consumers
    vector<double> Consumers;
    Consumers.resize(6);
    Consumers=aaverages(zHFound, zLFound, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
    double tthetaLbar=Consumers[3];
    double bbLaverage=Consumers[4];
    double invbbLaverage=1/bbLaverage;
    double Nl=Consumers[5];
    double VL=VariableCLow(Nl);
    double IL = (pL*Nl+(RMCL)-VL)/Nl;
    //
    if(Nl<0){
        IL=0;
    }
    
    double Answer=ZLPAR*pow(tthetaLbar,aalpha1)*pow(IL,aalpha2)*pow(invbbLaverage,1-aalpha1-aalpha2);
    if (isnan(Answer)==1){
        
        
        Answer=-100;
        //Answer=-10;
    }
    
    if(F[2]>=0.1){
        Answer=-50;
    }
    
    if(1==2){
        cout << " ------inOBJLOW------ " << endl;
        cout << zHFound << " zHfound " << endl;
        cout << zLFound << " zLFound " << endl;
        cout << bbgrid[0] << " bbgrid[0] " << endl;
        cout << ttgrid[0] << " ttgrid[0] " << endl;
        cout << ssigma << " ssigma " << endl;
        cout << bbeta << " bbeta " << endl;
        cout << w << " w " << endl;
        cout << pL << " pl " << endl;
        cout << pH << " ph " << endl;
        cout << ggamma << " ggamma " << endl;
        cout << Abar << " Abar " << endl;
        cout << r << " r " << endl;
        cout << tthetaLthres << " tthetaLthres " << endl;
        cout << tthetaHthres << " tthetaHthres " << endl;
        cout << tthetaRthreshold << " tthetaRthreshold " << endl;
        cout << bRthreshold << " bRthreshold " << endl;
        cout << Rsubsidy << " Rsubsidy " << endl;
        cout << tax << " tax " << endl;
        cout << IL << " I-L " << endl;
        cout << tthetaLbar << " tthetaLbar " << endl;
        cout << invbbLaverage << " invbbLaverage " << endl;
        cout << F[2] << " error in objlow" << endl;
        cout << Answer << " Answer " << endl;
    }
    
    
   
    return(Answer);
}

//Functions to optimize: Given a whole set of parameters, function of price and tthetabara.
//First I need to define a new structure
//Defining structure of parameters to use
typedef struct MAXZ_params{
    double PARaalpha1;
    double PARaalpha2;
    double PARRMCH;
    double PARRMCL;
    vector<double> PARbbgrid;
    vector<double> PARttgrid;
    double PARssigma;
    double PARbbeta;
    double PARw;
    double PARpOponnent;
    double PARggamma;
    double PARAbar;
    double PARr;
    double PARtthetathresOpponent;
    double PARtthetaRthreshold;
    double PARbRthreshold;
    double PARRsubsidy;
    double PARtax;
    double PARZHPAR;
    double PARZLPAR;
    double PART;
    double PARS;
    ;}optvalues;

//Objective of the high university
double zOBS_HIGH_maximizer(unsigned n, const double *x,
                           double *grad,
                           void *params){
    //0. Loading parameters
    struct MAXZ_params *p=(struct MAXZ_params *)params;
    double aalpha1=p->PARaalpha1;
    double aalpha2=p->PARaalpha2;
    double RMCH=p->PARRMCH;
    double RMCL=p->PARRMCL;
    vector<double> bbgrid=p->PARbbgrid;
    vector<double> ttgrid=p->PARttgrid;
    double ssigma=p->PARssigma;
    double bbeta=p->PARbbeta;
    double w=p->PARw;
    double pL=p->PARpOponnent;
    double ggamma=p->PARggamma;
    double Abar=p->PARAbar;
    double r=p->PARr;
    double tthetaLthres=p->PARtthetathresOpponent;
    double tthetaRthreshold=p->PARtthetaRthreshold;
    double bRthreshold=p->PARbRthreshold;
    double Rsubsidy=p->PARRsubsidy;
    double tax=p->PARtax;
    double ZHPAR=p->PARZHPAR;
    double ZLPAR=p->PARZLPAR;
    double T = p->PART;
    double S = p->PARS;
    
    //Loagind inputs
    //Loagind inputs
    double priceOWN=x[0];
    double tthetaOWNTHRESHOLD=x[1];
    
    //double priceOWN=exp(x[0]);
    //double tthetaOWNTHRESHOLD=exp(x[1])/(1+exp(x[1]));
    
    //Getting the functinon
    double OBJHIGH=zOBJ_HIGH(aalpha1, aalpha2, RMCH, RMCL, bbgrid,
                             ttgrid, ssigma, bbeta, w, pL, priceOWN,
                             ggamma, Abar, r, tthetaLthres,
                             tthetaOWNTHRESHOLD, tthetaRthreshold, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    
    //cout << "-----"<<endl;
    if(1==2){
        
        cout << " ------------ " << endl;
        cout << "------"<<endl;
        cout << " inobjhigh" << endl;
        cout << bbgrid[4] << " bbgrid[4] " << endl;
        cout << ttgrid[4] << " ttgrid[4] " << endl;
        cout << ssigma << " ssigma " << endl;
        cout << bbeta << " bbeta " << endl;
        cout << w << " w " << endl;
        cout << pL << " pl " << endl;
        cout << ggamma << " ggamma " << endl;
        cout << Abar << " Abar " << endl;
        cout << r << " r " << endl;
        cout << tthetaLthres << " tthetaLthres " << endl;
        cout << tthetaRthreshold << " tthetaRthreshold " << endl;
        cout << bRthreshold << " bRthreshold " << endl;
        cout << Rsubsidy << " Rsubsidy " << endl;
        cout << tax << " tax " << endl;
        cout << priceOWN << " priceOWN " << endl;
        cout << tthetaOWNTHRESHOLD << " tthetaOWNTHRESHOLD "<<endl;
        
        cout << OBJHIGH << " OBJHIGH "<< endl;
    }
    
    
    return(OBJHIGH);
}



//Objective of the low university
double zOBS_LOW_maximizer(unsigned n, const double *x,
                          double *grad,
                          void *params){
    //0. Loading parameters
    struct MAXZ_params *p=(struct MAXZ_params *)params;
    double aalpha1=p->PARaalpha1;
    double aalpha2=p->PARaalpha2;
    double RMCH=p->PARRMCH;
    double RMCL=p->PARRMCL;
    vector<double> bbgrid=p->PARbbgrid;
    vector<double> ttgrid=p->PARttgrid;
    double ssigma=p->PARssigma;
    double bbeta=p->PARbbeta;
    double w=p->PARw;
    double pH=p->PARpOponnent;
    double ggamma=p->PARggamma;
    double Abar=p->PARAbar;
    double r=p->PARr;
    double tthetaHthres=p->PARtthetathresOpponent;
    double tthetaRthreshold=p->PARtthetaRthreshold;
    double bRthreshold=p->PARbRthreshold;
    double Rsubsidy=p->PARRsubsidy;
    double tax=p->PARtax;
    double ZHPAR=p->PARZHPAR;
    double ZLPAR=p->PARZLPAR;
    double T = p->PART;
    double S = p->PARS;
    
    //Loagind inputs
    double priceOWN=x[0];
    double tthetaOWNTHRESHOLD=x[1];
    
    //double priceOWN=exp(x[0]);
    //double tthetaOWNTHRESHOLD=exp(x[1])/(1+exp(x[1]));
    
    //Getting the functinon
    double OBJLOW=zOBJ_LOW(aalpha1, aalpha2, RMCH, RMCL, bbgrid,
                           ttgrid, ssigma, bbeta, w, priceOWN, pH,
                           ggamma, Abar, r, tthetaOWNTHRESHOLD,
                           tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    if(1==2){
        
        cout << " ------------ " << endl;
        cout << "------"<<endl;
        cout << " inobjlow" << endl;
        cout << bbgrid[4] << " bbgrid[4] " << endl;
        cout << ttgrid[4] << " ttgrid[4] " << endl;
        cout << ssigma << " ssigma " << endl;
        cout << bbeta << " bbeta " << endl;
        cout << w << " w " << endl;
        cout << pH << " ph " << endl;
        cout << ggamma << " ggamma " << endl;
        cout << Abar << " Abar " << endl;
        cout << r << " r " << endl;
        cout << tthetaHthres << " tthetaHthres " << endl;
        cout << tthetaRthreshold << " tthetaRthreshold " << endl;
        cout << bRthreshold << " bRthreshold " << endl;
        cout << Rsubsidy << " Rsubsidy " << endl;
        cout << tax << " tax " << endl;
        cout << priceOWN << " priceOWN " << endl;
        cout << tthetaOWNTHRESHOLD << " tthetaOWNTHRESHOLD "<<endl;
        cout << OBJLOW << " OBJLOW "<< endl;
    }
    
    return(OBJLOW);
}

//Finding the optimal level of prices and ttheta for HIGH
vector<double> optHIGH(double aalpha1, double aalpha2,
                       double RMCH, double RMCL,
                       vector<double> bbgrid,
                       vector<double> ttgrid,
                       double ssigma,
                       double bbeta,
                       double w,
                       double pL,
                       double ggamma,
                       double Abar,
                       double r,
                       double tthetaLthres,
                       double tthetaRthreshold,
                       double bRthreshold,
                       double Rsubsidy,
                       double tax,
                       double ZHPAR,
                       double ZLPAR, double T, double S){
    
    //Loading the structure:
    optvalues paroptHIGH={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pL,ggamma,Abar,r,tthetaLthres,
        tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    
    //Set up the optimization algorrithm
    nlopt_opt opt;
    opt=nlopt_create(NLOPT_LN_NELDERMEAD,2);// Dimension 1. Algoritthm cobyla
    nlopt_set_max_objective(opt,zOBS_HIGH_maximizer,&paroptHIGH);
    nlopt_set_xtol_rel(opt, 1.0e-3); //Tolerance
    double LB[2]={0,0};
    double UB[2]={50,1};
    nlopt_set_lower_bounds(opt, LB);
    nlopt_set_upper_bounds(opt, UB);
    // const double TOL=1.0e-3;
    //nlopt_set_xtol_abs(opt,&TOL);
    //Initial guess

    double xtest[2]={};
    xtest[0]=0.3;
    xtest[1]=0.1;
    
    //xtest[0]=log(0.5);
    //xtest[1]=   .0;
    //Starting the optimization algorithm
    double minf;
    nlopt_optimize(opt, xtest, &minf);
    //if (nlopt_optimize(opt, xtest, &minf) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //    printf("found minimum at f(%g) = %0.10g\n", xtest[0], minf);
    //}
    vector<double> ANS;
    ANS.resize(3);
    ANS[0]=xtest[0];
    ANS[1]=xtest[1];
    ANS[2]=minf;
    return(ANS);
}



//Finding the optimal level of prices and ttheta for LOW
vector<double> optLOW(double aalpha1, double aalpha2,
                      double RMCH, double RMCL,
                      vector<double> bbgrid,
                      vector<double> ttgrid,
                      double ssigma,
                      double bbeta,
                      double w,
                      double pH,
                      double ggamma,
                      double Abar,
                      double r,
                      double tthetaHthres,
                      double tthetaRthreshold,
                      double bRthreshold,
                      double Rsubsidy,
                      double tax,
                      double ZHPAR,
                      double ZLPAR, double T, double S){
    
    //Loading the structure:
    optvalues paroptLOW={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pH,ggamma,Abar,r,tthetaHthres,
        tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    
    //Set up the optimization algorrithm
    nlopt_opt opt;
    opt=nlopt_create(NLOPT_LN_NELDERMEAD,2);// Dimension 1. Algoritthm cobyla
    nlopt_set_max_objective(opt,zOBS_LOW_maximizer,&paroptLOW);
    nlopt_set_xtol_rel(opt, 1.0e-3); //Tolerance
    double LB[2]={0,0};
    double UB[2]={pH,tthetaHthres};
    nlopt_set_lower_bounds(opt, LB);
    nlopt_set_upper_bounds(opt, UB);
    // const double TOL=1.0e-3;
    //nlopt_set_xtol_abs(opt,&TOL);
    //Initial guess


    double xtest[2]={};
    xtest[0]=0.3;
    xtest[1]= 0.1;
    //Starting the optimization algorithm
    double minf;
    nlopt_optimize(opt, xtest, &minf);
    //if (nlopt_optimize(opt, xtest, &minf) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //    printf("found minimum at f(%g) = %0.10g\n", xtest[0], minf);
    //}
    //cout << minf << endl;
    vector<double> ANS;
    ANS.resize(3);
    ANS[0]=xtest[0];
    ANS[1]=xtest[1];
    ANS[2]=minf;
    return(ANS);
}

//Finding the nash equilibrium
vector<double> NASHEQUILIBRIUM(double aalpha1, double aalpha2,
                               double RMCH, double RMCL,
                               vector<double> bbgrid,
                               vector<double> ttgrid,
                               double ssigma,
                               double bbeta,
                               double w,
                               double pL,
                               double pH,
                               double ggamma,
                               double Abar,
                               double r,
                               double tthetaLthres,
                               double tthetaHthres,
                               double tthetaRthreshold,
                               double bRthreshold,
                               double Rsubsidy,
                               double tax,
                               double ZHPAR,
                               double ZLPAR,
                               int WRITE, double T, double S){
    
    //Write is an integer. If set equal to one, it will write a csv file with the output of the NEQ.
    
    
    //Setting up the tolerance level
    double TOL=1.0e-5;
    double error=10;
    
    
    //Defining the structures
    //Defining the parameters
    optvalues paroptHIGH={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pL,ggamma,Abar,r,tthetaLthres,
        tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    optvalues paroptLOW={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pH,ggamma,Abar,r,tthetaHthres,
        tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    //Defining old and new guess vector
    double NewPH=0;
    double NewPL=0;
    double NewTH=0;
    double NewTL=0;
    
    double VH=0;
    double VL=0;
    double ZLPARUSED=ZLPAR;
    double ZHPARUSED=ZHPAR;
    
    //Loading up the parameters
    vector<double>NEQ;
    NEQ.resize(5);
    
    int iterator=1;
    
    while(error>TOL){
        cout << iterator << " -----iterator----- " << endl;
        //Finding optimal parameters:
        //if(iterator % 2 == 0){
        
        
        
        //Update zhparused
        //ZHPARUSED=ZHPAR;
        vector<double> optH=optHIGH(aalpha1, aalpha2, RMCH, RMCL, bbgrid,ttgrid, ssigma, bbeta, w, pL,ggamma, Abar, r, tthetaLthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax,ZHPARUSED,ZLPARUSED, T, S);
        
        NewPH= optH[0];
        NewTH= optH[1];
        VH=optH[2];
        //if(VH<0){
        //    ZHPARUSED=0;
        //}
        
        
        
        
        
        //}
        
        //if(iterator % 2 == 1){
        //ZLPARUSED=ZLPAR;
        vector<double> optL=optLOW(aalpha1, aalpha2, RMCH, RMCL, bbgrid,ttgrid, ssigma, bbeta, w, pH,ggamma, Abar, r, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax,ZHPARUSED,ZLPARUSED, T, S);
        
        NewPL=optL[0] ;
        NewTL=optL[1];
        VL=optL[2];
        
        
        
        //if(VL<0){
        //    ZLPARUSED=0;
        //}
        error=pow(NewPL-pL,2)+pow(NewPH-pH,2)+pow(NewTH-tthetaHthres,2)+pow(NewTL-tthetaLthres,2);
        
        cout << "--------------------------"<<endl;
        cout << error << " error " << endl;
        cout << pL << " pL old " << endl;
        cout << pH << " pH old " << endl;
        cout << NewTH <<  " New Th " << endl;
        cout << NewTL << " Newtl" << endl;
        cout << NewPL << " NEWPL"<<endl;
        cout << NewPH << " NEWPH"<<endl;
        cout << tthetaLthres << " tthetaLthres"<<endl;
        cout << tthetaHthres << " tthetaHthres"<<endl;
        cout << VL << " VL " << endl;
        cout << VH << " VH " << endl;
        cout << iterator << " iterator " << endl;
        cout << "--------------------------"<<endl;
        
        
        pL=NewPL;
        pH=NewPH;
        tthetaLthres=NewTL;
        tthetaHthres=NewTH;
        iterator++;
    }
    
    NEQ[0]=NewPL;
    NEQ[1]=NewTL;
    NEQ[2]=NewPH;
    NEQ[3]=NewTH;
    NEQ[4]=error;
    
    //I will write the NEQ in a file.
    
    if (WRITE==1){
        ofstream optparam("NashEquilibrium.csv");
        
        optparam<<NEQ[0] << endl;
        optparam<<NEQ[1] << endl;
        optparam<<NEQ[2] << endl;
        optparam<<NEQ[3] << endl;
        
        
        optparam.close();
    }
    
    
    
    
    return(NEQ);
}



//Finding the nash equilibrium->Using the guess from computer rather than contraction mapping style. 
double NASHERROR(double aalpha1, double aalpha2,
                               double RMCH, double RMCL,
                               vector<double> bbgrid,
                               vector<double> ttgrid,
                               double ssigma,
                               double bbeta,
                               double w,
                               double pL,
                               double pH,
                               double ggamma,
                               double Abar,
                               double r,
                               double tthetaLthres,
                               double tthetaHthres,
                               double tthetaRthreshold,
                               double bRthreshold,
                               double Rsubsidy,
                               double tax,
                               double ZHPAR,
                               double ZLPAR, double T, double S){
    


    //Defining the structures
    optvalues paroptHIGH={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pL,ggamma,Abar,r,tthetaLthres,
        tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    optvalues paroptLOW={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pH,ggamma,Abar,r,tthetaHthres,
        tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    //Defining old and new guess vector
    double NewPH=0;
    double NewPL=0;
    double NewTH=0;
    double NewTL=0;
    
    double VH=0;
    double VL=0;
    double ZLPARUSED=ZLPAR;
    double ZHPARUSED=ZHPAR;
    
    //Loading up the parameters
    //Finding optimized parameters for high university
    vector<double> optH=optHIGH(aalpha1, aalpha2, RMCH, RMCL, bbgrid,ttgrid, ssigma, bbeta, w, pL,ggamma, Abar, r, tthetaLthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax,ZHPARUSED,ZLPARUSED, T, S);
        
    NewPH= optH[0];
    NewTH= optH[1];
    VH=optH[2];

        
        
    //Finding optimized parameters for low university
    vector<double> optL=optLOW(aalpha1, aalpha2, RMCH, RMCL, bbgrid,ttgrid, ssigma, bbeta, w, pH,ggamma, Abar, r, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax,ZHPARUSED,ZLPARUSED, T, S);
        
    NewPL=optL[0] ;
    NewTL=optL[1];
    VL=optL[2];
    
    //And finally finding the error term between the guessed policy and the found one
    double error=0;
    error=pow(NewPL-pL,2)+pow(NewPH-pH,2)+pow(NewTH-tthetaHthres,2)+pow(NewTL-tthetaLthres,2);
    if(VL<0 || VH<0 ){
        error=10;
    }
    cout << "--------------------------"<<endl;
    cout << error << " error " << endl;
    cout << pL << " pL old " << endl;
    cout << pH << " pH old " << endl;
    cout << NewTH <<  " New Th " << endl;
    cout << NewTL << " Newtl" << endl;
    cout << NewPL << " NEWPL"<<endl;
    cout << NewPH << " NEWPH"<<endl;
    cout << tthetaLthres << " tthetaLthres"<<endl;
    cout << tthetaHthres << " tthetaHthres"<<endl;
    cout << VL << " VL " << endl;
    cout << VH << " VH " << endl;
    cout << "--------------------------"<<endl;
    
    return(error);
}


//Define the structure used to find the nash in the end

typedef struct NASHERROR_params{
    double PARaalpha1;
    double PARaalpha2;
    double PARRMCH;
    double PARRMCL;
    vector<double> PARbbgrid;
    vector<double> PARttgrid;
    double PARssigma;
    double PARbbeta;
    double PARw;
    double PARggamma;
    double PARAbar;
    double PARr;
    double PARtthetaRthreshold;
    double PARbRthreshold;
    double PARRsubsidy;
    double PARtax;
    double PARZHPAR;
    double PARZLPAR;
    double PART;
    double PARS;
    ;}NASHERRORSTRUCTURE;







//Objective of the low university

double NASHERROROPTIMIZERINTERM(unsigned n, const double *x,
                          double *grad,
                          void *params){
    //0. Loading parameters
    struct NASHERROR_params *p=(struct NASHERROR_params *)params;
    double aalpha1=p->PARaalpha1;
    double aalpha2=p->PARaalpha2;
    double RMCH=p->PARRMCH;
    double RMCL=p->PARRMCL;
    vector<double> bbgrid=p->PARbbgrid;
    vector<double> ttgrid=p->PARttgrid;
    double ssigma=p->PARssigma;
    double bbeta=p->PARbbeta;
    double w=p->PARw;
    double ggamma=p->PARggamma;
    double Abar=p->PARAbar;
    double r=p->PARr;
    double tthetaRthreshold=p->PARtthetaRthreshold;
    double bRthreshold=p->PARbRthreshold;
    double Rsubsidy=p->PARRsubsidy;
    double tax=p->PARtax;
    double ZHPAR=p->PARZHPAR;
    double ZLPAR=p->PARZLPAR;
    double T = p->PART;
    double S = p->PARS;
    
    //Loagind inputs
    double pL=x[0];
    double pH=x[1];
    double tthetaLthres=x[2];
    double tthetaHthres=x[3];
    
    //Finding the nash error
    double ans=0;
    cout << " here  " << endl;
    ans=NASHERROR(aalpha1, aalpha2, RMCH, RMCL,bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthreshold, bRthreshold,Rsubsidy, tax, ZHPAR, ZLPAR, T, S);
    return(ans);
}



//Finally, I will use the optimizer to find the NASH by minimizing the squared error
vector<double> NASHERRORFINDER(double aalpha1, double aalpha2,
                 double RMCH, double RMCL,
                 vector<double> bbgrid,
                 vector<double> ttgrid,
                 double ssigma,
                 double bbeta,
                 double w,
                 double ggamma,
                 double Abar,
                 double r,
                 double tthetaRthreshold,
                 double bRthreshold,
                 double Rsubsidy,
                 double tax,
                 double ZHPAR,
                 double ZLPAR, double T, double S){
    
    //Load the structure used for Nash
    NASHERRORSTRUCTURE NASHERROR_params={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,ggamma,Abar,r,tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    //Setting the optimizer
    nlopt_opt opt;
    opt=nlopt_create(NLOPT_LN_NELDERMEAD,4); //Algorithm and dimensions
    nlopt_set_min_objective(opt,NASHERROROPTIMIZERINTERM,&NASHERROR_params);
    nlopt_set_xtol_rel(opt, 1.0e-3); //Tolerance
    
    
    //Setting the initial levels
    double PLINIT=0.3;
    double PHINIT=0.3;
    double TTHETALINIT=0.1;
    double TTHETAHINIT=0.1;
    
    //Loading initial values
    double xtest[4]={};
    xtest[0]=PLINIT;
    xtest[1]=PHINIT;
    xtest[2]=TTHETALINIT;
    xtest[3]=TTHETAHINIT;
    
    //Setting the desired bounds
    double LB[4]={0,0,0,0};
    double UB[4]={PHINIT,50,1,1};
    nlopt_set_lower_bounds(opt, LB);
    nlopt_set_upper_bounds(opt, UB);
    
    //Start the optimization routine
    double minf;
    nlopt_optimize(opt,xtest,&minf);
    vector<double>NASHEQUILIBRIUM;
    NASHEQUILIBRIUM.resize(4);
    NASHEQUILIBRIUM[0]=xtest[0];
    NASHEQUILIBRIUM[1]=xtest[1];
    NASHEQUILIBRIUM[2]=xtest[2];
    NASHEQUILIBRIUM[3]=xtest[3];
    
    //Writing the neq if we want

    ofstream optparam("NashEquilibrium.csv");
        
    optparam<<NASHEQUILIBRIUM[0] << endl;
    optparam<<NASHEQUILIBRIUM[1] << endl;
    optparam<<NASHEQUILIBRIUM[2] << endl;
    optparam<<NASHEQUILIBRIUM[3] << endl;
    optparam.close();

    return(NASHEQUILIBRIUM);
}


//-----------------------------------------//
//FUNCIONES PARA VERIFICAR QUE MANDA ZARRUK//
//-----------------------------------------//


//Objective function of university
vector<double> QUALITIES(double aalpha1, double aalpha2,
                         double RMCH, double RMCL,
                         vector<double> bbgrid,
                         vector<double> ttgrid,
                         double ssigma,
                         double bbeta,
                         double w,
                         double pL,
                         double pH,
                         double zL,
                         double zH,
                         double ggamma,
                         double Abar,
                         double r,
                         double tthetaLthres,
                         double tthetaHthres,
                         double tthetaRthreshold,
                         double bRthreshold,
                         double Rsubsidy,
                         double tax,
                         double ZHPAR,
                         double ZLPAR, double T, double S){
    
    //Finding the elements from the consumers
    vector<double> Consumers;
    Consumers.resize(6);
    Consumers=aaverages(zH, zL, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
    double tthetaHbar=Consumers[0];
    double bbHaverage=Consumers[1];
    double Nh=Consumers[2];
    double VH=VariableCHigh(Nh);
    double IH = (pH*Nh+(RMCH)-VH)/Nh;
    //
    if(Nh<0){
        IH=0;
    }
    
    double tthetaLbar=Consumers[3];
    double bbLaverage=Consumers[4];
    double Nl = Consumers[5];
    double VL = VariableCLow(Nl);
    double IL = (pL*Nl+(RMCL)-VL)/Nl;
    //
    if(Nl < 0){
        IL = 0;
    }
    
    double Answer1=ZHPAR*pow(tthetaHbar,aalpha1)*pow(IH,aalpha2);
    double Answer2=ZLPAR*pow(tthetaLbar,aalpha1)*pow(IL,aalpha2);
    vector<double> Answer;
    Answer.resize(2);
    Answer[0] = Answer1;
    Answer[1] = Answer2;
    
    return(Answer);
    
}


//Finding the optimal level of prices and ttheta for HIGH
vector<double> optHIGHverify(double aalpha1, double aalpha2,
                             double RMCH, double RMCL,
                             vector<double> bbgrid,
                             vector<double> ttgrid,
                             double ssigma,
                             double bbeta,
                             double w,
                             double pL,
                             double ggamma,
                             double Abar,
                             double r,
                             double tthetaLthres,
                             double tthetaRthreshold,
                             double bRthreshold,
                             double Rsubsidy,
                             double tax,
                             double ZHPAR,
                             double ZLPAR, double PHinitial, double tHinitial, double T, double S){
    
    //Loading the structure:
    optvalues paroptHIGH={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pL,ggamma,Abar,r,tthetaLthres,
        tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    
    //Set up the optimization algorrithm
    nlopt_opt opt;
    opt=nlopt_create(NLOPT_LN_NELDERMEAD,2);// Dimension 1. Algoritthm cobyla
    nlopt_set_max_objective(opt,zOBS_HIGH_maximizer,&paroptHIGH);
    nlopt_set_xtol_rel(opt, 1.0e-3); //Tolerance
    double LB[2]={0,0};
    double UB[2]={50,1};
    nlopt_set_lower_bounds(opt, LB);
    nlopt_set_upper_bounds(opt, UB);
    // const double TOL=1.0e-3;
    //nlopt_set_xtol_abs(opt,&TOL);
    //Initial guess
    
    double xtest[2]={};
    xtest[0]=PHinitial;
    xtest[1]=tHinitial;
    
    //xtest[0]=log(0.5);
    //xtest[1]=   .0;
    //Starting the optimization algorithm
    double minf;
    nlopt_optimize(opt, xtest, &minf);
    //if (nlopt_optimize(opt, xtest, &minf) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //    printf("found minimum at f(%g) = %0.10g\n", xtest[0], minf);
    //}
    vector<double> ANS;
    ANS.resize(3);
    ANS[0]=xtest[0];
    ANS[1]=xtest[1];
    ANS[2]=minf;
    return(ANS);
}

//Finding the optimal level of prices and ttheta for LOW
vector<double> optLOWverify(double aalpha1, double aalpha2,
                            double RMCH, double RMCL,
                            vector<double> bbgrid,
                            vector<double> ttgrid,
                            double ssigma,
                            double bbeta,
                            double w,
                            double pH,
                            double ggamma,
                            double Abar,
                            double r,
                            double tthetaHthres,
                            double tthetaRthreshold,
                            double bRthreshold,
                            double Rsubsidy,
                            double tax,
                            double ZHPAR,
                            double ZLPAR, double PLinitial, double tLinitial, double T, double S){
    
    //Loading the structure:
    optvalues paroptLOW={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pH,ggamma,Abar,r,tthetaHthres,
        tthetaRthreshold,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    
    //Set up the optimization algorrithm
    nlopt_opt opt;
    opt=nlopt_create(NLOPT_LN_NELDERMEAD,2);// Dimension 1. Algoritthm cobyla
    nlopt_set_max_objective(opt,zOBS_LOW_maximizer,&paroptLOW);
    nlopt_set_xtol_rel(opt, 1.0e-3); //Tolerance
    double LB[2]={0,0};
    double UB[2]={pH,tthetaHthres};
    nlopt_set_lower_bounds(opt, LB);
    nlopt_set_upper_bounds(opt, UB);
    // const double TOL=1.0e-3;
    //nlopt_set_xtol_abs(opt,&TOL);
    //Initial guess
    
    
    double xtest[2]={};
    xtest[0]=PLinitial;
    xtest[1]= tLinitial;
    //Starting the optimization algorithm
    double minf;
    nlopt_optimize(opt, xtest, &minf);
    //if (nlopt_optimize(opt, xtest, &minf) < 0) {
    //    printf("nlopt failed!\n");
    //}
    //else {
    //    printf("found minimum at f(%g) = %0.10g\n", xtest[0], minf);
    //}
    //cout << minf << endl;
    vector<double> ANS;
    ANS.resize(3);
    ANS[0]=xtest[0];
    ANS[1]=xtest[1];
    ANS[2]=minf;
    return(ANS);
}


//Function to compute the moments of the model
vector<double>MOMENTS(double aalpha1, double aalpha2,
                      double RMCH, double RMCL,
                      vector<double> bbgrid,
                      vector<double> ttgrid,
                      double ssigma,
                      double bbeta,
                      double w,
                      double pL,
                      double pH,
                      double ggamma,
                      double Abar,
                      double r,
                      double tthetaLthres,
                      double tthetaHthres,
                      double tthetaRthreshold,
                      double bRthreshold,
                      double Rsubsidy,
                      double tax,
                      double ZHPAR,
                      double ZLPAR, double T, double S){
    
    
    int iterator=0; ///Variable used multiple times in this function
    
    //It should input values in equilibrium. The output will be a vector of theoretical moments derived from the model
    vector<double>MOMENTS;
    MOMENTS.resize(6);
    
    //First compute the fixed point and find ZH and ZL
    vector<double> F;
    F.resize(3);
    
    
    F=ZFP(aalpha1, aalpha2, RMCH, RMCL, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    double zLFound=F[0];
    double zHFound=F[0]+F[1];
    //Finding the elements from the consumers
    vector<double> Consumers;
    Consumers.resize(6);
    Consumers=aaverages(zHFound, zLFound, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
    double tthetaHbar=Consumers[0];
    double bbHaverage=Consumers[1];
    double invbbHaverage=1/bbHaverage;
    double Nh=Consumers[2];
    double VH=VariableCHigh(Nh);
    double IH = (pH*Nh+(RMCH)-VH)/Nh;
    double tthetaLbar=Consumers[3];
    double bbLaverage=Consumers[4];
    double invbbLaverage=1/bbLaverage;
    double Nl=Consumers[5];
    double VL=VariableCHigh(Nl);
    double IL = (pL*Nl+(RMCL)-VL)/Nl;
    
    
    
    //Conditional moments for decisions of the household
    int bbnumber=bbgrid.size(); //Size along b dim
    int ttnumber=ttgrid.size(); //Size along ttheta dimension.
    
    
    //Matrices MATGRIDh and MATGRIDl will be used to define conditional moments on a particular attendance decisions of the households
    
    vector<vector<double> > MATGRIDh;
    vector<vector<double> > MATGRIDl;
    MATGRIDh.resize(bbnumber);
    MATGRIDl.resize(bbnumber);
    vector<double> res;
    res.resize(2);
    for (int bb=0; bb<bbnumber; bb=bb+1){
        for (int tt=0; tt<ttnumber; tt=tt+1){
            res = decision(bbgrid[bb], ttgrid[tt], ssigma, bbeta, w, w*(1+zLFound), w*(1+zHFound), pL, pH, tthetaLthres, tthetaHthres, ggamma, r, Abar, tthetaRthreshold, bRthreshold, Rsubsidy, tax, T, S);
            MATGRIDh[bb][tt] = res[0];
            MATGRIDl[bb][tt] = res[1];
        }
    }
    
    //In order to compute the proportino of people that attend university h for a given subset of ttheta bbgrid:
    double CONDM1H=0;
    double CONDM1L=0;
    iterator=0;
    for (int bb=0; bb<bbnumber; bb++){
        for(int tt=0; tt<ttnumber;tt++){
            CONDM1H+=MATGRIDh[bb][tt];
            CONDM1L+=MATGRIDl[bb][tt];
            iterator++;
        }
    }
    
    CONDM1H=CONDM1H/iterator;
    CONDM1L=CONDM1L/iterator;
    
    //Loading moments
    MOMENTS[0]=IH;
    MOMENTS[1]=IL;
    MOMENTS[2]=Nh/(bbnumber*ttnumber); //PROPORTION OF STUDENTS IN HIGH UNIVERSITY
    MOMENTS[3]=Nl/(bbnumber*ttnumber); //Proportion of students in low university
    MOMENTS[4]=CONDM1H; //CONDITIONAL PROPORTION 1 OF STUDENTS IN HIGH UNIVERSITY
    MOMENTS[5]=CONDM1L; //CONDITIONAL PROPORTION 1 OF STUDENTS IN HIGH UNIVERSITY
    
    
    
    return(MOMENTS);
}

int main(int argc, const char * argv[])
{
    
    
    //0. Testing the functions
    
    //0.0 Setting the parameters to test the functions
    double aalpha1=0.211; //Initial 0.25 ; 0.25
    double aalpha2=0.358;
    double RMCH=-12;  //For high I did RMCH=-1500,RMCL=-2000. Initial 3 and 2.
    double RMCL=-7;
    double zHof=1.4;
    double zLof=1.2;
    
    //Generating bbgrid:
    vector<double> bbgrid;
    int bgridsize=100;
    bbgrid.resize(bgridsize);
    double bbgridmax=10;
    double bbgridmin=0;
    double bbstep=(bbgridmax-bbgridmin)/(bgridsize-1);
    double it=0;
    bbgrid.resize(bgridsize);
    it=0;
    
    for(int b=0; b<bgridsize; b++){
        
        bbgrid[b]=bbgridmin+it*bbstep;
        //cout << bbgrid[b] << " bbgrid[b]"<<endl;
        it++;
        
    }
    
    //Ttheta grid
    vector<double> ttgrid;
    int ttgridsize=200;
    ttgrid.resize(ttgridsize);
    double ttgridmax=1;
    double ttgridmin=0;
    double ttstep=(ttgridmax-ttgridmin)/(ttgridsize-1);
    it=0;
    for(int b=0; b<ttgridsize; b++){
        ttgrid[b]=ttgridmin+it*ttstep;
        //cout << ttgrid[b] << " ttgrid[b]"<<endl;
        it++;
    }
    
    //Verify parallelization
    
    
    int xtestparallel=50;
#pragma omp parallel for
    for (int x=0; x<xtestparallel; x=x+1){
       // cout << x << " x test parallel"<< endl;
    }
    
    double T = 80.0;
    double S = 5.0;
    double ssigma=2;
    double bbeta=0.971;
    double w=2;
    double pL=0.3;
    double pH=0.55;
    double ggamma=0.00;
    double Abar=2;
    double r=0.02;
    double tthetaLthres=0.1;
    double tthetaHthres=0.6;
    double tthetaRthres=0.7;
    double bRthreshold=5;
    double Rsubsidy=0.5; //Rsubsidy=0.5 and tax of 0.1 works fine
    double tax=0.1;
    double ZHPAR=1.4;  //90 AND 80 WORKED FINE
    double ZLPAR=1.2;
    //Defining the structure
    pricesolver parstructest={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pL,pH,ggamma,Abar,r,tthetaLthres,tthetaHthres,
        tthetaRthres,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    //0.1 Testing the functions
    
    
    //1 Aaverages
    vector<double> x;
    x.resize(6);
    //pL=1.43;
    //tthetaLthres=0.8;
    //pH=1.41;
    //tthetaHthres=0.8;
    
    x = aaverages(zHof, zLof, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthres, bRthreshold, Rsubsidy, tax, T, S);
    cout << " ======Averages======= " << endl;
    cout << zHof << " zHof " << endl;
    cout << zLof << " zLof " << endl;
    cout << " Test aaverages " << endl;
    cout << x[0]<< " tthetaHaverage"<< endl;
    cout << x[1]<< " bbHaverage"<< endl;
    cout << x[2]<< " Htotal"<< endl;
    cout << x[3]<< " tthetaLaverage"<< endl;
    cout << x[4]<< " BBLAVERAGE"<< endl;
    cout << x[5]<< " Ltotal"<< endl;
    cout << " ============= " << endl;
    

    //2 Zs
    vector<double> testZs=Zs(aalpha1, aalpha2, RMCH, RMCL, zHof, zLof, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    cout << " ============= " << endl;
    cout << " Test Zs " << endl;
    cout << testZs[0] << " testZs[0] " << endl;
    cout << testZs[1] << " testZs[1] " << endl;
    
    
    //3 test
    
    
    double testERZ=errZ(aalpha1, aalpha2, RMCH, RMCL, zHof, zLof, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    cout << " ============= " << endl;
    cout << " Test Zs " << endl;
    cout << testERZ<< " errZ " << endl;
    
    //4 Squared error for minimizer
    cout << " ============= " << endl;
    cout << " FPTEST" << endl;
    double FPTEST=FP_rootminsquarederror(1, &zLof, &zHof, &parstructest);
    
    cout << FPTEST<< " FP_rootminsquarederror " << endl;
    //5 Fixed point
    
    
    
    vector<double> F;
    F.resize(3);
    cout << " ============= " << endl;
    cout << " Test ZFP" << endl;
    cout << " ===INSIDE OF ZFP========== " << endl;

    
    F=ZFP(aalpha1, aalpha2, RMCH, RMCL, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    
    
    
    
    //Re normalizing
    double zLFound=F[0];
    double zHFound=F[0]+F[1];
    cout << zLFound <<"zLFound" << endl;
    cout << zHFound << "zHFound"<< endl;
    
    
    //Verifying the solution
    cout << " =========== " << endl;
    cout << " Verifying solution " << endl;
    double errV=errZ(aalpha1, aalpha2, RMCH, RMCL, zHFound, zLFound, bbgrid, ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    
    cout << errV << " Error verified " << endl;
    
    //Verifying objective function of high firm
    cout << "==========="<< endl;
    cout << " Obj High " << endl;
    double GainHigh=zOBJ_HIGH(aalpha1, aalpha2, RMCH, RMCL, bbgrid,ttgrid, ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    
    
    cout << GainHigh << " Objective of high university " << endl;
    
    cout << " ------ " << endl;
    double GainLow=zOBJ_LOW(aalpha1, aalpha2, RMCH, RMCL, bbgrid,ttgrid, ssigma, bbeta, w, pL, pH,ggamma, Abar, r, tthetaLthres,tthetaHthres, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    cout << GainLow << " Objective of low university "<< endl;
    
    
    //cout << "======= " << endl;
    //cout << " zOBS_high_maximizer"<< endl;
    vector<double> xZobsHigh;
    //Defining the parameters
    optvalues paroptHIGH={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pL,ggamma,Abar,r,tthetaLthres,
        tthetaRthres,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    optvalues paroptLOW={aalpha1,aalpha2,RMCH,RMCL,bbgrid,ttgrid,
        ssigma,bbeta,w,pH,ggamma,Abar,r,tthetaHthres,
        tthetaRthres,bRthreshold,Rsubsidy,tax,ZHPAR,ZLPAR, T, S};
    
    //Defining the initial inputs
    xZobsHigh.resize(2);
    xZobsHigh[0]=pH;
    xZobsHigh[1]=tthetaHthres;
    
    
    //Defining the output
    //cout << tthetaHthres << " tthetaHthres" << endl;
    //cout << &tthetaHthres << " &tthetaHthres" << endl;
    
    //This output not need to coincide with GainHigh:  &pH, &tthetaHthres
    //cout << " zobjoutputhigh  " << endl;
    //double PHTEST=log(0.5);
    //double THTHRES=50;
    //const double A[2]={log(0.5),-90000000};
    //double zobjOutputHIGH=zOBS_HIGH_maximizer(1, A, &PHTEST , &paroptHIGH);
    //cout << zobjOutputHIGH << " zobjOutputHIGH" << endl;
    
    //double zobjOutputLOW=zOBS_LOW_maximizer(2, &pL, &tthetaLthres, &paroptLOW);
    //cout << zobjOutputLOW << " zobjOutputLOW" << endl;
    
    //Optimal price and ttheta for low
    if(1==2){
    cout << " ------------------------" << endl;
    cout << " Start verifying optLow  " << endl;
    cout << " ------------------------" << endl;

    
    vector<double> optL=optLOW(aalpha1, aalpha2, RMCH, RMCL, bbgrid,ttgrid, ssigma, bbeta, w, pH,ggamma, Abar, r, tthetaHthres, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    

    
    cout << " optimization found " << endl;
    cout << optL[0]<< " price low opt " << endl;
    cout << optL[1] << " ttheta low opt " << endl;
    cout << optL[2]<< " value low " << endl;
    
    
    //Optimal price and ttheta for high
    
    cout << " ------------------------" << endl;
    cout << " Start verifying optHigh  " << endl;
    cout << " ------------------------" << endl;
    //double plTEST=5.58;
    //double tlTltest=0.985587;
    
    //double plTEST=7.38;
    //double tlTltest=0.985587;
    //vector<double> optH=optHIGH(aalpha1, aalpha2, RMCH, RMCL, bbgrid,ttgrid, ssigma, bbeta, w, pL,ggamma, Abar, r, tthetaLthres, tthetaRthres, bRthreshold, Rsubsidy, tax,30,10);
    //vector<double> optH=optHIGH(aalpha1, aalpha2, RMCH, RMCL, bbgrid,ttgrid, ssigma, bbeta, w, pL,ggamma, Abar, r, tthetaLthres, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR);
    
    vector<double>optH;
    optH.resize(3);
    optH=optHIGH(aalpha1, aalpha2, RMCH, RMCL, bbgrid,ttgrid, ssigma, bbeta, w, pL,ggamma, Abar, r, tthetaLthres, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    cout << " optimization found " << endl;
    cout << optH[0] << " price high opt " << endl;
    cout << optH[1] << " ttheta high opt " << endl;
    cout << optH[2] << " value high " << endl;
    
    }//Finish 1==2
    
    //Checking error nash
    cout << "--------------------- " << endl;
    cout << " Checking error Nash " << endl;
    cout <<  " -------------------- " << endl;
    double PLNASH=0.3;
    double PHNASH=0.3;
    double tthetaHnash=0.1;
    double tthetaLnash=0.1;

    cout <<  NASHERROR(aalpha1, aalpha2, RMCH, RMCL, bbgrid, ttgrid, ssigma, bbeta, w, PLNASH,PHNASH,ggamma, Abar, r, tthetaLnash,tthetaHnash,tthetaRthres, bRthreshold, Rsubsidy, tax, ZHPAR, ZLPAR, T, S) << " error " << endl;
    
    //Nash Equilibrium
    
    cout << "------------------------------------------ " << endl;
    cout << " Start solving NE via optimization routine " << endl;
    cout << " ----------------------------------------- " << endl;
    vector<double> NASHMACHINE=NASHERRORFINDER(aalpha1, aalpha2, RMCH, RMCL, bbgrid, ttgrid, ssigma, bbeta, w, ggamma, Abar, r, tthetaRthres, bRthreshold, Rsubsidy, tax, ZHPAR, ZLPAR, T, S);
    
    
    // -----------------------------------------------------------------------------------------//
    //          ESTA PARTE ES PARA CHEQUEAR QUE NO HAYA ERRORES EN LAS COMPUTACIONES            //
    // -----------------------------------------------------------------------------------------//
    
    double PL;
    double PH;
    double TTHETAL;
    double TTHETAH;
    
    PL = NASHMACHINE[0];
    PH = NASHMACHINE[1];
    TTHETAL = NASHMACHINE[2];
    TTHETAH = NASHMACHINE[3];
    
    cout << "  " << endl;
    cout << " ----------------------------------------------- " << endl;
    cout << " ----------    ERROR CHECKING    --------------- " << endl;
    cout << " ----------------------------------------------- " << endl;
    
    cout << "  " << endl;
    cout << " EQUILIBRIUM POLICIES ARE: " << endl;
    cout << "  " << endl;
    cout << "Precio de baja Pl: " << PL << endl;
    cout << "Precio de alta Ph: " << PH << endl;
    cout << "  " << endl;
    cout << "Theta de baja Thetal: " << TTHETAL << endl;
    cout << "Theta de alta Thetah: " << TTHETAH << endl;
    cout << "  " << endl;
    cout << "------------------------------------------------- " << endl;
    cout << "  " << endl;
    
    // Verifico manualmente que de lo mismo
    vector<double> optLow = optLOWverify(aalpha1, aalpha2, RMCH, RMCL, bbgrid, ttgrid, ssigma, bbeta, w, PH, ggamma, Abar, r, TTHETAH, tthetaRthres, bRthreshold, Rsubsidy, tax, ZHPAR, ZLPAR, PL, TTHETAL, T, S);
    vector<double> optHigh = optHIGHverify(aalpha1, aalpha2, RMCH, RMCL, bbgrid, ttgrid, ssigma, bbeta, w, PL, ggamma, Abar, r, TTHETAL, tthetaRthres, bRthreshold, Rsubsidy, tax, ZHPAR, ZLPAR, PH, TTHETAH, T, S);
    
    cout << " Manualmente: " << endl;
    cout << "  " << endl;
    cout << "Precio de baja Pl: " << optLow[0] << endl;
    cout << "Precio de alta Ph: " << optHigh[0] << endl;
    cout << "  " << endl;
    cout << "Theta de baja Thetal: " << optLow[1] << endl;
    cout << "Theta de alta Thetah " << optHigh[1] << endl;
    cout << "  " << endl;
    cout << "------------------------------------------------- " << endl;
    cout << "  " << endl;
    
    
    F=ZFP(aalpha1, aalpha2, RMCH, RMCL, bbgrid, ttgrid, ssigma, bbeta, w, PL, PH, ggamma, Abar, r, TTHETAL, TTHETAH, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    
    double errorZ;
    double ZL;
    double ZH;
    double ZL2;
    double ZH2;
    
    ZL = F[0];
    ZH = F[0] + F[1];
    errorZ = F[2];
    
    F = Zs(aalpha1, aalpha2, RMCH, RMCL, ZH, ZL, bbgrid, ttgrid, ssigma, bbeta, w, PL, PH, ggamma, Abar, r, TTHETAL, TTHETAH, tthetaRthres, bRthreshold, Rsubsidy, tax,ZHPAR,ZLPAR, T, S);
    
    ZH2 = F[0];
    ZL2 = F[1];
    
    cout << "  " << endl;
    cout << " QUALITIES ARE: " << endl;
    cout << "  " << endl;
    cout << "Calidad Z de universidad baja: " << ZL << endl;
    cout << "Calidad Z de universidad ALTA: " << ZH << endl;
    cout << "  " << endl;
    cout << "Error es: " << errorZ << endl;
    cout << "  " << endl;
    
    cout << "  " << endl;
    cout << "Calidad Z de universidad baja - usando Zs: " << ZL2 << endl;
    cout << "Calidad Z de universidad ALTA - usando Zs: " << ZH2 << endl;
    cout << "  " << endl;
    
    vector<double> Q;
    Q.resize(2);
    Q = QUALITIES(aalpha1, aalpha2, RMCH, RMCL, bbgrid, ttgrid, ssigma, bbeta, w, PL, PH, ZL, ZH, ggamma, Abar, r, TTHETAL, TTHETAH, tthetaRthres, bRthreshold, Rsubsidy, tax, ZHPAR, ZLPAR, T, S);
    
    double Ql;
    double Qh;
    
    Ql = Q[1];
    Qh = Q[0];
    
    cout << "  " << endl;
    cout << "Calidad Z de universidad baja - calculado manualmente: " << Ql << endl;
    cout << "Calidad Z de universidad ALTA - calculado manualmente: " << Qh << endl;
    cout << "  " << endl;
    
    vector<double> Consumers;
    Consumers.resize(6);
    Consumers=aaverages(ZH, ZL, bbgrid, ttgrid, ssigma, bbeta, w, PL, PH, ggamma, Abar, r, TTHETAL, TTHETAH, tthetaRthres, bRthreshold, Rsubsidy, tax, T, S);
    
    cout << "  " << endl;
    cout << " OTHER STATISTICS ARE: " << endl;
    cout << "  " << endl;
    cout << "Chiquillos en la baja: " << Consumers[5] << " (" << Consumers[5]/(bgridsize*ttgridsize) << ")" << endl;
    cout << "Theta promedio de chiquillos en la baja: " << Consumers[3] << endl;
    cout << "Riqueza promedio de chiquillos en la baja: " << Consumers[4] << endl;
    cout << "  " << endl;
    cout << "Chiquillos en la alta: " << Consumers[2] << " (" << Consumers[2]/(bgridsize*ttgridsize) << ")"  << endl;
    cout << "Theta promedio de quicos en la alta: " << Consumers[0] << endl;
    cout << "Riqueza promedio de quicos en la alta: " << Consumers[1] << endl;
    cout << "  " << endl;
    
    
    
    if(1==2){
    cout << "------------------------------------ "<< endl;
    cout << " Start solving NE via Manual inputs " << endl;
    cout << " ---------------------------------- " << endl;
    vector<double>NASEQSOL=NASHEQUILIBRIUM(aalpha1, aalpha2,RMCH, RMCL, bbgrid, ttgrid,ssigma, bbeta, w, pL, pH, ggamma, Abar, r, tthetaLthres, tthetaHthres, tthetaRthres, bRthreshold, Rsubsidy, tax, ZHPAR, ZLPAR,1, T, S);
    
    cout << NASEQSOL[0] << " NASEQSOL0"<< endl;
    cout << NASEQSOL[1] << " NASEQSOL1"<< endl;
    cout << NASEQSOL[2] << " NASEQSOL2"<< endl;
    cout << NASEQSOL[3] << " NASEQSOL3"<< endl;
    }
    if(1==2){
    }
    return 0;
}

