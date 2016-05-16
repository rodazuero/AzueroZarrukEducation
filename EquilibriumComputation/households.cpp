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
                 double tax){
  // Relevant Abar
  Abar = AbarC(ttheta, tthetaRthreshold, Abar);
  
  //Computing subsidized interest rate
  
  double subR=rsub(ttheta, tthetaRthreshold, b, bRthreshold, Rsubsidy);
  
  //Computing thefinal income
  b=bfinal(b, tax);
  
  
  
  //If person's ttheta is below the threshold of admission, utility of -10000
  if (ttheta<tthetaadmission){
    cout << "pailas" << endl;
    return(-10000);
  }
  else{
    //Consumption in t=0
    double num0=bbeta*(1+r);
    num0=pow(num0,-1/ssigma);
    
    double num1=w*ttheta+b*(1+r)-P*(1+r)
      +min(P,Abar)*r*subR;
    
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
              double tax){
  
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
    
    //Consumption in t=0
    double cons0=b-P+Abar;
    
    
    //Consumption in t=1
    double cons1=w*ttheta-(1+r)*Abar+r*subR*min(Abar,P);
    
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
                   double tax){
  
  // Relevant Abar
  Abar = 0;
  
  //Computing final income
  b=bfinal(b,tax);
  
  //Consumption in t=0
  double num0=bbeta*(1+r);
  num0=pow(num0,-1/ssigma);
  
  double num1=w*ttheta*(2+r)+b*(1+r);
  
  double den0=bbeta*(1+r);
  den0=pow(den0,-1/ssigma);
  den0=(1+r)*den0;
  den0=1+den0;
  
  double cons0=num0*num1/den0;
  
  
  //Consumption in t=1
  double cons1=num1/den0;
  
  //Computing utility if studies
  
  double utility=pow(cons0,1-ssigma)/(1-ssigma)+bbeta*pow(cons1,1-ssigma)/(1-ssigma);
  
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
                 double tax){
  
  // Relevant Abar
  Abar = 0;
  
  //Final income
  b=bfinal(b, tax);
  
  //Consumption in t=0
  double cons0=b+Abar+w;
  
  
  //Consumption in t=1
  double cons1=w*ttheta-(1+r)*Abar;
  
  //Computing utility:
  double utility=0;
  
  
  //Consumption negative,bad utility
  if (cons0 <0 || cons1<0){
    utility=pow(-10.0,5.0);
  }
  
  
  //Consumption positive, compute usual utility
  else{
    utility=pow(cons0,1-ssigma)/(1-ssigma)+
      bbeta*pow(cons1,1-ssigma)/(1-ssigma);
  }
  
  return(utility);
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
                         double tax){
  
  vector<double> decision;
  decision.resize(2);
  
  double uncstudy_h = UNCSTUDY(b, ttheta, ssigma, bbeta, wh, Ph, ggamma, r, Abar, tthetaadmissionh, tthetaRthreshold, bRthreshold, Rsubsidy, tax);
  double uncstudy_l = UNCSTUDY(b, ttheta, ssigma, bbeta, wl, Pl, ggamma, r, Abar, tthetaadmissionl, tthetaRthreshold, bRthreshold, Rsubsidy, tax);
  
  double cstudy_h = CSTUDY(b, ttheta, ssigma, bbeta, wh, Ph, ggamma, r, Abar, tthetaadmissionh, tthetaRthreshold, bRthreshold, Rsubsidy, tax);
  double cstudy_l = CSTUDY(b, ttheta, ssigma, bbeta, wl, Pl, ggamma, r, Abar, tthetaadmissionl, tthetaRthreshold, bRthreshold, Rsubsidy, tax);
  
  double uncnotstudy = UNCNOTSTUDY(b, ttheta, ssigma, bbeta, w, ggamma, r, Abar, tax);
  double cnotstudy = CNOTSTUDY(b, ttheta, ssigma, bbeta, w, ggamma, r, Abar, tax);
  

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
                         double tax){
  
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
  res.resize(2);
#pragma omp parallel for shared(MATGRIDh, MATGRIDl) private(res)
  for (int bb=0; bb<bbnumber; bb=bb+1){
    for (int tt=0; tt<ttnumber; tt=tt+1){
      res = decision(bbgrid[bb], ttgrid[tt], ssigma, bbeta, w, w*(1+zl), w*(1+zh), Pl, Ph, tthetaadmissionl, tthetaadmissionh, ggamma, r, Abar, tthetaRthreshold, bRthreshold, Rsubsidy, tax);
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