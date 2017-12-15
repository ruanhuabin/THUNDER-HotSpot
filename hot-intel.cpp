/*******************************************************************
 *       Filename:  hot-intel.cpp                                     
 *                                                                 
 *    Description:                                        
 *                                                                 
 *        Version:  1.0                                            
 *        Created:  12/06/2017 06:30:05 PM                                 
 *       Revision:  none                                           
 *       Compiler:  gcc                                           
 *                                                                 
 *         Author:  Ruan Huabin                                      
 *          Email:  ruanhuabin@tsinghua.edu.cn                                        
 *        Company:  Dep. of CS, Tsinghua Unversity                                      
 *                                                                 
 *******************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <immintrin.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <stdarg.h>
#include <unistd.h>
#include <memory>
#include "Macro.h"
#include "Typedef.h"
#include <sys/time.h>
using namespace std;

#define ABS(a) gsl_complex_abs(a)
#define ABS2(a) gsl_complex_abs2(a)
typedef gsl_complex Complex;

typedef struct  _complex{
    /*
     *double dat[2] __attribute__((aligned(32)));
     */
    double dat[2];
} Complex2;

inline Complex operator*(const Complex a, const double x)
{
    return gsl_complex_mul_real(a, x);
}

inline Complex operator*(const double x, const Complex a)
{
    return a * x;
}

inline Complex operator+(const Complex a, const Complex b)
{
    return gsl_complex_add(a, b);
}

inline Complex operator-(const Complex a, const Complex b)
{
    return gsl_complex_sub(a, b);
}


/*
*m: number of input graphs
*n: pix number of every input graph
*/
vec logDataVSPrior_OrigScanningHot_Oldest(const Complex* dat,
                   const Complex* pri,
                   const double* ctf,
                   const double* sigRcp,
                   const int n,
                   const int m)
{
    vec result = vec::Zero(n);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
        {
            int idx = i * n + j;
            result(j) += ABS2(dat[idx] - ctf[idx] * pri[i]) * sigRcp[idx];
        }

    return result;
}

vec logDataVSPrior_OrigScanningHot(const Complex* dat,
                   const Complex* pri,
                   const double* ctf,
                   const double* sigRcp,
                   const int n,
                   const int m)
{
    vec result = vec::Zero(n);

    Complex temp;
    double temp1;
    double temp2;
    for (int i = 0; i < m; i++)
    {   
        for (int j = 0; j < n; j++)
        {
            int idx = i * n + j;
            /*
             *result(j) += ABS2(dat[idx] - ctf[idx] * pri[i]) * sigRcp[idx];
             */
            temp = dat[idx] - ctf[idx] * pri[i];
            temp1 = ABS2(temp);
            temp2 = temp1 * sigRcp[idx];
            result(j) += temp2;
        }
    }
    return result;
}


double logDataVSPrior_OrigSearchHot(const Complex* dat, const Complex* pri, const double* ctf, const double* sigRcp, const int m)
{
    
    Complex temp;
    double temp1;
    double temp2;
    double result2 = 0.0;
    for (int i = 0; i < m; i++)
    {
        temp = dat[i] - ctf[i] * pri[i];
        temp1 = ABS2(temp);
        temp2 = temp1 * sigRcp[i];
        result2 += temp2;
    }
    return result2;
}

double logDataVSPrior_OptSearchHot(Complex2* dat, const Complex2* pri, const double* ctf, const double* sigRcp, const int m)
{


    double result2 = 0.0;
    double tmpReal = 0.0;
    double tmpImag = 0.0;
    double tmp1Real = 0.0;
    double tmp1Imag = 0.0;
    double tmp2;
    double tmp3;
    for (int i = 0; i < m; i++)
    {
            
        tmpReal = ctf[i] * pri[i].dat[0];
        tmpImag = ctf[i] * pri[i].dat[1];
        tmp1Real = dat[i].dat[0] - tmpReal;
        tmp1Imag = dat[i].dat[1] - tmpImag;

        tmp2 = tmp1Real * tmp1Real + tmp1Imag * tmp1Imag;
        tmp3 = tmp2 * sigRcp[i];
        
        result2 += tmp3;

    }



    return result2;
}




double logDataVSPrior_SIMDOptSearchHot(Complex2* dat, const Complex2* pri, const double* ctf, const double* sigRcp, const int m)
{
    __m256d ymm3, ymm4, ymm5,ymm6,ymm7,ymm8,ymm9,ymm10,ymm11,ymm12,ymm13,ymm14,ymm15,ymm16,ymm17;
    ymm17 = _mm256_setzero_pd();
    for(int i = 0; i <= (m - 4); i +=4)
    {
        ymm3 = _mm256_set_pd(ctf[i+3], ctf[i+2], ctf[i+1], ctf[i]); //ctf[i]
        ymm4 = _mm256_set_pd(dat[i+3].dat[0], dat[i+2].dat[0], dat[i+1].dat[0],dat[i].dat[0]);//dat[i].dat[0]
        ymm5 = _mm256_set_pd(dat[i+3].dat[1], dat[i+2].dat[1], dat[i+1].dat[1],dat[i].dat[1]);//dat[i].dat[1]
        ymm6 = _mm256_set_pd(pri[i+3].dat[0], pri[i+2].dat[0], pri[i+1].dat[0],pri[i].dat[0]);//pri[i].dat[0]
        ymm7 = _mm256_set_pd(pri[i+3].dat[1], pri[i+2].dat[1], pri[i+1].dat[1],pri[i].dat[1]);//pri[i].dat[1]

        /*
         *ymm3 = _mm256_loadu_pd(&ctf[i]);
         *ymm4 = _mm256_loadu_pd(&(dat[i].dat[0]));
         *ymm5 = _mm256_loadu_pd(&(dat[i].dat[1]));
         *ymm6 = _mm256_loadu_pd(&(pri[i].dat[0]));
         *ymm7 = _mm256_loadu_pd(&(pri[i].dat[1]));
         */

        ymm8 = _mm256_mul_pd(ymm3, ymm6); //tmpReal
        ymm9 = _mm256_mul_pd(ymm3, ymm7);//tmpImag
        /*
         *ymm10 = _mm256_sub_pd(ymm6, ymm8);//tmp1Real
         *ymm11 = _mm256_sub_pd(ymm7, ymm9); //tmp1Imag
         */
        
        ymm10 = _mm256_sub_pd(ymm4, ymm8);//tmp1Real
        ymm11 = _mm256_sub_pd(ymm5, ymm9); //tmp1Imag

        ymm12 = _mm256_mul_pd(ymm10, ymm10);
        ymm13 = _mm256_mul_pd(ymm11, ymm11);

        ymm14 = _mm256_add_pd(ymm12, ymm13); //tmp2

        ymm15 = _mm256_set_pd(sigRcp[i+3], sigRcp[i+2], sigRcp[i+1], sigRcp[i]); //sigRcp

        ymm16 = _mm256_mul_pd(ymm14, ymm15);//tmp3

        ymm17 = _mm256_add_pd(ymm17, ymm16);//result2

    }

    double tmp[4] __attribute__((aligned(32)));
    _mm256_store_pd(tmp, ymm17);

    double result2 = tmp[0] + tmp[1] + tmp[2] + tmp[3];
    return result2;

}

double logDataVSPrior_2SIMDOptSearchHot(Complex2* dat, const Complex2* pri, const double* ctf, const double* sigRcp, const int m)
{
    __m256d ymm3, ymm4, ymm5,ymm6,ymm7,ymm8,ymm9,ymm10,ymm11,ymm12,ymm13,ymm14,ymm15,ymm16,ymm17;
    ymm17 = _mm256_setzero_pd();
    for(int i = 0; i <= (m - 4); i +=4)
    {
/*
 *        ymm3 = _mm256_set_pd(ctf[i+3], ctf[i+2], ctf[i+1], ctf[i]); //ctf[i]
 *        ymm4 = _mm256_set_pd(dat[i+3].dat[0], dat[i+2].dat[0], dat[i+1].dat[0],dat[i].dat[0]);//dat[i].dat[0]
 *        ymm5 = _mm256_set_pd(dat[i+3].dat[1], dat[i+2].dat[1], dat[i+1].dat[1],dat[i].dat[1]);//dat[i].dat[1]
 *        ymm6 = _mm256_set_pd(pri[i+3].dat[0], pri[i+2].dat[0], pri[i+1].dat[0],pri[i].dat[0]);//pri[i].dat[0]
 *        ymm7 = _mm256_set_pd(pri[i+3].dat[1], pri[i+2].dat[1], pri[i+1].dat[1],pri[i].dat[1]);//pri[i].dat[1]
 *
 */
        ymm3 = _mm256_loadu_pd(&ctf[i]);
        ymm4 = _mm256_loadu_pd(&(dat[i].dat[0]));
        ymm5 = _mm256_loadu_pd(&(dat[i].dat[1]));
        ymm6 = _mm256_loadu_pd(&(pri[i].dat[0]));
        ymm7 = _mm256_loadu_pd(&(pri[i].dat[1]));

        ymm8 = _mm256_mul_pd(ymm3, ymm6); //tmpReal
        ymm9 = _mm256_mul_pd(ymm3, ymm7);//tmpImag
        /*
         *ymm10 = _mm256_sub_pd(ymm6, ymm8);//tmp1Real
         *ymm11 = _mm256_sub_pd(ymm7, ymm9); //tmp1Imag
         */
        
        ymm10 = _mm256_sub_pd(ymm4, ymm8);//tmp1Real
        ymm11 = _mm256_sub_pd(ymm5, ymm9); //tmp1Imag

        ymm12 = _mm256_mul_pd(ymm10, ymm10);
        ymm13 = _mm256_mul_pd(ymm11, ymm11);

        ymm14 = _mm256_add_pd(ymm12, ymm13); //tmp2

        ymm15 = _mm256_set_pd(sigRcp[i+3], sigRcp[i+2], sigRcp[i+1], sigRcp[i]); //sigRcp

        ymm16 = _mm256_mul_pd(ymm14, ymm15);//tmp3

        ymm17 = _mm256_add_pd(ymm17, ymm16);//result2

    }

    double tmp[4] __attribute__((aligned(32)));
    _mm256_store_pd(tmp, ymm17);

    double result2 = tmp[0] + tmp[1] + tmp[2] + tmp[3];
    return result2;

}


int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

int main(int argc, char *argv[])
{ 
    
    char *opt_str = "ht:f:";
    int opt = 0;
    int runTimes = 1000000000;
    int m = 14000;
    int n = 2500;
    int functionFlag = 0;
    
    while ( (opt = getopt(argc, argv, opt_str)) != -1)
    {
        switch (opt)
        {
            case 't':
                runTimes = atoi(optarg);
                break;
            case 'f':
                functionFlag = atoi(optarg);
                break;
            case 'h':
                fprintf(stderr, "%s -t <runTimes> -f <function flag>", argv[0]); 
                return 0;

            default:
                printf("Invalid command line parameters: %s!, type -h for help\n", optarg);
                return -1;
        }
    }



    Complex *dat = new Complex[m];
    Complex *pri = new Complex[m];
    Complex2 *dat2 = new Complex2[m];
    Complex2 *pri2 = new Complex2[m];
    double *ctf = new double[m];
    /*
     *__attribute__ ((aligned (64))) double *sigRcp = new double[m];
     */
    double *sigRcp = new double[m];

    for(int i = 0; i < m; i ++)
    {
        GSL_SET_COMPLEX(&dat[i], (i + 1)%100 + 0.28, (i + 2) % 100 + 0.34);
        GSL_SET_COMPLEX(&pri[i], (i + 2) % 50 + 0.59, (i + 4) % 50 + 0.78);

        dat2[i].dat[0] = (i + 1) % 100 + 0.28;
        dat2[i].dat[1] = (i + 2) % 100 + 0.34;
        pri2[i].dat[0] = (i + 2) % 50 + 0.59;
        pri2[i].dat[1] = (i + 4) % 50 + 0.78;
        ctf[i] = i % 100 + 0.78;
        sigRcp[i] = i % 50 + 0.26;
    }


    Complex *dat0 = new Complex[m * n];
    Complex *pri0 = new Complex[m];
    double *ctf0 = new double[m * n];
    double *sigRcp0 = new double[m * n];

    for(int i = 0; i < m; i ++)
    {
        GSL_SET_COMPLEX(&pri0[i], (i + 2) % 50 + 0.59, (i + 4) % 50 + 0.78);
        for(int j = 0; j < n; j ++)
        {
            int index = i * n + j;
            GSL_SET_COMPLEX(&dat0[index], (index + 1)%100 + 0.28, (index + 2) % 100 + 0.34);
            ctf0[index] = index % 100 + 0.78;
            sigRcp0[index] = index % 50 + 0.26;
        }
    }

    char *functionName[] = {"logDataVSPrior_OrigScanningHot", "logDataVSPrior_OrigSearchHot", "logDataVSPrior_OptSearchHot", "logDataVSPrior_SIMDOptSearchHot", "logDataVSPrior_2SIMDSearchHot", "logDataVSPrior_OrigScanningHot_Oldest"};
    printf("m = %d, runTimes = %ld, functionFlag = %d, functionName = %s\n", m, runTimes, functionFlag, functionName[functionFlag]);
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    double result = 0.0;

    vec resultVec;

    if(functionFlag == 0)
    {
        for(int i = 0; i < runTimes; i ++)
        {
            resultVec = logDataVSPrior_OrigScanningHot(dat0, pri0, ctf0, sigRcp0, n, m);
        }
    }
    else if(functionFlag == 1)
    {
        for(int i = 0; i < runTimes; i ++)
        {
            result = logDataVSPrior_OrigSearchHot(dat, pri, ctf, sigRcp, m);
        }
    }
    else if(functionFlag == 2)
    {
        for(int i = 0; i < runTimes; i ++)
        {
            result = logDataVSPrior_OptSearchHot(dat2, pri2, ctf, sigRcp, m);
        }
    }
    else if(functionFlag == 3)
    {
        for(int i = 0; i < runTimes; i ++)
        {    
            result = logDataVSPrior_SIMDOptSearchHot(dat2, pri2, ctf, sigRcp, m);
        }

    }
    else if(functionFlag == 4)
    {
        for(int i = 0; i < runTimes; i ++)
        {    
            result = logDataVSPrior_2SIMDOptSearchHot(dat2, pri2, ctf, sigRcp, m);
        }

    }
    else if(functionFlag == 5)
    {
        for(int i = 0; i < runTimes; i ++)
        {    
            resultVec = logDataVSPrior_OrigScanningHot_Oldest(dat0, pri0, ctf0, sigRcp0, n, m);
        }

    }


    gettimeofday(&stop, NULL);
    
    struct timeval timeElapse;
    timeval_subtract(&timeElapse, &stop, &start);
    
    float totalRunSecs = (float)timeElapse.tv_sec + (float)timeElapse.tv_usec / (float)1000000;
    if(functionFlag != 0 && functionFlag != 5)
        printf("result: %f\n\n", result);
    else
    {
        for(int i = 0; i < 5; i ++)
        {
            printf("resultVec[%d]: %f\n", i, resultVec(i));
        }
    }
    printf("Elapse Time: %f seconds\n", totalRunSecs);



    
    delete[] dat;
    delete[] pri;
    delete[] dat2;
    delete[] pri2;
    delete[] ctf;
    delete[] sigRcp;
    delete[] dat0;
    delete[] pri0;
    delete[] ctf0;
    delete[] sigRcp0;
   /*
    *free(sigRcp);
    */
    return EXIT_SUCCESS;
}
