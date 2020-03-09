#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h> 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include "nmath.h"
//#include "functions.h"

#define repeat for(;;)

//extern "C" {
	
double rbinom(double nin, double pp)
{
    /* FIXME: These should become THREAD_specific globals : */
	
    static double c, fm, npq, p1, p2, p3, p4, qn;
    static double xl, xll, xlr, xm, xr;
	
    static double psave = -1.0;
    static int nsave = -1;
    static int m;
	
    double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
    double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
    int i, ix, k, n;
	
    if (!R_FINITE(nin)) ML_ERR_return_NAN;
    r = floor(nin + 0.5);
    if (r != nin) ML_ERR_return_NAN;
    if (!R_FINITE(pp) ||
		/* n=0, p=0, p=1 are not errors <TSL>*/
		r < 0 || pp < 0. || pp > 1.)	ML_ERR_return_NAN;
	
    if (r == 0 || pp == 0.) return 0;
    if (pp == 1.) return r;
	
    if (r >= INT_MAX)/* evade integer overflow,
                      and r == INT_MAX gave only even values */
		return qbinom(unif_rand(), r, pp, /*lower_tail*/ 0, /*log_p*/ 0);
    /* else */
    n = (int) r;
	
    p = fmin2(pp, 1. - pp);
    q = 1. - p;
    np = n * p;
    r = p / q;
    g = r * (n + 1);
	
    /* Setup, perform only when parameters change [using static (globals): */
	
    /* FIXING: Want this thread safe
	 -- use as little (thread globals) as possible
	 */
    if (pp != psave || n != nsave) {
		psave = pp;
		nsave = n;
		if (np < 30.0) {
			/* inverse cdf logic for mean less than 30 */
			qn = pow(q, (double) n);
			goto L_np_small;
		} else {
			ffm = np + p;
			m = (int) ffm;
			fm = m;
			npq = np * q;
			p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
			xm = fm + 0.5;
			xl = xm - p1;
			xr = xm + p1;
			c = 0.134 + 20.5 / (15.3 + fm);
			al = (ffm - xl) / (ffm - xl * p);
			xll = al * (1.0 + 0.5 * al);
			al = (xr - ffm) / (xr * q);
			xlr = al * (1.0 + 0.5 * al);
			p2 = p1 * (1.0 + c + c);
			p3 = p2 + c / xll;
			p4 = p3 + c / xlr;
		}
    } else if (n == nsave) {
		if (np < 30.0)
			goto L_np_small;
    }
	
    /*-------------------------- np = n*p >= 30 : ------------------- */
    repeat {
		u = unif_rand() * p4;
		v = unif_rand();
		/* triangular region */
		if (u <= p1) {
			ix = (int)(xm - p1 * v + u);
			goto finis;
		}
		/* parallelogram region */
		if (u <= p2) {
			x = xl + (u - p1) / c;
			v = v * c + 1.0 - fabs(xm - x) / p1;
			if (v > 1.0 || v <= 0.)
				continue;
			ix = (int) x;
		} else {
			if (u > p3) {	/* right tail */
				ix = (int)(xr - log(v) / xlr);
				if (ix > n)
					continue;
				v = v * (u - p3) * xlr;
			} else {/* left tail */
				ix = (int)(xl + log(v) / xll);
				if (ix < 0)
					continue;
				v = v * (u - p2) * xll;
			}
		}
		/* determine appropriate way to perform accept/reject test */
		k = abs(ix - m);
		if (k <= 20 || k >= npq / 2 - 1) {
			/* explicit evaluation */
			f = 1.0;
			if (m < ix) {
				for (i = m + 1; i <= ix; i++)
					f *= (g / i - r);
			} else if (m != ix) {
				for (i = ix + 1; i <= m; i++)
					f /= (g / i - r);
			}
			if (v <= f)
				goto finis;
		} else {
			/* squeezing using upper and lower bounds on log(f(x)) */
			amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
			ynorm = -k * k / (2.0 * npq);
			alv = log(v);
			if (alv < ynorm - amaxp)
				goto finis;
			if (alv <= ynorm + amaxp) {
				/* stirling's formula to machine accuracy */
				/* for the final acceptance/rejection test */
				x1 = ix + 1;
				f1 = fm + 1.0;
				z = n + 1 - fm;
				w = n - ix + 1.0;
				z2 = z * z;
				x2 = x1 * x1;
				f2 = f1 * f1;
				w2 = w * w;
				if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
					goto finis;
			}
		}
	}
	
L_np_small:
    /*---------------------- np = n*p < 30 : ------------------------- */
	
	repeat {
		ix = 0;
		f = qn;
		u = unif_rand();
		repeat {
			if (u < f)
				goto finis;
			if (ix > 110)
				break;
			u -= f;
			ix++;
			f *= (g / ix - r);
		}
	}
finis:
    if (psave > 0.5)
		ix = n - ix;
	return (double)ix;
}

double* drift(int* Ne, double* P_freq, int* length)
{
	int i;
	 
    double N = *Ne;  //this needs to be double for division
	//printf("after drift /n");
	for (i=0; i<*length; i++) {
        if (N < 1) {
            P_freq[i]=0;
            printf("N is less than 1");
        }
        if (N > 0){
			GetRNGstate();
            P_freq[i] = rbinom(*Ne, P_freq[i]) / N;
			PutRNGstate();
		//	printf("%d %f %f\n", *Ne, N, P_freq[i]);
        }    
	}
	return P_freq;
}

double* mutation(double* u, double* P_freq, int* length)
{
	int i;
	for (i=0; i<*length; i++) {
		P_freq[i] = (1-*u) * P_freq[i] + *u * (1-P_freq[i]);
	}	
    return P_freq;
}


double* migration_IM(double* m, double* P_freq, int* length)
{
	int i;
	double mean_p;
	double sum=0;
	for (i=0; i<*length; i++){
		sum = sum + P_freq[i];
	}
	mean_p = sum / *length;
	//printf("after migration /n");
	//printf("%f\n", mean_p);
    
	for (i=0; i<*length; i++){
		P_freq[i] = P_freq[i] * (1- *m) + *m * (mean_p);
	}	
    
    return P_freq;
}

double* migrationHierIM(double* P_freq, double* m_demes, double* m_groups,  int* totNumDemes, int* numGroups, int* numDemesPerGroup)
{
	//m_demes is migration among demes within groups
	//m_groups is migration among groups
	int i;
	double mean_p_Overall;
	double mean_p_WithinGroup[*numGroups];
	int GroupID=0;
	double sum=0;
	int j=0;
	int gr;
	double meanp_thisgroup;
	
	
	for (i=0; i<*totNumDemes; i++){
		sum = sum + P_freq[i];
	}
	mean_p_Overall = sum / *totNumDemes;
	//printf("-------------------yo------------------------ \n");
	
	sum=0;

	for (i=0; i<(*totNumDemes); i++){
		gr = (i+1) % *numDemesPerGroup;
		
		sum = sum + P_freq[i];
		//printf("%d %d %d %f %f %d \n", i, j, gr, sum, mean_p_WithinGroup[j], *numDemesPerGroup);
		
		if (gr==0 && i>0){
			mean_p_WithinGroup[j] = sum/ *numDemesPerGroup;
		//	printf("%d %d %d %f %f %d \n", i, j, gr, sum, mean_p_WithinGroup[j], *numDemesPerGroup);
			j= j + 1;
		}
		if (gr==0){
			sum=0;
		}
	}
	
	/*for (j=0; j< *numGroups; j++){
		printf("Loop through j %f \n", mean_p_WithinGroup[j]);
	}*/
		 
	j=0;
	
	for (i=0; i<*totNumDemes; i++){
		gr = (i+1) % *numDemesPerGroup;
		meanp_thisgroup = mean_p_WithinGroup[j];
		
		
		P_freq[i] = *m_groups * (mean_p_Overall * *numGroups - meanp_thisgroup)/ (*numGroups-1) + *m_demes * (meanp_thisgroup * *numDemesPerGroup - P_freq[i]) / (*numDemesPerGroup-1) + P_freq[i]*(1- *m_demes- *m_groups);

		//P_freq[i] = *m_groups * (mean_p_Overall) + *m_demes * (meanp_thisgroup) + P_freq[i]*(1- *m_demes- *m_groups);
		//#switched m_demes and m_groups?
		if (gr==0 && i>0){
			//printf("Loop through i %d %f \n", j, meanp_thisgroup);
			j= j + 1;
		}
		
	}
    
    return P_freq;
}

void SimulateIM(double* P_freq, int* length, double* u, double* m, int* Ne, int* gens)
{
	int i;
	for (i=0; i < *gens; i++){
	  	P_freq = mutation(u, P_freq, length);
        P_freq = migration_IM(m, P_freq, length);
		P_freq = drift(Ne, P_freq, length);
	}

}

void SimulateHierIM(double* P_freq, int* totNumDemes, int* numGroups, int* numDemesPerGroup, double* u, double* m_demes, double* m_groups, int* Ne, int* gens)
//m_demes is migration between/among groups
//m_groups is migration within groups
{
	int i;
	for (i=0; i < *gens; i++){
	  	P_freq = mutation(u, P_freq, totNumDemes);
        P_freq = migrationHierIM(P_freq, m_demes, m_groups, totNumDemes, numGroups, numDemesPerGroup);
		P_freq = drift(Ne, P_freq, totNumDemes);
		//printf("%d \n", i);
	}
	
}



void thetacal(int* gen, int* psum, int* length_gen, int* noall, double* sample_size, double* no_of_samples, double* het0, double* het1, double* fst)
//int *gen[],noall,sample_size[],no_of_samples;
//float *het0,*het1,*fst;
{
	int i,j, ptot, row, col;
	double xx,yy,nbar,nc,q2,q3,nbar2;
	
	//psum = (int *)malloc(*noall*sizeof(int));
	
	//for(i=0;i< *noall;++i)psum[i] = 0;
	
	nbar = nbar2 = 0;
	
	/* for(j=0,xx=0.0;j< *no_of_samples;++j) {
		nbar += sample_size[j];
		nbar2 += (double)sample_size[j]*(double)sample_size[j];
		for(i=0;i< *noall;++i) {
			psum[i] += gen[j][i];
			if(sample_size[j] > 0)
				xx += ((double)gen[j][i] * (double)gen[j][i])/(double)sample_size[j];
		}
	}*/
	
	for (j=0; j< *no_of_samples; ++j){
		nbar += sample_size[j];
		nbar2 += (double)sample_size[j]*(double)sample_size[j];
	}
	
	for (j=0, xx=0.0; j< *length_gen; ++j){
		row=  floor(j / *no_of_samples);
		col=   j - row* *no_of_samples;
		//psum[i] += gen[j];  //this is a column sum of allele counts in all populations; input from R
		//printf("%d %d %d", row, col, psum);
		if(sample_size[row] > 0)
			xx += ((double)gen[j] * (double)gen[j])/(double)sample_size[col];
		printf("%d %d %d %d %f %f \n", j, row, col, gen[j], sample_size[col], xx);
	}
	
	nc = 1.0/(*no_of_samples - 1.0)*(nbar - nbar2/nbar);
	nbar /= *no_of_samples;
	
	for(i=0,yy=0.0; i<*noall; ++i) {
		yy += psum[i]*psum[i];
		//ptot += psum[i];
	}
	q2 = (xx- *no_of_samples)/(*no_of_samples*(nbar - 1.0));
	q3 = 1.0/(*no_of_samples*(*no_of_samples-1.0)*nbar*nc)*
	(yy - nbar*(nc-1.0)/(nbar-1.0)*xx) +
	(nbar-nc)/(nc*(nbar-1.0))*(1.0-1.0/
							   (*no_of_samples - 1.0)*xx);
	
	
	*het0 = 1.0 - q2;
	*het1 = 1.0 - q3;
	if(*het1 < 1.0e-10)*fst = -100.0;
	else *fst = 1.0 - (*het0)/(*het1);
	printf("%f %f %f %f \n", nc, nbar, nbar2, *no_of_samples);
	printf("%f %f\n", xx, yy);
	printf("%f %f\n", q2, q3);
	printf("%f %f %f \n",  *het0, *het1, *fst);
		   
	//free(psum);
	
}



//} // extern "C"
