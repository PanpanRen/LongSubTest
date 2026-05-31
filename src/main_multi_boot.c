#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define pi 3.1415926

static int MatrixInvSymmetric(double *a,int n){
	int i,j,k,m;
    double w,g,*b;
    b = (double*)malloc(n*sizeof(double));

    for (k=0; k<=n-1; k++){
        w=a[0];
        if (fabs(w)+1.0==1.0){
            free(b); return(-2);
        }
        m=n-k-1;
        for (i=1; i<=n-1; i++){
            g=a[i*n]; b[i]=g/w;
            if (i<=m) b[i]=-b[i];
            for (j=1; j<=i; j++)
                a[(i-1)*n+j-1]=a[i*n+j]+g*b[j];
        }
        a[n*n-1]=1.0/w;
        for (i=1; i<=n-1; i++)
            a[(n-1)*n+i-1]=b[i];
    }
    for (i=0; i<=n-2; i++)
        for (j=i+1; j<=n-1; j++)
            a[i*n+j]=a[j*n+i];
    free(b);
    return(2);
}

static void AbyB(double *outVector, double *A, double *v, int n, int p, int q){
	int i,j,k;
	double tmp;
	for (i=0;i<n;i++){
		for(k=0;k<q;k++){
			tmp = 0;
			for(j=0;j<p;j++)
				tmp += A[j*n + i]*v[k*p + j];
			outVector[k*n+i] = tmp;
		}
	}
}

static void ATbyB(double *outVector, double *A, double *v, int n, int p, int q){
	int i,j,k;
	double tmp;
	for (j=0;j<p;j++){
		for(k=0;k<q;k++){
			tmp = 0;
			for(i=0;i<n;i++)
				tmp += A[j*n + i]*v[k*n + i];
			outVector[k*p+j] = tmp;
		}
	}
}

static void Initialest(double *alpha,int n,int m,int p,double *x,double *y,double *hess){
    int nm = n*m;
    int i,j,k, temp, temp1, temp2;
    double tmp, *xy;
    xy 	= (double*)malloc(sizeof(double)*p);

    for(j=0; j < p; j++){
        temp = j*nm;
        temp2 = j*p;
		for(k=j; k < p; k++){
            temp1 = k*nm;
			tmp = 0.0;
			for(i=0; i<nm; i++){
				tmp += x[temp+i]*x[temp1+i];
			}
			hess[temp2+k] = tmp;
		}
	}
    for(j=1; j < p; j++){
        temp2 = j*p;
		for(k=0; k < j; k++){
			hess[temp2+k] = hess[k*p+j]; // p*p
		}
	}
    MatrixInvSymmetric(hess,p);
    for(j=0;j<p;j++) xy[j] = 0.0;
	for(i=0;i<nm;i++){
		for(j=0;j<p;j++){
			xy[j] 	+= x[j*nm+i]*y[i]; // p*1
		}
	}
    AbyB(alpha,hess,xy,p,p,1);
    free(xy);
}

// the inverse of sigma
static void sigmatilde(double *sigma,double *alpha,double *xalpha,int n,int m,int p,double *x,double *y){
    int nm = n*m;
    int i,j,k, temp, temp1,temp2;
    double tmp;

    AbyB(xalpha,x,alpha,nm,p,1);
    for(j=0;j<m;j++){
        temp = j*m;
        for(k=j;k<m;k++){
            tmp = 0.0;
            for(i=0;i<n;i++){
                temp1 = i*m+j;
                temp2 = i*m+k;
                tmp += (y[temp1]-xalpha[temp1])*(y[temp2]-xalpha[temp2]);
            }
            sigma[temp+k] = tmp/(n); // m*m
        }
    }
    for(j=1;j<m;j++){
        temp = j*m;
        for(k=0;k<j;k++){
            sigma[temp+k] = sigma[k*m+j];
        }
    }
    MatrixInvSymmetric(sigma,m);
}

static void Finalest(double *alpha,int n,int m,int p,double *x,double *y,double *sigma){
    int nm = n*m;
    int i,j,k;
    double *tem,*tem1,*tem2;
    tem = (double*)malloc(sizeof(double)*(p*nm));
    tem1 = (double*)malloc(sizeof(double)*(p*p));
    tem2 = (double*)malloc(sizeof(double)*(p));
    for(i=0;i<nm;i++){
        int res = i%m;
        int res1 = i/m;
        for(k=0;k<p;k++){
            double tmp = 0.0;
            for(j=0;j<m;j++){
                tmp += sigma[j*m+res]*x[k*nm+j+res1*m];
            }
            tem[i*p+k] = tmp;
        }
    }
    AbyB(tem1,tem,x,p,nm,p);
    MatrixInvSymmetric(tem1,p);
    AbyB(tem2,tem,y,p,nm,1);
    AbyB(alpha,tem1,tem2,p,p,1);

    free(tem);
    free(tem1);
    free(tem2);
}

// \hat{J}
void Jhat(double *J,double *sigma,double *x,int n,int m,int p){
    int nm = n*m;
    int i, j, k, l1, l2, temp1, temp2, temp_sigma;
    double tmp;

    // (1/n) * sum(X_i^T * Sigma^{-1} * X_i)
    for(j = 0; j < p; j++){
        temp1 = j*nm;
        for(k = j; k < p; k++){
            temp2 = k*nm;
            tmp = 0.0;
            for(i = 0; i < n; i++){
                for(l1 = 0; l1 < m; l1++){
                    temp_sigma = l1*m;
                    for(l2 = 0; l2 < m; l2++){
                        tmp += x[temp1 + i*m + l1] * sigma[temp_sigma + l2] * x[temp2 + i*m + l2];
                    }
                }
            }
            J[j*p + k] = tmp / n;
        }
    }
    
    for(j = 1; j < p; j++){
        for(k = 0; k < j; k++){
            J[j*p + k] = J[k*p + j];
        }
    }
    
    MatrixInvSymmetric(J, p);
}

void omega(double *omegaij, int n, int r, double *u)
{
    int i, j, k, temp;
    double *normij;
    double tmp1, tmp2, rho, coef;

    normij = (double*)malloc(sizeof(double) * n);

    for (i = 0; i < n * n; i++) omegaij[i] = 0.0;

    for(i=0;i<n;i++){
        tmp1 = 0.0;
        for(j=0;j<r;j++){
            temp = j*n+i;
            tmp1 += u[temp]*u[temp];
        }
        normij[i] = 1.0 / sqrt(tmp1);
    }

    for(i=0;i<n;i++){
        for(k=0;k<n;k++){
            if(i==k){
                omegaij[i*n+k] = 0.5;
                continue;
            }
            else{
            tmp2 = 0.0;
            for(j=0;j<r;j++){
                tmp2 += u[j*n+k]*u[j*n+i];
            }
            rho = tmp2*normij[i]*normij[k];
            coef = rho / sqrt(1.0 - rho * rho);
            omegaij[i*n+k] = 1.0/4.0+atan(coef)/(2.0*pi);
        }
        }
    }

    free(normij);
}

void Sigbyres(double *Sigbyre, double *sigma, int n, int m, int p,
                 double *x, double *y, double *alpha)
{
    int nm  = n*m, i;
    double *hess,*xalpha,*tmp;
    hess = (double*)malloc(sizeof(double)*(p*p));
    xalpha = (double*)malloc(sizeof(double)*(nm));
    tmp = (double*)malloc(sizeof(double)*(nm));

    for (i = 0; i < nm; i++) Sigbyre[i] = 0.0; // m*n

    Initialest(alpha,n,m,p,x,y,hess);
    sigmatilde(sigma,alpha,xalpha,n,m,p,x,y);
    Finalest(alpha,n,m,p,x,y,sigma);
    AbyB(xalpha,x,alpha,nm,p,1);
    for(i=0; i<nm; i++){
        tmp[i] = y[i] - xalpha[i]; // m*n
    }
    AbyB(Sigbyre,sigma,tmp,m,m,n);

    free(hess);
    free(xalpha);
    free(tmp);
}

void kernel(double *hij, double *tn, int n, int m, int p, int q, int r,
                 double *x, double *y, double *z, double *u)
{
    int nm  = n*m, i, k, j, temp, temp1, temp2, temp3, l1, l2;
    double *alpha,*sigma,*Sigbyre,*J,*omegaij,*scoreij,*score_innerij,
    *rhoij,*S0,*MT,*A,*K_mat,*Hmat,tmp1;
    alpha         = (double*)malloc(sizeof(double)*(p));
    sigma         = (double*)malloc(sizeof(double)*(m*m));
    Sigbyre       = (double*)malloc(sizeof(double)*(m*n));
    J             = (double*)malloc(sizeof(double)*(p*p));
    omegaij       = (double*)malloc(sizeof(double)*(n*n));
    scoreij       = (double*)malloc(sizeof(double)*(q*n));
    score_innerij = (double*)malloc(sizeof(double)*(n*n));
    rhoij         = (double*)malloc(sizeof(double)*(n*n));
    S0       = (double*)malloc(sizeof(double) * p * n);
    MT       = (double*)malloc(sizeof(double) * p * q * n);
    A        = (double*)malloc(sizeof(double) * p * q * n);
    K_mat    = (double*)malloc(sizeof(double) * p * n);
    Hmat     = (double*)malloc(sizeof(double) * p * p);
    // row_mean = (double*)malloc(sizeof(double) * n);
    // col_mean = (double*)malloc(sizeof(double) * n);

    for (i = 0; i < n * n; i++) hij[i] = 0.0;

    Sigbyres(Sigbyre, sigma, n, m, p, x, y, alpha);
    Jhat(J,sigma,x,n,m,p);
    omega(omegaij, n, r, u);

    int mq = m*q;
    for(i=0;i<n;i++){
        temp1 = i*mq;
        temp = i*m;
        temp3 = i*q;
        for(j=0;j<q;j++){
            temp2 = temp1 + j;
            tmp1 = 0.0;
            for(k=0;k<m;k++){
                tmp1 += z[temp2+k*q]*Sigbyre[temp + k];
            }
            scoreij[temp3+j] = tmp1;
        }
    }

    *tn = 0.0;
    for(i=0;i<n;i++){
        for(k=0;k<n;k++){
            score_innerij[i*n+k] = 0.0;
            for(j=0;j<q;j++){
                score_innerij[i*n+k] += scoreij[k*q+j]*scoreij[i*q+j];
            }
            rhoij[i*n+k] = score_innerij[i*n+k] * omegaij[i*n+k];
            if (i!=k) { *tn += rhoij[i*n+k]; }
        }
    }
    *tn /= (n * (n - 1.0));

    // S0 (X_i^T * \tilde{\Sigma}^{-1} (Y_i - X_i\hat{\alpha}))
    for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) {
            double tmp_S0 = 0.0;
            for (k = 0; k < m; k++) {
                tmp_S0 += x[j*nm + i*m + k] * Sigbyre[i*m + k];
            }
            S0[i*p + j] = tmp_S0;
        }
    }

    // M_k^T = X_k^T * \tilde{\Sigma}^{-1} * Z_k
    for (int idx = 0; idx < n; idx++) {
        for (int row_x = 0; row_x < p; row_x++) {
            for (int col_z = 0; col_z < q; col_z++) {
                double tmp_MT = 0.0;
                for (l1 = 0; l1 < m; l1++) {
                    for (l2 = 0; l2 < m; l2++) {
                        tmp_MT += x[row_x * nm + idx * m + l1] * sigma[l1 * m + l2] * z[idx * m * q + l2 * q + col_z];
                    }
                }
                MT[idx * p * q + col_z * p + row_x] = tmp_MT;
            }
        }
    }

    // A_k = \hat{J} * M_k^T
    for (k = 0; k < n; k++) {
        for (j = 0; j < q; j++) {
            for (i = 0; i < p; i++) {
                double tmp_A = 0.0;
                for (l1 = 0; l1 < p; l1++) {
                    tmp_A += J[l1 * p + i] * MT[k * p * q + j * p + l1];
                }
                A[k * p * q + j * p + i] = tmp_A;
            }
        }
    }

    // K_i (\hat{K}_i = (1/n) * \sum_{k=1}^n A_k * score * \omega_{ki})
    for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) { 
            double tmp_K = 0.0;
            for (k = 0; k < n; k++) {
                double Ak_Di_j = 0.0;
                for (l1 = 0; l1 < q; l1++) {
                    Ak_Di_j += A[k * p * q + l1 * p + j] * scoreij[i * q + l1];
                }
                tmp_K += Ak_Di_j * omegaij[k * n + i]; 
            }
            K_mat[i * p + j] = tmp_K / n;
        }
    }

    // Hmat (\hat{H} = (1/n^2) * \sum_{k=1}^n \sum_{l=1}^n A_k * A_l^T * \omega_{kl})
    for (i = 0; i < p * p; i++) Hmat[i] = 0.0;
    
    for (k = 0; k < n; k++) {
        for (int l = 0; l < n; l++) {
            double w_kl = omegaij[k * n + l];
            for (i = 0; i < p; i++) {
                for (j = 0; j < p; j++) {
                    double tmp_H = 0.0;
                    for (l1 = 0; l1 < q; l1++) {
                        tmp_H += A[k * p * q + l1 * p + i] * A[l * p * q + l1 * p + j];
                    }
                    Hmat[j * p + i] += tmp_H * w_kl;
                }
            }
        }
    }
    for (i = 0; i < p * p; i++) Hmat[i] /= (n * n);



    ////////////////////////////////////////////////////////////
    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            double term2 = 0.0, term3 = 0.0, term4 = 0.0;
            
            // term2 = E_i^T * \hat{K}_k
            // term3 = \hat{K}_i^T * E_k
            for (j = 0; j < p; j++) {
                term2 += S0[i * p + j] * K_mat[k * p + j];
                term3 += K_mat[i * p + j] * S0[k * p + j];
            }
            
            // term4 = E_i^T * \hat{H} * E_k
            for (j = 0; j < p; j++) {
                double HE_kj = 0.0;
                for (l1 = 0; l1 < p; l1++) {
                    HE_kj += Hmat[l1 * p + j] * S0[k * p + l1];
                }
                term4 += S0[i * p + j] * HE_kj;
            }
            
            hij[i * n + k] = rhoij[i * n + k] - term2 - term3 + term4;
        }
    }
    ////////////////////////////////////////////////////////////

    // double grand_mean = 0.0;
    // for (i = 0; i < n; i++) { row_mean[i] = 0.0; col_mean[i] = 0.0; }
    
    // for (i = 0; i < n; i++) {
    //     for (k = 0; k < n; k++) {
    //         double term2 = 0.0, term3 = 0.0, term4 = 0.0;
            
    //         // term2 = E_i^T * \hat{K}_k
    //         // term3 = \hat{K}_i^T * E_k
    //         for (j = 0; j < p; j++) {
    //             term2 += S0[i * p + j] * K_mat[k * p + j];
    //             term3 += K_mat[i * p + j] * S0[k * p + j];
    //         }
            
    //         // term4 = E_i^T * \hat{H} * E_k
    //         for (j = 0; j < p; j++) {
    //             double HE_kj = 0.0;
    //             for (l1 = 0; l1 < p; l1++) {
    //                 HE_kj += Hmat[l1 * p + j] * S0[k * p + l1];
    //             }
    //             term4 += S0[i * p + j] * HE_kj;
    //         }
            
    //         hij[i * n + k] = rhoij[i * n + k] - term2 - term3 + term4;
            
    //         row_mean[i] += hij[i * n + k];
    //         col_mean[k] += hij[i * n + k];
    //         grand_mean  += hij[i * n + k];
    //     }
    // }
    
    // // Double Centering
    // for (i = 0; i < n; i++) { row_mean[i] /= n; col_mean[i] /= n; }
    // grand_mean /= (n * n);
    
    // for (i = 0; i < n; i++) {
    //     for (k = 0; k < n; k++) {
    //         hij[i * n + k] = hij[i * n + k] - row_mean[i] - col_mean[k] + grand_mean;
    //     }
    // }

    free(alpha);
    free(sigma);
    free(Sigbyre);
    free(J);
    free(omegaij);
    free(scoreij);
    free(score_innerij);
    free(rhoij);

    free(S0);
    free(MT);
    free(A);
    free(K_mat);
    free(Hmat);
    // free(row_mean);
    // free(col_mean);
}

void pval(double *result, int n, int m, int p, int q, int r, double *x, double *y, double *z, double *u, int B, double *zeta)
{
    int i, k, b;
    double tn = 0.0;
    double *hij;
    int count = 0;

    hij = (double*)malloc(sizeof(double) * (n * n));

    kernel(hij, &tn, n, m, p, q, r, x, y, z, u);

    // U_n
    double h_mean = 0.0;
    for (i = 0; i < n - 1; i++) {
        for (k = i + 1; k < n; k++) {
            h_mean += hij[i * n + k];
        }
    }
    h_mean = h_mean * 2.0 / (n * (n - 1.0));

    // Bootstrap
    for (b = 0; b < B; b++) {
        double Tn_star = 0.0;
        
        double *current_zeta = zeta + b * n;

        for (i = 0; i < n - 1; i++) {
            for (k = i + 1; k < n; k++) {
                Tn_star += current_zeta[i] * current_zeta[k] * (hij[i * n + k] - h_mean);
            }
        }
        
        Tn_star = Tn_star * 2.0 / (n * (n - 1.0));

        /////////////////////////////////////////////////
        // if (Tn_star > tn) {
        //     count++;
        // }
        if (Tn_star > h_mean) {
            count++;
        }
        /////////////////////////////////////////////////
    }

    result[0] = tn;
    result[1] = (count*1.0)/(B*1.0);

    free(hij);

}

SEXP PVAL(SEXP N, SEXP M, SEXP P, SEXP Q, SEXP R, SEXP X, SEXP Y, SEXP Z, SEXP U, 
    SEXP BB, SEXP Zeta){
    
    int n = INTEGER(N)[0];
    int m = INTEGER(M)[0];
    int p = INTEGER(P)[0];
    int q = INTEGER(Q)[0];
    int r = INTEGER(R)[0];
    int B = INTEGER(BB)[0];

    double *x = REAL(X);
    double *y = REAL(Y);
    double *z = REAL(Z);
    double *u = REAL(U);
    double *zeta = REAL(Zeta);

    double result_c[2] = {0.0, 0.0};
    pval(result_c, n, m, p, q, r, x, y, z, u, B, zeta);

    SEXP tn0, pval_res, list, list_names;
    PROTECT(tn0        = allocVector(REALSXP, 1));
    PROTECT(pval_res   = allocVector(REALSXP, 1));
    PROTECT(list       = allocVector(VECSXP,  2));
    PROTECT(list_names = allocVector(STRSXP,  2));

    REAL(tn0)[0]      = result_c[0];
    REAL(pval_res)[0] = result_c[1];

    SET_STRING_ELT(list_names, 0, mkChar("tn0"));
    SET_STRING_ELT(list_names, 1, mkChar("pval"));

    SET_VECTOR_ELT(list, 0, tn0);
    SET_VECTOR_ELT(list, 1, pval_res);

    setAttrib(list, R_NamesSymbol, list_names);
    
    UNPROTECT(4);
    return list;
}