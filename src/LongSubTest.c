#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#define pi 3.1415926

void printArrayDouble(const double *arr, int n, int ncol){
    for (int i=0; i<n; i++){
        Rprintf("%f   ",arr[i]);
        if(!((i+1)%ncol)) Rprintf("\n");
    }
    printf("\n");
}

void printArrayDouble2(double **arr, int n, int ncol){
    for (int i=0; i<n; i++){
        for (int j = 0; j < ncol; j++ ) {
            Rprintf("%f   ",arr[i][j]);
        }
        Rprintf("\n");
    }
    printf("\n");
}

void printArrayDoubleInt(const int *arr, int n, int ncol){
    for (int i=0; i<n; i++){
        printf("%d   ",arr[i]);
        if(!((i+1)%ncol)) printf("\n");
    }
    printf("\n");
}

void printArrayDouble1(const double *arr, int n, int nrow, int ncol){
    if(nrow>n) fprintf(stderr, "nrow must not be greater than n!");
    for (int i=0; i<nrow; i++){
        for(int j=0; j<ncol; j++){
            printf("%f   ",arr[j*n+i]);
        }
        printf("\n");
    }
}

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
    // input t(y)
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
    // tem = (double*)malloc(sizeof(double)*nm);

    // Initialest(alpha,n,m,p,x,y,hess);
    AbyB(xalpha,x,alpha,nm,p,1);
    // for(i=0; i<nm; i++){
    //     tem[i] = y[i] - xalpha[i]; // m*n
    // }
    for(j=0;j<m;j++){
        temp = j*m;
        for(k=j;k<m;k++){
            tmp = 0.0;
            for(i=0;i<n;i++){
                temp1 = i*m+j;
                temp2 = i*m+k;
                tmp += (y[temp1]-xalpha[temp1])*(y[temp2]-xalpha[temp2]);
                // tmp += tem[temp1+j]*tem[temp1+k];
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
    // printArrayDouble(sigma,m*m,m);
    MatrixInvSymmetric(sigma,m);
    // free(tem);
}

static void Finalest(double *alpha,int n,int m,int p,double *x,double *y,double *sigma){
    int nm = n*m;
    int i,j,k;
    double *tem,*tem1,*tem2;
    // big_matrix = (double*)malloc(sizeof(double)*(nm*nm));
    tem = (double*)malloc(sizeof(double)*(p*nm));
    tem1 = (double*)malloc(sizeof(double)*(p*p));
    tem2 = (double*)malloc(sizeof(double)*(p));
    // x1 = (double*)malloc(sizeof(double)*(p*nm));
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
    } // t(X)%*%solve(sigma)
    // for (i = 0; i < nm; i++) {
    //     temp = i*nm;
    //     for (k = 0; k < nm; k++) {
    //         big_matrix[temp+k] = 0.0;
    //     }
    //     temp = i*p;
    //     for(j=0;j<p;j++){
    //         x1[temp+j] = x[j*nm+i];
    //     }
    // }
    // for (i = 0; i < n; i++) {
    //     temp = i*m;
    //     for (j = 0; j < m; j++) {
    //         temp1 = temp + (temp+j)*nm;
    //         // temp += ((temp+j)*nm);
    //         temp2 = j*m;
    //         for (k = 0; k < m; k++) {
    //             big_matrix[temp1+k] = sigma[temp2+k];
    //         }
    //     }
    // }
    // for(i=0;i<nm;i++){
    //     for(j=0;j<p;j++){
    //         x1[i*p+j] = x[j*nm+i];
    //     }
    // }
    // AbyB(tem,x1,big_matrix,p,nm,nm);
    AbyB(tem1,tem,x,p,nm,p);
    MatrixInvSymmetric(tem1,p);
    AbyB(tem2,tem,y,p,nm,1);
    AbyB(alpha,tem1,tem2,p,p,1);

    // free(big_matrix);
    free(tem);
    free(tem1);
    free(tem2);
    // free(x1);
}

double st0(double *tnb,int n,int m,int p,int q,double *x,double *y,double *z,double *ind,int B,double *alpha,double *xi){
    int nm = n*m;
    int i,j,k,b,temp,temp1,temp2,temp4;
    double *tem,*s1,*ztu,*ktilde,*KJ,*z_i,*x_i,*tem_i,*tem1_i,*S0,*KJS0,*S1star,*S1star2,*var,*s1_boot,*xalpha,*hess;
    double *tn0,*tnboot,temp3;
    // tem1 = (double*)malloc(sizeof(double)*(nm));
    tem = (double*)malloc(sizeof(double)*(nm));
    s1 = (double*)malloc(sizeof(double)*(q));
    ztu = (double*)malloc(sizeof(double)*(q*nm));
    ktilde = (double*)malloc(sizeof(double)*(q*p));
    KJ = (double*)malloc(sizeof(double)*(q*p));
    z_i = (double*)malloc(sizeof(double)*(q*m));
    x_i = (double*)malloc(sizeof(double)*(p*m));
    tem_i = (double*)malloc(sizeof(double)*(m));
    tem1_i = (double*)malloc(sizeof(double)*(m));
    // S1 = (double*)malloc(sizeof(double)*(q));
    S0 = (double*)malloc(sizeof(double)*(p));
    KJS0 = (double*)malloc(sizeof(double)*(q));
    S1star = (double*)malloc(sizeof(double)*(q));
    S1star2 = (double*)malloc(sizeof(double)*(q*q));
    var = (double*)malloc(sizeof(double)*(q*q));
    s1_boot = (double*)malloc(sizeof(double)*(q*B));
    // s1var = (double*)malloc(sizeof(double)*(q));
    // s1_boot_b = (double*)malloc(sizeof(double)*(q));
    // s1var_boot = (double*)malloc(sizeof(double)*(q));
    tn0 = (double*)malloc(sizeof(double)*(1));
    tnboot = (double*)malloc(sizeof(double)*(1));
    xalpha = (double*)malloc(sizeof(double)*(nm));
    hess = (double*)malloc(sizeof(double)*(p*p));

    for(j=0;j<q*q;j++){var[j]=0.0;}
    for(j=0;j<q*B;j++){s1_boot[j]=0.0;}

    Initialest(alpha,n,m,p,x,y,hess);
    AbyB(xalpha,x,alpha,nm,p,1);
    // for(i=0; i<nm; i++){
    //     tem1[i] = y[i] - xalpha[i]; // m*n
    // }
    for(k=0;k<m;k++){
        for(i=0;i<n;i++){
            temp = i*m+k;
            tem[temp] = (y[temp]-xalpha[temp])*ind[i];
            // tem[temp] = tem1[temp] * ind[i];
        }
    }
    
    // ztu:q*nm, z*ind
    int mq = m*q;
    for(i=0;i<n;i++){
        temp3 = ind[i]/n;
        temp = i*mq;
        for(k=0;k<m;k++){
            temp1 = temp + k*q;
            for(j=0;j<q;j++){
                temp2 = temp1 + j;
                ztu[temp2] = -z[temp2]*temp3;
            }
        }
    }
    AbyB(ktilde,ztu,x,q,nm,p);
    for(j=0;j<p;j++){
        temp = j*p;
        for(k=0;k<p;k++){
            hess[temp+k] *= (-n);
        }
    } // jtilde
    AbyB(KJ,ktilde,hess,q,p,p);
    // z_i:q*m
    // x_i:p*m
    for(i=0;i<n;i++){
        temp1 = i*mq;
        temp2 = i*m;
        for(k=0;k<m;k++){
            temp4 = temp2 + k;
            tem_i[k] = tem[temp4];
            tem1_i[k] = y[temp4] - xalpha[temp4];
            // tem1_i[k] = tem1[temp4];
            temp = k*p;
            for(j=0;j<p;j++){
                x_i[temp+j] = x[j*nm+temp4];
            }
            temp = k*q;
            temp4 = temp1 + temp;
            for(j=0;j<q;j++){
                z_i[temp+j] = z[temp4+j];
            }
        }
        //tem_i: m*1
        // temp = i*m;
        // for(k=0;k<m;k++){
        //     temp1 = temp + k;
        //     tem_i[k] = tem[temp1];
        //     tem1_i[k] = tem1[temp1];
        // }
        AbyB(s1,z_i,tem_i,q,m,1); // q*1
        AbyB(S0,x_i,tem1_i,p,m,1);
        AbyB(KJS0,KJ,S0,q,p,1);
        for(j=0;j<q;j++){
            S1star[j] = s1[j] - KJS0[j];
        }
        AbyB(S1star2,S1star,S1star,q,1,q);
        for(j=0;j<q;j++){
            temp = j*q;
            for(k=0;k<q;k++){
                temp1 = temp + k;
                var[temp1] += S1star2[temp1]; // q*q
            }
        }
        for(b=0;b<B;b++){
            temp = b*q;
            temp1 = b*n+i;
            for(j=0;j<q;j++){
                s1_boot[temp+j] += S1star[j]*xi[temp1]; // q*B
            }
        }
    }
    MatrixInvSymmetric(var,q);
    // z:q*nm
    AbyB(s1,z,tem,q,nm,1);
    // AbyB(s1var,s1,var,1,q,q);
    AbyB(KJS0, s1, var, 1, q, q);
    AbyB(tn0, KJS0, s1, 1, q, 1);
    // AbyB(tn0,s1var,s1,1,q,1);
    for(b=0;b<B;b++){
        temp = b*q;
        for(j=0;j<q;j++){
            KJS0[j] = s1_boot[temp+j];
            // s1_boot_b[j] = s1_boot[temp+j];
        }
        AbyB(s1, KJS0, var, 1, q, q);
        // AbyB(s1var_boot,s1_boot_b,var,1,q,q);
        AbyB(tnboot, s1, KJS0, 1, q, 1);
        // AbyB(tnboot,s1var_boot,s1_boot_b,1,q,1);
        tnb[b] = tnboot[0];
    }
    // free(tem1);
    free(tem);
    free(s1);
    free(ztu);
    free(ktilde);
    free(KJ);
    free(z_i);
    free(x_i);
    free(tem_i);
    free(tem1_i);
    // free(S1);
    free(S0);
    free(KJS0);
    free(S1star);
    free(S1star2);
    free(var);
    free(s1_boot);
    // free(s1var);
    // free(s1_boot_b);
    // free(s1var_boot);
    free(tnboot);
    free(xalpha);
    free(hess);
    return(tn0[0]);
}

void pval0(double *result,int n,int m,int p,int q,double *x,double *y,double *z,double *ugamma,int B,int K,double *alpha,double *xi){
    double *tn0_asym,*ind,*tnb;
    double tn0 = 0.0,tn0_k,temp1;
    int count = 0;
    int i,k,b,temp;
    tn0_asym = (double*)malloc(sizeof(double)*B);
    ind = (double*)malloc(sizeof(double)*n);
    tnb = (double*)malloc(sizeof(double)*B);
    for(b=0;b<B;b++){
        tn0_asym[b] = 0.0;
    }
    // input U%*%Gamma>=0: n*K
    for(k=0;k<K;k++){
        temp = k*n;
        for(i=0;i<n;i++){
            ind[i] = ugamma[temp+i];
        }
        tn0_k = st0(tnb,n,m,p,q,x,y,z,ind,B,alpha,xi);
        if(tn0<tn0_k){tn0 = tn0_k;}
        for(b=0;b<B;b++){
            temp1 = tnb[b];
            if(tn0_asym[b]<temp1){tn0_asym[b]=temp1;}
        }
    }
    for(b=0;b<B;b++){
        if(tn0_asym[b]>tn0){count++;}
    }
    result[0] = tn0;
    result[1] = (count*1.0)/(B*1.0);
    free(tn0_asym);
    free(ind);
    free(tnb);
}

SEXP _ST0(SEXP N, SEXP M, SEXP P, SEXP Q, SEXP X, SEXP Y, SEXP Z, SEXP UGAMMA, SEXP B, SEXP K, SEXP XI){
    double *result,*alpha;
    int p = INTEGER(P)[0];
    result = (double*)malloc(sizeof(double)*2);
    alpha = (double*)malloc(sizeof(double)*p);

    SEXP tn0, pval, list, list_names;
    PROTECT(tn0 		= allocVector(REALSXP, 	1));
	PROTECT(pval 		= allocVector(REALSXP, 	1));
	PROTECT(list 		= allocVector(VECSXP, 	2));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));

    pval0(result,INTEGER(N)[0],INTEGER(M)[0],p,INTEGER(Q)[0],REAL(X),REAL(Y),REAL(Z),REAL(UGAMMA),INTEGER(B)[0],INTEGER(K)[0],alpha,REAL(XI));
    REAL(tn0)[0] = result[0];
    REAL(pval)[0] = result[1];

    SET_STRING_ELT(list_names, 	0,	mkChar("tn0"));
	SET_STRING_ELT(list_names, 	1,	mkChar("pval"));
	SET_VECTOR_ELT(list, 		0, 	tn0);
	SET_VECTOR_ELT(list, 		1, 	pval);
	setAttrib(list, R_NamesSymbol, 	list_names);

    UNPROTECT(4);
	return list;
}

double st1(double *tnb, int n, int m, int p, int q, double *x, double *y, double *z, double *ind, int B, double *alpha, double *xi){
    int i, nm = n*m, k, j, b, temp, temp1, temp3, temp4;
    double *sigma, *xalpha, *hess, *tem1, *s1, *ztu, *xsigma, *khat, *z_i, *x_i, *tem_i, *tem2_i, *S0, *KJ, *KJS0, *S1star, *S1star2, *var, *s1_boot, *tn0, *tnboot;
    double temp2;
    sigma = (double*)malloc(sizeof(double)*(m*m));
    xalpha = (double*)malloc(sizeof(double)*(nm));
    hess = (double*)malloc(sizeof(double)*(p*p));
    tem1 = (double*)malloc(sizeof(double)*(nm));
    // tem2 = (double*)malloc(sizeof(double)*(nm));
    // tem = (double*)malloc(sizeof(double)*(nm));
    s1 = (double*)malloc(sizeof(double)*(q));
    ztu = (double*)malloc(sizeof(double)*(q*nm));
    // sigma_kron = (double*)malloc(sizeof(double)*(nm*nm));
    xsigma = (double*)malloc(sizeof(double)*(nm*p));
    khat = (double*)malloc(sizeof(double)*(q*p));
    // jhat = (double*)malloc(sizeof(double)*(p*p));
    z_i = (double*)malloc(sizeof(double)*(q*m));
    x_i = (double*)malloc(sizeof(double)*(p*m));
    tem_i = (double*)malloc(sizeof(double)*(m));
    tem2_i = (double*)malloc(sizeof(double)*(m));
    // S1 = (double*)malloc(sizeof(double)*(q));
    S0 = (double*)malloc(sizeof(double)*(p));
    KJ = (double*)malloc(sizeof(double)*(q*p));
    KJS0 = (double*)malloc(sizeof(double)*(q));
    S1star = (double*)malloc(sizeof(double)*(q));
    S1star2 = (double*)malloc(sizeof(double)*(q*q));
    var = (double*)malloc(sizeof(double)*(q*q));
    s1_boot = (double*)malloc(sizeof(double)*(q*B));
    // s1var = (double*)malloc(sizeof(double)*(q));
    tn0 = (double*)malloc(sizeof(double)*(1));
    // s1_boot_b = (double*)malloc(sizeof(double)*(q));
    // s1var_boot = (double*)malloc(sizeof(double)*(q));
    tnboot = (double*)malloc(sizeof(double)*(1));

    for(j=0;j<q*q;j++){var[j]=0.0;}
    for(j=0;j<q*B;j++){s1_boot[j]=0.0;}

    Initialest(alpha,n,m,p,x,y,hess);
    sigmatilde(sigma,alpha,xalpha,n,m,p,x,y);
    Finalest(alpha,n,m,p,x,y,sigma);
    AbyB(xalpha,x,alpha,nm,p,1);
    for(i=0; i<nm; i++){
        tem1[i] = y[i] - xalpha[i]; // m*n
    }
    // AbyB(tem2, sigma, tem1, m, m, n);
    AbyB(xalpha, sigma, tem1, m, m, n);
    for(k=0;k<m;k++){
        for(i=0;i<n;i++){
            temp = i*m+k;
            // tem[temp] = xalpha[temp] * ind[i];
            tem1[temp] = xalpha[temp] * ind[i];
        }
    }
    // ztu:q*nm, z*ind
    int mq = m*q;
    for(i=0;i<n;i++){
        temp2 = ind[i];
        temp = i*mq;
        for(k=0;k<m;k++){
            temp1 = temp + k*q;
            for(j=0;j<q;j++){
                temp3 = temp1 + j;
                ztu[temp3] = z[temp3]*temp2;
            }
        }
    }
    // kronecker
    // for(i=0;i<nm*nm;i++){sigma_kron[i]=0.0;}
    // int nm1 = (nm+1)*m;
    // for(i=0;i<n;i++){
    //     temp = nm1*i;
    //     for(k=0;k<m;k++){
    //         temp1 = temp + k*nm;
    //         temp3 = k*m;
    //         for(j=0;j<m;j++){
    //             sigma_kron[temp1+j] = sigma[temp3+j];
    //         }
    //     }
    // }
    // AbyB(xsigma,sigma_kron,x,nm,nm,p);
    for(i=0;i<nm;i++){
        int res = i%m;
        int res1 = i/m;
        for(k=0;k<p;k++){
            double tmp = 0.0;
            for(j=0;j<m;j++){
                tmp += sigma[j*m+res]*x[k*nm+j+res1*m];
            }
            xsigma[k*nm+i] = tmp;
        }
    }
    AbyB(khat,ztu,xsigma,q,nm,p);
    // ATbyB(jhat,x,xsigma,nm,p,p);
    ATbyB(hess,x,xsigma,nm,p,p);
    // MatrixInvSymmetric(jhat,p);
    MatrixInvSymmetric(hess,p);
    for(i=0;i<n;i++){
        temp1 = i*mq;
        temp3 = i*m;
        for(k=0;k<m;k++){
            temp4 = temp3 + k;
            tem_i[k] = tem1[temp4];
            tem2_i[k] = xalpha[temp4];
            temp = k*p;
            for(j=0;j<p;j++){
                x_i[temp+j] = x[j*nm+temp4];
            }
            temp = k*q;
            temp4 = temp1 + temp;
            for(j=0;j<q;j++){
                z_i[temp+j] = z[temp4+j];
            }
        }
        //tem_i: m*1
        // temp3 = i*m;
        // for(k=0;k<m;k++){
        //     temp4 = temp3 + k;
        //     // tem_i[k] = tem[temp4];
        //     tem_i[k] = tem1[temp4];
        //     tem2_i[k] = xalpha[temp4];
        // }
        // AbyB(S1,z_i,tem_i,q,m,1); // q*1
        AbyB(s1,z_i,tem_i,q,m,1);
        AbyB(S0,x_i,tem2_i,p,m,1);
        // AbyB(KJ,khat,jhat,q,p,p);
        AbyB(KJ,khat,hess,q,p,p);
        AbyB(KJS0,KJ,S0,q,p,1);
        for(j=0;j<q;j++){
            S1star[j] = s1[j] - KJS0[j];
        }
        AbyB(S1star2,S1star,S1star,q,1,q);
        for(j=0;j<q;j++){
            temp = j*q;
            for(k=0;k<q;k++){
                temp1 = temp + k;
                var[temp1] += S1star2[temp1]; // q*q
            }
        }
        for(b=0;b<B;b++){
            temp = b*q;
            temp1 = b*n+i;
            for(j=0;j<q;j++){
                s1_boot[temp+j] += S1star[j]*xi[temp1]; // q*B
            }
        }
    }
    MatrixInvSymmetric(var,q);
    // z:q*nm
    // AbyB(s1,z,tem,q,nm,1);
    AbyB(s1,z,tem1,q,nm,1);
    // AbyB(s1var,s1,var,1,q,q);
    AbyB(KJS0,s1,var,1,q,q);
    AbyB(tn0,KJS0,s1,1,q,1);
    for(b=0;b<B;b++){
        temp = b*q;
        for(j=0;j<q;j++){
            // s1_boot_b[j] = s1_boot[temp+j];
            KJS0[j] = s1_boot[temp+j];
        }
        // AbyB(s1var_boot,KJS0,var,1,q,q);
        AbyB(s1,KJS0,var,1,q,q);
        // AbyB(tnboot,s1var_boot,KJS0,1,q,1);
        AbyB(tnboot,s1,KJS0,1,q,1);
        tnb[b] = tnboot[0];
    }
    free(sigma);
    free(xalpha);
    free(hess);
    free(tem1);
    // free(tem2);
    // free(tem);
    free(s1);
    free(ztu);
    // free(sigma_kron);
    free(xsigma);
    free(khat);
    // free(jhat);
    free(z_i);
    free(x_i);
    free(tem_i);
    free(tem2_i);
    // free(S1);
    free(S0);
    free(KJ);
    free(KJS0);
    free(S1star);
    free(S1star2);
    free(var);
    free(s1_boot);
    // free(s1var);
    // free(s1_boot_b);
    // free(s1var_boot);
    free(tnboot);
    return(tn0[0]);
}

void pval1(double *result,int n,int m,int p,int q,double *x,double *y,double *z,double *ugamma,int B,int K,double *alpha,double *xi){
    double *tn0_asym,*ind,*tnb;
    double tn0 = 0.0,tn0_k,temp1;
    int count = 0;
    int i,k,b,temp;
    tn0_asym = (double*)malloc(sizeof(double)*B);
    ind = (double*)malloc(sizeof(double)*n);
    tnb = (double*)malloc(sizeof(double)*B);
    for(b=0;b<B;b++){
        tn0_asym[b] = 0.0;
    }
    // input U%*%Gamma>=0: n*K
    for(k=0;k<K;k++){
        temp = k*n;
        for(i=0;i<n;i++){
            ind[i] = ugamma[temp+i];
        }
        tn0_k = st1(tnb,n,m,p,q,x,y,z,ind,B,alpha,xi);
        if(tn0<tn0_k){tn0 = tn0_k;}
        for(b=0;b<B;b++){
            temp1 = tnb[b];
            if(tn0_asym[b]<temp1){tn0_asym[b]=temp1;}
        }
    }
    for(b=0;b<B;b++){
        if(tn0_asym[b]>tn0){count++;}
    }
    result[0] = tn0;
    result[1] = (count*1.0)/(B*1.0);
    free(tn0_asym);
    free(ind);
    free(tnb);
}

SEXP _ST1(SEXP N, SEXP M, SEXP P, SEXP Q, SEXP X, SEXP Y, SEXP Z, SEXP UGAMMA, SEXP B, SEXP K, SEXP XI){
    double *result,*alpha;
    int p = INTEGER(P)[0];
    result = (double*)malloc(sizeof(double)*2);
    alpha = (double*)malloc(sizeof(double)*p);

    SEXP tn0, pval, list, list_names;
    PROTECT(tn0 		= allocVector(REALSXP, 	1));
	PROTECT(pval 		= allocVector(REALSXP, 	1));
	PROTECT(list 		= allocVector(VECSXP, 	2));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));

    pval1(result,INTEGER(N)[0],INTEGER(M)[0],p,INTEGER(Q)[0],REAL(X),REAL(Y),REAL(Z),REAL(UGAMMA),INTEGER(B)[0],INTEGER(K)[0],alpha,REAL(XI));
    REAL(tn0)[0] = result[0];
    REAL(pval)[0] = result[1];

    SET_STRING_ELT(list_names, 	0,	mkChar("tn0"));
	SET_STRING_ELT(list_names, 	1,	mkChar("pval"));
	SET_VECTOR_ELT(list, 		0, 	tn0);
	SET_VECTOR_ELT(list, 		1, 	pval);
	setAttrib(list, R_NamesSymbol, 	list_names);

    UNPROTECT(4);
	return list;
}

double st2(int n,int m,int p, int q, int r, double *x, double *y, double *z, double *u, double *alpha){
    int nm = n*m, i, k, j, temp, temp1, temp2, temp3, temp4;
    double tn2, *hess, *xalpha, *scoreij, tmp, *normij, tmp1, tmp2;
    hess = (double*)malloc(sizeof(double)*(p*p));
    xalpha = (double*)malloc(sizeof(double)*(nm));
    // tem = (double*)malloc(sizeof(double)*(nm));
    // z_i = (double*)malloc(sizeof(double)*(q*m));
    // tem_i = (double*)malloc(sizeof(double)*(m));
    // scoreij_i = (double*)malloc(sizeof(double)*(q));
    scoreij = (double*)malloc(sizeof(double)*(q*n));
    // score = (double*)malloc(sizeof(double)*(n*n));
    normij = (double*)malloc(sizeof(double)*(n));
    // normu = (double*)malloc(sizeof(double)*(n*n));
    // utu = (double*)malloc(sizeof(double)*(n*n));

    tn2 = 0.0;
    Initialest(alpha,n,m,p,x,y,hess);
    
    AbyB(xalpha,x,alpha,nm,p,1);
    // for(i=0; i<nm; i++){
    //     tem[i] = y[i] - xalpha[i]; // m*n
    // }
    int mq = m*q;
    for(i=0;i<n;i++){
        // temp1 = i*mq;
        // for(k=0;k<m;k++){
        //     temp = k*q;
        //     temp2 = temp1 + temp;
        //     for(j=0;j<q;j++){
        //         z_i[temp+j] = z[temp2+j];
        //     }
        // }
        // //tem_i: m*1
        // temp = i*m;
        // for(k=0;k<m;k++){
        //     // tem_i[k] = tem[temp+k];
        //     temp1 = temp + k;
        //     tem_i[k] = y[temp1] - xalpha[temp1];
        // }
        // AbyB(scoreij_i,z_i,tem_i,q,m,1); // q*1
        // temp = i*q;
        // for(j=0;j<q;j++){
        //     scoreij[temp+j] = scoreij_i[j];
        // }
        temp1 = i*mq;
        temp = i*m;
        temp4 = i*q;
        for(j=0;j<q;j++){
            // scoreij_i[j] = 0.0;
            temp2 = temp1 + j;
            tmp = 0.0;
            for(k=0;k<m;k++){
                temp3 = temp + k;
                tmp += z[temp2+k*q]*(y[temp3]-xalpha[temp3]);
            }
            scoreij[temp4+j] = tmp;
        }
    }

    for(i=0;i<n;i++){
        tmp = 0.0;
        for(j=0;j<r;j++){
            temp = j*n+i;
            tmp += u[temp]*u[temp];
        }
        normij[i] = 1.0 / sqrt(tmp);
    }
    for(i=0;i<n;i++){
        temp = i*n;
        for(k=0;k<n;k++){
            temp1 = temp + k;
            // utu[temp1] = 0.0;
            tmp2 = 0.0;
            for(j=0;j<r;j++){
                // utu[temp1] += u[j*n+k]*u[j*n+i];
                tmp2 += u[j*n+k]*u[j*n+i];
            }
            tmp = tmp2*normij[i]*normij[k];
            if(i==k){tmp1=0.0;}
            else{
                tmp1 = 1.0/4.0+atan(tmp/sqrt(1.0-tmp*tmp))/(2.0*pi);
            }
            // score[temp1] = 0.0;
            tmp2 = 0.0;
            for(j=0;j<q;j++){
                // score[temp1] += scoreij[k*q+j]*scoreij[i*q+j];
                tmp2 += scoreij[k*q+j]*scoreij[i*q+j];
            }
            // tn2 += score[temp1]*tmp1;
            tn2 += tmp2*tmp1;
        }
    }
    // AbyB(normu,normij,normij,n,1,n);
    // ATbyB(score,scoreij,scoreij,q,n,n);
    // for(i=0;i<n*n;i++){
    //     // varrho[i] = utu[i]/normu[i];
    //     // tmp = utu[i]/normu[i];
    //     int res = i/n;
    //     int res1 = i%n;
    //     tmp = utu[i]*normij[res]*normij[res1];
    //     // tmp = utu[i]*normu[i];
    //     if(i%(n+1)==0){tmp1=0.0;}
    //     else{
    //         tmp1 = 1.0/4.0+atan(tmp/sqrt(1.0-tmp*tmp))/(2.0*pi);
    //     }
    //     tn2 += score[i]*tmp1;
    // }
    tn2/=(n*(n-1));
    free(hess);
    free(xalpha);
    // free(tem);
    // free(z_i);
    // free(tem_i);
    // free(scoreij_i);
    free(scoreij);
    // free(score);
    free(normij);
    // free(normu);
    // free(utu);
    return(tn2);
}

void pval2(double *result, int n, int m, int p, int q, int r, double *x, double *y, double *z, double *u, int family, int B, double *alpha, double *v){
    int nm = n*m, i, b, count = 0, temp;
    double tn2, *xalpha, *tn2_asym, temp1, *yb;
    xalpha = (double*)malloc(sizeof(double)*(nm));
    // residual = (double*)malloc(sizeof(double)*(nm));
    yb = (double*)malloc(sizeof(double)*(nm));
    tn2_asym = (double*)malloc(sizeof(double)*(B));

    tn2 = st2(n,m,p,q,r,x,y,z,u,alpha);
    AbyB(xalpha,x,alpha,nm,p,1);
    // for(i=0; i<nm; i++){
    //     residual[i] = y[i] - xalpha[i]; // m*n
    // }
    for(b=0;b<B;b++){
        // temp = b*nm;
        temp = b*n;
        if(family == 1){ // gauss
            // v: mn*B
            for(i=0; i<nm; i++){
                // yb[i] = xalpha[i] + v[temp+i]*residual[i];
                temp1 = xalpha[i];
                // yb[i] = temp1 + v[temp+i]*(y[i]-temp1);
                yb[i] = temp1 + v[temp+i/m]*(y[i]-temp1);
            }
        }
        else if(family==2){ // binomial
            // v: mn*B
            for(i=0;i<nm;i++){
                // yb[i] = v[temp+i] < 1.0/(1.0+exp(-xalpha[i]));
                yb[i] = v[temp+i/m] < 1.0/(1.0+exp(-xalpha[i]));
            }
        }
        else{
            // v: mn*B
            for(i=0; i<nm; i++){
                // yb[i] = xalpha[i] + v[temp+i]*residual[i];
                temp1 = xalpha[i];
                // yb[i] = temp1 + v[temp+i]*(y[i]-temp1);
                yb[i] = temp1 + v[temp+i/m]*(y[i]-temp1);
            }
        }
        tn2_asym[b] = st2(n,m,p,q,r,x,yb,z,u,alpha);
        if(tn2_asym[b]>tn2){count++;}
    }
    // printArrayDouble(tn2_asym,20,20);
    result[0] = tn2;
    result[1] = (count*1.0)/(B*1.0);
    free(xalpha);
    // free(residual);
    free(yb);
    free(tn2_asym);
}

SEXP _WAST0(SEXP N, SEXP M, SEXP P, SEXP Q, SEXP R, SEXP X, SEXP Y, SEXP Z, SEXP U, SEXP FAMILY, SEXP B, SEXP V){
    double *result,*alpha;
    int p = INTEGER(P)[0];
    result = (double*)malloc(sizeof(double)*2);
    alpha = (double*)malloc(sizeof(double)*p);

    SEXP tn0, pval, list, list_names;
    PROTECT(tn0 		= allocVector(REALSXP, 	1));
	PROTECT(pval 		= allocVector(REALSXP, 	1));
	PROTECT(list 		= allocVector(VECSXP, 	2));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));

    pval2(result,INTEGER(N)[0],INTEGER(M)[0],p,INTEGER(Q)[0],INTEGER(R)[0],REAL(X),REAL(Y),REAL(Z),REAL(U),INTEGER(FAMILY)[0],INTEGER(B)[0],alpha,REAL(V));
    REAL(tn0)[0] = result[0];
    REAL(pval)[0] = result[1];

    SET_STRING_ELT(list_names, 	0,	mkChar("tn0"));
	SET_STRING_ELT(list_names, 	1,	mkChar("pval"));
	SET_VECTOR_ELT(list, 		0, 	tn0);
	SET_VECTOR_ELT(list, 		1, 	pval);
	setAttrib(list, R_NamesSymbol, 	list_names);

    UNPROTECT(4);
	return list;
}

double st3(int n, int m, int p, int q, int r, double *x, double *y, double *z, double *u, double *alpha){
    int nm = n*m, i, k, j, temp, temp1, temp2, temp3;
    double *hess,*sigma,*xalpha,*tmp,*tem,*scoreij,tmp1,*normij,tmp2,tn2=0.0,tmp3;
    hess = (double*)malloc(sizeof(double)*(p*p));
    sigma = (double*)malloc(sizeof(double)*(m*m));
    xalpha = (double*)malloc(sizeof(double)*(nm));
    tmp = (double*)malloc(sizeof(double)*(nm));
    tem = (double*)malloc(sizeof(double)*(m*n));
    // z_i = (double*)malloc(sizeof(double)*(q*m));
    // tem_i = (double*)malloc(sizeof(double)*(m));
    // scoreij_i = (double*)malloc(sizeof(double)*(q));
    scoreij = (double*)malloc(sizeof(double)*(q*n));
    // score = (double*)malloc(sizeof(double)*(n*n));
    normij = (double*)malloc(sizeof(double)*(n));
    // normu = (double*)malloc(sizeof(double)*(n*n));
    // utu = (double*)malloc(sizeof(double)*(n*n));

    Initialest(alpha,n,m,p,x,y,hess);
    sigmatilde(sigma,alpha,xalpha,n,m,p,x,y);
    // printArrayDouble(sigma,m*m,m);
    // sigma[0] = 1.098901;
    // sigma[1] = -0.3296703;
    // sigma[2] = -1.77e-17;
    // sigma[3] = 8.36e-18;
    // sigma[4] = -0.3296703;
    // sigma[5] = 1.197802;
    // sigma[6] = -0.3296703;
    // sigma[7] = -1.52e-17;
    // sigma[8] = -1.38e-17;
    // sigma[9] = -0.3296703;
    // sigma[10] = 1.197802;
    // sigma[11] = -0.3296703;
    // sigma[12] = 1.15e-17;
    // sigma[13] = -1.52e-17;
    // sigma[14] = -0.3286703;
    // sigma[15] = 1.098901;
    Finalest(alpha,n,m,p,x,y,sigma);
    AbyB(xalpha,x,alpha,nm,p,1);
    for(i=0; i<nm; i++){
        tmp[i] = y[i] - xalpha[i]; // m*n
    }
    AbyB(tem,sigma,tmp,m,m,n);
    int mq = m*q;
    for(i=0;i<n;i++){
        temp1 = i*mq;
        temp = i*m;
        temp3 = i*q;
        for(j=0;j<q;j++){
            temp2 = temp1 + j;
            tmp1 = 0.0;
            for(k=0;k<m;k++){
                tmp1 += z[temp2+k*q]*tem[temp + k];
            }
            scoreij[temp3+j] = tmp1;
        }

        // temp1 = i*mq;
        // for(k=0;k<m;k++){
        //     temp = k*q;
        //     temp2 = temp1 + temp;
        //     for(j=0;j<q;j++){
        //         z_i[temp+j] = z[temp2+j];
        //     }
        // }
        // //tem_i: m*1
        // temp = i*m;
        // for(k=0;k<m;k++){
        //     tem_i[k] = tem[temp+k];
        // }
        // AbyB(scoreij_i,z_i,tem_i,q,m,1); // q*1
        // temp = i*q;
        // for(j=0;j<q;j++){
        //     scoreij[temp+j] = scoreij_i[j];
        // }
    }

    for(i=0;i<n;i++){
        tmp1 = 0.0;
        for(j=0;j<r;j++){
            temp = j*n+i;
            tmp1 += u[temp]*u[temp];
        }
        normij[i] = 1.0 / sqrt(tmp1);
    }

    for(i=0;i<n;i++){
        temp = i*n;
        for(k=0;k<n;k++){
            temp1 = temp + k;
            // utu[temp1] = 0.0;
            tmp2 = 0.0;
            for(j=0;j<r;j++){
                // utu[temp1] += u[j*n+k]*u[j*n+i];
                tmp2 += u[j*n+k]*u[j*n+i];
            }
            tmp3 = tmp2*normij[i]*normij[k];
            if(i==k){tmp1=0.0;}
            else{
                tmp1 = 1.0/4.0+atan(tmp3/sqrt(1.0-tmp3*tmp3))/(2.0*pi);
            }
            // score[temp1] = 0.0;
            tmp2 = 0.0;
            for(j=0;j<q;j++){
                // score[temp1] += scoreij[k*q+j]*scoreij[i*q+j];
                tmp2 += scoreij[k*q+j]*scoreij[i*q+j];
            }
            // tn2 += score[temp1]*tmp1;
            tn2 += tmp2*tmp1;
        }
    }

    // AbyB(normu,normij,normij,n,1,n);
    // ATbyB(score,scoreij,scoreij,q,n,n);
    // // AbyB(utu, u, u, n, r, n);
    // for(i=0;i<n;i++){
    //     temp = i*n;
    //         for(k=0;k<n;k++){
    //             temp1 = temp + k;
    //             utu[temp1] = 0.0;
    //             for(j=0;j<r;j++){
    //                 utu[temp1] += u[j*n+k]*u[j*n+i];
    //         }
    //     }
    // }
    // for(i=0;i<n*n;i++){
    //     // varrho[i] = utu[i]/normu[i];
    //     // tmp = utu[i]/normu[i];
    //     tmp1 = utu[i]*normu[i];
    //     if(i%(n+1)==0){tmp2=0.0;}
    //     else{
    //         tmp2 = 1.0/4.0+atan(tmp1/sqrt(1.0-tmp1*tmp1))/(2.0*pi);
    //     }
    //     tn2 += score[i]*tmp2;
    // }
    tn2/=(n*(n-1.0));
    free(hess);
    free(sigma);
    free(xalpha);
    free(tmp);
    free(tem);
    // free(z_i);
    // free(tem_i);
    // free(scoreij_i);
    free(scoreij);
    // free(score);
    free(normij);
    // free(normu);
    // free(utu);
    return(tn2);
}

double st3_res_boot(int n, int m, int p, int q, int r, double *x, double *y, double *z, double *u, double *alpha, double *sigma){
    int nm = n*m, i, k, j, temp, temp1, temp2, temp3;
    double *hess,*xalpha,*tmp,*tem,*scoreij,tmp1,*normij,tmp2,tn2=0.0,tmp3;
    hess = (double*)malloc(sizeof(double)*(p*p));
    // sigma = (double*)malloc(sizeof(double)*(m*m));
    xalpha = (double*)malloc(sizeof(double)*(nm));
    tmp = (double*)malloc(sizeof(double)*(nm));
    tem = (double*)malloc(sizeof(double)*(m*n));
    // z_i = (double*)malloc(sizeof(double)*(q*m));
    // tem_i = (double*)malloc(sizeof(double)*(m));
    // scoreij_i = (double*)malloc(sizeof(double)*(q));
    scoreij = (double*)malloc(sizeof(double)*(q*n));
    // score = (double*)malloc(sizeof(double)*(n*n));
    normij = (double*)malloc(sizeof(double)*(n));
    // normu = (double*)malloc(sizeof(double)*(n*n));
    // utu = (double*)malloc(sizeof(double)*(n*n));

    Initialest(alpha,n,m,p,x,y,hess);
    sigmatilde(sigma,alpha,xalpha,n,m,p,x,y);
    // printArrayDouble(sigma,m*m,m);
    // sigma[0] = 1.098901;
    // sigma[1] = -0.3296703;
    // sigma[2] = -1.77e-17;
    // sigma[3] = 8.36e-18;
    // sigma[4] = -0.3296703;
    // sigma[5] = 1.197802;
    // sigma[6] = -0.3296703;
    // sigma[7] = -1.52e-17;
    // sigma[8] = -1.38e-17;
    // sigma[9] = -0.3296703;
    // sigma[10] = 1.197802;
    // sigma[11] = -0.3296703;
    // sigma[12] = 1.15e-17;
    // sigma[13] = -1.52e-17;
    // sigma[14] = -0.3286703;
    // sigma[15] = 1.098901;
    Finalest(alpha,n,m,p,x,y,sigma);
    AbyB(xalpha,x,alpha,nm,p,1);
    for(i=0; i<nm; i++){
        tmp[i] = y[i] - xalpha[i]; // m*n
    }
    AbyB(tem,sigma,tmp,m,m,n);
    int mq = m*q;
    for(i=0;i<n;i++){
        temp1 = i*mq;
        temp = i*m;
        temp3 = i*q;
        for(j=0;j<q;j++){
            temp2 = temp1 + j;
            tmp1 = 0.0;
            for(k=0;k<m;k++){
                tmp1 += z[temp2+k*q]*tem[temp + k];
            }
            scoreij[temp3+j] = tmp1;
        }

        // temp1 = i*mq;
        // for(k=0;k<m;k++){
        //     temp = k*q;
        //     temp2 = temp1 + temp;
        //     for(j=0;j<q;j++){
        //         z_i[temp+j] = z[temp2+j];
        //     }
        // }
        // //tem_i: m*1
        // temp = i*m;
        // for(k=0;k<m;k++){
        //     tem_i[k] = tem[temp+k];
        // }
        // AbyB(scoreij_i,z_i,tem_i,q,m,1); // q*1
        // temp = i*q;
        // for(j=0;j<q;j++){
        //     scoreij[temp+j] = scoreij_i[j];
        // }
    }

    for(i=0;i<n;i++){
        tmp1 = 0.0;
        for(j=0;j<r;j++){
            temp = j*n+i;
            tmp1 += u[temp]*u[temp];
        }
        normij[i] = 1.0 / sqrt(tmp1);
    }

    for(i=0;i<n;i++){
        temp = i*n;
        for(k=0;k<n;k++){
            temp1 = temp + k;
            // utu[temp1] = 0.0;
            tmp2 = 0.0;
            for(j=0;j<r;j++){
                // utu[temp1] += u[j*n+k]*u[j*n+i];
                tmp2 += u[j*n+k]*u[j*n+i];
            }
            tmp3 = tmp2*normij[i]*normij[k];
            if(i==k){tmp1=0.0;}
            else{
                tmp1 = 1.0/4.0+atan(tmp3/sqrt(1.0-tmp3*tmp3))/(2.0*pi);
            }
            // score[temp1] = 0.0;
            tmp2 = 0.0;
            for(j=0;j<q;j++){
                // score[temp1] += scoreij[k*q+j]*scoreij[i*q+j];
                tmp2 += scoreij[k*q+j]*scoreij[i*q+j];
            }
            // tn2 += score[temp1]*tmp1;
            tn2 += tmp2*tmp1;
        }
    }

    // AbyB(normu,normij,normij,n,1,n);
    // ATbyB(score,scoreij,scoreij,q,n,n);
    // // AbyB(utu, u, u, n, r, n);
    // for(i=0;i<n;i++){
    //     temp = i*n;
    //         for(k=0;k<n;k++){
    //             temp1 = temp + k;
    //             utu[temp1] = 0.0;
    //             for(j=0;j<r;j++){
    //                 utu[temp1] += u[j*n+k]*u[j*n+i];
    //         }
    //     }
    // }
    // for(i=0;i<n*n;i++){
    //     // varrho[i] = utu[i]/normu[i];
    //     // tmp = utu[i]/normu[i];
    //     tmp1 = utu[i]*normu[i];
    //     if(i%(n+1)==0){tmp2=0.0;}
    //     else{
    //         tmp2 = 1.0/4.0+atan(tmp1/sqrt(1.0-tmp1*tmp1))/(2.0*pi);
    //     }
    //     tn2 += score[i]*tmp2;
    // }
    tn2/=(n*(n-1.0));
    free(hess);
    // free(sigma);
    free(xalpha);
    free(tmp);
    free(tem);
    // free(z_i);
    // free(tem_i);
    // free(scoreij_i);
    free(scoreij);
    // free(score);
    free(normij);
    // free(normu);
    // free(utu);
    return(tn2);
}

SEXP ST3_RES_BOOT(SEXP N, SEXP M, SEXP P, SEXP Q, SEXP R, SEXP X, SEXP Y, SEXP Z, SEXP U){
    int p = INTEGER(P)[0];
    int m = INTEGER(M)[0];
    // double *alpha, *sigma;
    // alpha = (double*)malloc(sizeof(double)*p);
    // sigma = (double*)malloc(sizeof(double)*(m*m));
    SEXP tn, alpha, sigma, list, list_names;
    PROTECT(tn = allocVector(REALSXP, 1));
    PROTECT(alpha = allocVector(REALSXP, p));
    PROTECT(sigma = allocVector(REALSXP, m*m));
    PROTECT(list = allocVector(VECSXP, 3));
    PROTECT(list_names = allocVector(STRSXP, 3));
    REAL(tn)[0] = st3_res_boot(INTEGER(N)[0],m,p,INTEGER(Q)[0],INTEGER(R)[0],REAL(X),REAL(Y),REAL(Z),REAL(U),REAL(alpha),REAL(sigma));

    SET_STRING_ELT(list_names, 	0,	mkChar("tn"));
	SET_STRING_ELT(list_names, 	1,	mkChar("alphahat"));
	SET_STRING_ELT(list_names, 	2,	mkChar("sigmainvhat"));
	SET_VECTOR_ELT(list, 		0, 	tn);
	SET_VECTOR_ELT(list, 		1, 	alpha);
	SET_VECTOR_ELT(list, 		2, 	sigma);
	setAttrib(list, R_NamesSymbol, 	list_names);

    UNPROTECT(5);
	return list;
}

void pval3(double *result, int n, int m, int p, int q, int r, double *x, double *y, double *z, double *u, int family, int B, double *alpha, double *v){
    int nm = n*m, i, b, count = 0, temp;
    double tn2, *xalpha, *yb, *tn2_asym, temp1;
    xalpha = (double*)malloc(sizeof(double)*(nm));
    // residual = (double*)malloc(sizeof(double)*(nm));
    yb = (double*)malloc(sizeof(double)*(nm));
    tn2_asym = (double*)malloc(sizeof(double)*(B));

    tn2 = st3(n,m,p,q,r,x,y,z,u,alpha);
    AbyB(xalpha,x,alpha,nm,p,1);
    // for(i=0; i<nm; i++){
    //     residual[i] = y[i] - xalpha[i]; // m*n
    // }
    for(b=0;b<B;b++){
        if(family == 1){ // gauss
            // v: mn*B
            // temp = b*nm;
            temp = b*n;
            for(i=0; i<nm; i++){
                // yb[i] = xalpha[i] + v[temp+i]*residual[i];
                temp1 = xalpha[i];
                // yb[i] = temp1 + v[temp+i]*(y[i]-temp1);   // v:nm*B
                yb[i] = temp1 + v[temp + i/m]*(y[i]-temp1); // v: n*B
            }
        }
        else if(family==2){ // binomial
            // v: mn*B
            // temp = b*nm;
            temp = b*n;
            for(i=0;i<nm;i++){
                // yb[i] = v[temp+i] < 1.0/(1.0+exp(-xalpha[i]));
                // yb[i] = v[temp+i] < 1.0/(1.0+exp(-xalpha[i]));
                yb[i] = v[temp+i/m] < 1.0/(1.0+exp(-xalpha[i]));
            }
        }
        else{
            // v: mn*B
            // temp = b*nm;
            temp = b*n;
            for(i=0; i<nm; i++){
                // yb[i] = xalpha[i] + v[temp+i]*residual[i];
                temp1 = xalpha[i];
                // yb[i] = temp1 + v[temp+i]*(y[i]-temp1);
                yb[i] = temp1 + v[temp+i/m]*(y[i]-temp1);
            }
        }
        tn2_asym[b] = st3(n,m,p,q,r,x,yb,z,u,alpha);
        if(tn2_asym[b]>tn2){count++;}
    }
    // printArrayDouble(tn2_asym,20,20);
    result[0] = tn2;
    result[1] = (count*1.0)/(B*1.0);
    free(xalpha);
    // free(residual);
    free(yb);
    free(tn2_asym);
}

SEXP _WAST1(SEXP N, SEXP M, SEXP P, SEXP Q, SEXP R, SEXP X, SEXP Y, SEXP Z, SEXP U, SEXP FAMILY, SEXP B, SEXP V){
    double *result,*alpha;
    int p = INTEGER(P)[0];
    result = (double*)malloc(sizeof(double)*2);
    alpha = (double*)malloc(sizeof(double)*p);

    SEXP tn0, pval, list, list_names;
    PROTECT(tn0 		= allocVector(REALSXP, 	1));
	PROTECT(pval 		= allocVector(REALSXP, 	1));
	PROTECT(list 		= allocVector(VECSXP, 	2));
	PROTECT(list_names 	= allocVector(STRSXP, 	2));

    pval3(result,INTEGER(N)[0],INTEGER(M)[0],p,INTEGER(Q)[0],INTEGER(R)[0],REAL(X),REAL(Y),REAL(Z),REAL(U),INTEGER(FAMILY)[0],INTEGER(B)[0],alpha,REAL(V));
    REAL(tn0)[0] = result[0];
    REAL(pval)[0] = result[1];

    SET_STRING_ELT(list_names, 	0,	mkChar("tn0"));
	SET_STRING_ELT(list_names, 	1,	mkChar("pval"));
	SET_VECTOR_ELT(list, 		0, 	tn0);
	SET_VECTOR_ELT(list, 		1, 	pval);
	setAttrib(list, R_NamesSymbol, 	list_names);

    UNPROTECT(4);
	return list;
}

