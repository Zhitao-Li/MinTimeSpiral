#include <mex.h>

#include <stdio.h>
#include <sys/time.h>
#ifdef userusage
#include <sys/resource.h>
#else
#include <sys/times.h>
#endif

#include <math.h>
#include <float.h>
#include <cblas.h>


#ifndef CLK_TCK
#define CLK_TCK 60
#endif


/* basic macros */
#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))


/* constants for socp algorithm */
#define MAX_ITER_PLANE 20
#define MAX_LAMBDA2 1e-2 /* plane search stopping crit. */
#define DIV_ALPHA 2
#define MIN_ALPHA 1e-6 /* max. of about 20 line search iterations */


#ifdef nounderscores
#define ddot_ ddot
#define dcopy_ dcopy
#define daxpy_ daxpy
#define dscal_ dscal
#define dgemv_ dgemv
#define dsyr_ dsyr
#define dsyrk_ dsyrk
#define dposvx_ dposvx
#define dgelss_ dgelss
#endif


/* BLAS 1 */
double ddot_( );
void dcopy_( );
void daxpy_( );
void dscal_( );

/* BLAS 2 */
void dgemv_( );
void dsyr_( );

/* BLAS 3 */
void dsyrk_( );

/* LAPACK */
void dposvx_( );
void dgelss_( );


/*
*  set vector of n doubles to zero
*/
void dzero(mwSize n, double* p)
{
    /* faster version:
    n*=sizeof(double)/sizeof(char);
    memset(p,0,n);
    */
    /* alternative: (using BLAS to avoid need for memory.h) */
    mwSize int1 = 1;
    double double0 = 0.0;
    dscal_(&n, &double0, p, &int1);
}


/*
* returns sum of elementwise division of x by y
*   in matlab notation: sum(x./y)
*/
double dsumdiv(mwSize n, double *x, double *y)
{
    mwSize i;
    double sumdiv=0.0;
    for (i=0; i<n; ++i)
        sumdiv += x[i]/y[i];
    
    return sumdiv;
}


/*
*  Convert symmetric matrix from upper storage to full storage
*/
void dupge(mwSize n, double *A)
{
    mwSize i, j, k;
    mwSize int1=1;

    for (i=n-1, j=1, k=n; i>0; --i, j+=n+1, k+=n+1)
        dcopy_(&i,A+k,&n,A+j,&int1);
}


/*
*  note: this fct currently not used
*
*  compute duality gap and deviation from centrality
*
*  output:
*   *pgap = duality gap
*   *pdev = deviation from centrality
*/
void dgapdev(mwSize m, mwSize L, mwSize* N,
             double *u, double *z, double *pgap, double *pdev)
     
{
    mwSize i, j, k;
    double fu, fz;

    mwSize int1 = 1;

    /* gap = u'*z; */
    *pgap = ddot_(&m, u, &int1, z, &int1);

    /* dev = -sum([log(SF'*(u.^2));log(SF'*(z.^2))]) + 2*L*(log(gap)-log(L)); */
    *pdev=2*L*(log(*pgap)-log(L));
    for (i=0, k=0; i<L; ++i)
    {
        for (j=0, fu=0.0, fz=0.0; j<N[i]-1; ++j, ++k)
        {
            fu-=SQR(u[k]);
            fz-=SQR(z[k]);
        }
        fu+=SQR(u[k]);
        fz+=SQR(z[k]);
        ++k;
        *pdev-=log(fu)+log(fz);
    }
}


/*
*  compute gradient of potential fct wrt. p and q
*   at u+p*du, z+q*dz
*
*  output:
*   *pgp = gradient wrt. p
*   *pgq = gradient wrt. q
*/
void dgrad(
    double w,
    mwSize L,
    double c1, /* pre-computed constants, dependent on: u, z, du, dz */
    double c2,
    double c3,
    double* d1,
    double* d2,
    double* d3,
    double* e1,
    double* e2,
    double* e3,
    double p, /* du scaling */
    double q, /* dz scaling */
    /* output args. */
    double* pgp, /* return gradient */
    double* pgq,
    double* t1, /* intermediate values that will be re-used in dhess() */
    double* t2,
    double* t3,
    double* t4
)
{
    double ptwo=2*p;
    double psqr=SQR(p);
    double qtwo=2*q;
    double qsqr=SQR(q);
    
    mwSize int1=1;
    
    /* t1 = d2 + p*d3 */
    dcopy_(&L,d2,&int1,t1,&int1);
    daxpy_(&L,&p,d3,&int1,t1,&int1);
    
    /* t2 = d1 + 2*p*d2 + p*p*d3 */
    dcopy_(&L,d1,&int1,t2,&int1);
    daxpy_(&L,&ptwo,d2,&int1,t2,&int1);
    daxpy_(&L,&psqr,d3,&int1,t2,&int1);
    
    /* t3 = e2 + q*e3 */
    dcopy_(&L,e2,&int1,t3,&int1);
    daxpy_(&L,&q,e3,&int1,t3,&int1);
    
    /* t4 = e1 + 2*q*e2 + q*q*e3 */
    dcopy_(&L,e1,&int1,t4,&int1);
    daxpy_(&L,&qtwo,e2,&int1,t4,&int1);
    daxpy_(&L,&qsqr,e3,&int1,t4,&int1);
    
    *pgp=w*c2/(c1+p*c2+q*c3)-2*dsumdiv(L,t1,t2);
    *pgq=w*c3/(c1+p*c2+q*c3)-2*dsumdiv(L,t3,t4);
}


/*
*  compute hessian of primal barrier wrt. p and q
*   at u+p*du, z+q*dz
*
*  output:
*   *php = hessian wrt. p
*   *phq = hessian wrt. q
*/
void dhess(
    /* input args. */
    mwSize L,
    double* d3, /* pre-computed constants */
    double* e3,
    double* t1, /* intermediate values from dgrad() */
    double* t2,
    double* t3,
    double* t4,
    /* output args. */
    double* php,
    double* phq
)    
{
    mwSize i;

    for (i=0, *php=0.0, *phq=0.0; i<L; ++i)
    {
        *php+=(d3[i]*t2[i]-2*SQR(t1[i]))/SQR(t2[i]);
        *phq+=(e3[i]*t4[i]-2*SQR(t3[i]))/SQR(t4[i]);
    }
    *php *= -2;
    *phq *= -2;
}


/*
*  returns workspace size in number of doubles, pointers and ints
*  required by socp()
*  (use a macro in socp.h instead?)
*
*  input arguments:
*   use the same values as for socp()
*
*  output arguments:
*   *mhist-- number of lines in *hist (elements are doubles)
*   *nhist-- number of columns in *hist
*   *ndbl -- number of doubles in *dblwork
*   *nptr -- number of pointers in *ptrwork
*   *nint -- number of integers in *intwork
*/
void socp_getwork(
    /* input args.: */
    mwSize L,
    mwSize* N,
    mwSize n,
    mwSize max_iter,
    mwSize out_mode,
    /* output args.: */
    mwSize* mhist,
    mwSize* nhist,
    mwSize* ndbl,
    mwSize* nint
)
{
    mwSize i;
    mwSize m;
    
    /* compute m */
    for (i=0, m=0; i<L; m+=N[i++]);
    
    *mhist = out_mode;   /* 0: none;  1: gap only;  2: gap and deviation */
    *nhist = max_iter+1; /* for initial point and after each iteration */
    
    *ndbl = 7*m + 2*n + 10*n + 2*SQR(n) + 11*L;
    *nint = n;
}


/*
*  solve second order cone problem
*
*/
int socp(
    mwSize L,
    mwSize* N,
    mwSize n,
    
    double* f,
    double* A,
    double* b,
    
    double* x,
    double* z,
    
    double abs_tol,
    double rel_tol,
    double target,
    mwSize* iter,
    
    double Nu,
    
    mwSize* info,
    mwSize out_mode,
    double* hist,
    
    double* dblwork,
    mwSize* intwork
)
{
    mwSize m;
    double w;
    double gap;
    double dev;
    
    mwSize iter_out;
    double lbound;
    double *u, *fc, *gup, *gu; /* u[m], fc[L], gup[m], gu[m] */
    double *du, *dz; /* du[m], dz[m] */
    double *gx, *dx; /* gx[n], dx[n] */
    double *Hx; /* Hx[n*(n+1)] */
    double s;
    
    /* to record reason for exiting main loop */
    mwSize xabs=0;
    mwSize xrel=0;
    mwSize xtarget=0;
    mwSize xiter=0;
    
    /* for plane search */
    double c1, c2, c3;
    double *d1, *d2, *d3, *e1, *e2, *e3; /* [L] */
    double p, q, gp, gq, hp, hq, dp, dq;
    double *t1, *t2, *t3, *t4; /* [L] */
    double lambda2;
    mwSize iter_plane;
    
    /* for line search */
    double p0, q0;
    double alpha, galpha;
    double *ua, *za; /* ua[m], za[m] */
    double fu, fz;
    mwSize out_barrier;
    
    /* general */
    mwSize i, j, k, l;

    /* for linear system solver */
    /* for both dposvx_() and dgelss_() */
    mwSize use_ls = 0; /* will be changed to 1 if dgelss_ is to be used */
    mwSize la_info; /* exit info from dposvx_ and dgelss_ */
    double *Hwork; /* Hwork[SQR(n)] (dgelss_ only needs n) */
    double rcond;
    /* for dposvx_() */
    mwSize equed;
    double scale, ferr, berr;
    /* for dgelss_() */
    mwSize rank;
    mwSize lwork = 10*n; /* size of workspace, minimum is 5*n */
    
    /* constant variables for BLAS */
    mwSize int1=1;
    double double1 = 1.0, double0 = 0.0, double_1 = -1.0;



    /* compute m: size of dual variable */
    for (i=0, m=0; i<L; m+=N[i++]);
    
    /*
    * organize workspace
    *
    * total space needed =
    *   doubles:   7*m + 2*n + max(3*n,lwork) + 2*SQR(n) + 11*L
    *   ints:      n
    */

    /* dblwork is used for dposvx_() and dgelss_(), */
    /*   need 3*n and lwork doubles respectively. */
    u = dblwork + lwork;    /* lwork=10*n */
    fc = u + m;
    gup = fc + L;
    gu = gup + m;
    du = gu + m;
    dz = du + m;
    gx = dz + m,
    dx = gx + n;
    Hx = dx + n;
    d1 = Hx + SQR(n);
    d2 = d1 + L;
    d3 = d2 + L;
    e1 = d3 + L;
    e2 = e1 + L;
    e3 = e2 + L;
    t1 = e3 + L;
    t2 = t1 + L;
    t3 = t2 + L;
    t4 = t3 + L;
    ua = t4 + L;
    za = ua + m;
    Hwork = za + m;
    /* Hwork needs SQR(n) doubles */


    /* gap reduction vs. centering */
    w=2*L + Nu*sqrt(2*L);
    
    /* u=A*x+b */
    dcopy_(&m, b, &int1, u, &int1); /* u=b */
    dgemv_("N", &m, &n, &double1, A, &m, x, &int1, &double1, u, &int1); /* u=A*x+u */

    /* compute gap (and deviation from centrality and store in hist) */
    if (out_mode == 2)
    {
        dgapdev(m, L, N, u, z, &gap, &dev);
        hist[0]=gap;
        hist[1]=dev;
    }
    else
    {
        gap=ddot_(&m, u, &int1, z, &int1); /* gap = u'*z; */
        if (out_mode == 1)
            hist[0]=gap;
    }

    /* outer loop */
    iter_out=0;
    while (!((rel_tol < 0.0 &&
              (xtarget = (ddot_(&n, f, &int1, x, &int1) < target ||
              -ddot_(&m, b, &int1, z, &int1) >= target))) ||
             (abs_tol > 0.0 && (xabs = (gap <= abs_tol))) ||
             (rel_tol > 0.0 &&
             (xrel = (((lbound =- ddot_(&m, b, &int1, z, &int1)) > 0.0 &&
              gap/lbound <= rel_tol) ||
              ((lbound =- ddot_(&n, f, &int1, x, &int1)) > 0.0 &&
              gap/lbound<=rel_tol)))) ||
             (xiter = (iter_out >= *iter))))
    {
        ++iter_out;
        
        /* compute gup (gradient of primal barrier wrt. u) */
        /* also, compute fc(i)=2/(t(i)^2-u(i)'*u(i)) */
        for (i=0, k=0; i<L; ++i)
        {
            for (j=0, fc[i]=0.0; j<N[i]-1; ++j)
                fc[i]-=SQR(u[k+j]);
            fc[i]+=SQR(u[k+j]);
            fc[i]=2.0/fc[i];
            for (j=0; j<N[i]-1; ++j, ++k)
                gup[k]=fc[i]*u[k];
            gup[k]=-fc[i]*u[k]; ++k;
        }

        /* compute gu (gradient of potential wrt. u) */
        s = w/gap;
        dcopy_(&m, gup, &int1, gu, &int1); /* gu=gup */
        daxpy_(&m, &s, z, &int1, gu, &int1); /* gu=gu+(w/gap)*z */

        /* compute Hx = A'*Hu*A */
        /* where   Hu(i) = fc(i)*diag([1 1 ... 1 -1]) + gup(i)*gup(i)' */
        /* Hx=0 */
        dzero(SQR(n), Hx);

        for (i=0, k=0; i<L; k+=N[i++])
        { /* for each constraint */
            /* n. of rows of A(i) */
            j = N[i]-1;
            
            /* Hx = Hx + fc(i)*A(i)'*A(i) */
            if (j>0)
                dsyrk_("U", "T", &n, &j, &(fc[i]), &(A[k]), &m, &double1, Hx, &n);
            
            /* Hx = Hx - fc(i)*c(i)*c(i)'  (rank one update) */
            s = -fc[i];
            dsyr_("U", &n, &s, &(A[k+j]), &m, Hx, &n);
            
            /* gx = [A(i);c(i)']'*gup(i)  (this is not gx, just used as aux.) */
            dgemv_("T", &(N[i]), &n, &double1, &(A[k]), &m, &(gup[k]), &int1, &double0, gx, &int1);
            
            /* Hx = Hx + gx*gx'  (rank one update) */
            dsyr_("U", &n, &double1, gx, &int1, Hx, &n);
        }

        /* solve linear system: dx = -Hx\(A'*gu) */
        /* gx = -A'*gu */
        dgemv_("T", &m, &n, &double_1, A, &m, gu, &int1, &double0, gx, &int1);
        /* dx = Hx\gx */
        if (!use_ls)
        { 
            /* solve linear system by QR fact. */
            dposvx_("N", "U", &n, &int1, Hx, &n, Hwork, &n, &equed, &scale,
                    gx, &n, dx, &n,
                    &rcond, &ferr, &berr, dblwork, intwork, &la_info);
            if (la_info>0) /* from 1 to n, Hessian not positive def.; */
                use_ls = 1; /* n+1, Hessian badly conditioned; */
                            /* do SVD now and switch to SVD for all */
                            /* future iterations */
        }
        if (use_ls)
        {
            /* solve linear system in least squares sense using SVD */
            dupge(n,Hx); /* convert to general storage */
            rcond =- 1; /* keep singular values down to machine precision */
            dgelss_(&n, &n, &int1, Hx, &n, gx, &n,
                    Hwork, &rcond, &rank, /* (only first n of Hwork used, for S) */
                    dblwork, &lwork, &la_info); /* dblwork: lwork doubles (>= 5*n) */
            dcopy_(&n, gx, &int1, dx, &int1); /* dx=gx (dgelss_ overwrites gx) */
        }
        if (la_info)
            return la_info; /* abort: return error in lapack solver */

        /* du = A*dx */
        dgemv_("N", &m, &n, &double1, A, &m, dx, &int1, &double0, du, &int1);
        
        /* dz = -(gu+Hu*du) */
        /* computed one constraint at a time: */
        /* dz(i)= -(gu(i)+fc(i)*[du(i);-dt(i)]+gup(i)*(gup(i)'*[du(i);dt(i)])) */
        /* dz = gu */
        dcopy_(&m, gu, &int1, dz, &int1);
        for (i=0, k=0; i<L; k+=N[i++])
        {
            j = N[i]-1;
            /* dz(i) = dz(i) + fc[i]*[du(i);-dt(i)] */
            daxpy_(&j, &(fc[i]), &(du[k]), &int1, &(dz[k]), &int1);
            dz[k+j] -= fc[i]*du[k+j];
            /* s = gup(i)'*du(i) */
            s = ddot_(&(N[i]), &(gup[k]), &int1, &(du[k]), &int1);
            /* dz(i) = dz(i) + s*gup(i) */
            daxpy_(&(N[i]), &s, &(gup[k]), &int1, &(dz[k]), &int1);
        }
        /* dz = - dz */
        dscal_(&m, &double_1, dz, &int1);
        /* optional: scale dz by 1/rho=gap/w
          s=gap/w;
          dscal_(&m,&s,dz,&int1);
        */

        /*
         *  constants for plane search
         */
        c1=gap;
        c2=ddot_(&m,du,&int1,z,&int1);
        c3=ddot_(&m,u,&int1,dz,&int1);
        for (i=0, k=0; i<L; ++i)
        {
            d1[i]=0.0;
            d2[i]=0.0;
            d3[i]=0.0;
            e1[i]=0.0;
            e2[i]=0.0;
            e3[i]=0.0;
            for (j=0; j<N[i]-1; ++j, ++k)
            {
                d1[i]-=SQR(u[k]);
                d2[i]-=u[k]*du[k];
                d3[i]-=SQR(du[k]);
                e1[i]-=SQR(z[k]);
                e2[i]-=z[k]*dz[k];
                e3[i]-=SQR(dz[k]);
            }
            d1[i]+=SQR(u[k]);
            d2[i]+=u[k]*du[k];
            d3[i]+=SQR(du[k]);
            e1[i]+=SQR(z[k]);
            e2[i]+=z[k]*dz[k];
            e3[i]+=SQR(dz[k]); ++k;
        }

        /* plane search loop */
        p=0.0;
        q=0.0;
        for (lambda2=2*MAX_LAMBDA2, iter_plane=0;
             lambda2>MAX_LAMBDA2 && iter_plane<MAX_ITER_PLANE;
             ++iter_plane)
        {
            /* compute gradient and Hessian wrt. p and q */
            /*   at u+p*du, z+q*dz */
            dgrad(w, L, c1, c2, c3, d1, d2, d3, e1, e2, e3, p, q, &gp, &gq, t1, t2, t3, t4);
            dhess(L, d3, e3, t1, t2, t3, t4, &hp, &hq);
            
            /* Newton step */
            dp =- gp/hp;
            dq =- gq/hq;

            /* line search loop: scale down step */
            /*   until inside feasible region */
            /*   (and until between previous point and line minimum) */
            alpha=1;  /* scaling factor for dp, dq */
            p0=p;
            q0=q;
            while (1)
            {
                p=p0+alpha*dp;
                q=q0+alpha*dq;
                
                /* ua=u+p*du */
                dcopy_(&m,u,&int1,ua,&int1);
                daxpy_(&m,&p,du,&int1,ua,&int1);
                /* za=z+q*dz */
                dcopy_(&m,z,&int1,za,&int1);
                daxpy_(&m,&q,dz,&int1,za,&int1);

                /* check constraints: */
                /*   fu = t(i)^2 - u(i)'*u(i) > 0  and  t(i) > 0 */
                /*   fz = s(i)^2 - z(i)'*z(i) > 0  and  s(i) > 0 */
	            for (i=0, k=0, out_barrier=0; i<L && !out_barrier; ++i)
                {
	                for (j=0, fu=0.0, fz=0.0; j<N[i]-1; ++j, ++k)
                    {
	                    fu-=SQR(ua[k]);
	                    fz-=SQR(za[k]);
	                }
	                fu+=SQR(ua[k]);
	                fz+=SQR(za[k]);
	                if (fu<=0 || fz<=0 || ua[k]<=0 || za[k]<=0)
                    {
	                    out_barrier = 1; /* (replace <=0 with <TINY ?) */
	                    break;
                    }
	                ++k;
	            }

	            if (!out_barrier)
                {
	                /* compute gradient along search line wrt. alpha */
	                dgrad(w, L, c1, c2, c3, d1, d2, d3, e1, e2, e3, p, q, &gp, &gq, t1, t2, t3, t4);
	                galpha = dp*gp + dq*gq;
	                /* exit if: all barriers ok (i.e. out_barrier==0) */
	                /*  and gradient negative (=> between initial and minimum) */
	                if (galpha<=0)
	                    break; /* EXIT line search */
	            }
                alpha/=DIV_ALPHA; /* scale down the p, q step */
                if (alpha<MIN_ALPHA)
                {
                    /* line search failed */
                    alpha=0.0;
                    p=p0;
                    q=q0;
                    break; /* EXIT line search */
                }
            } /* end of line search loop */

            /* plane search convergence criterium */
            lambda2=hp*SQR(dp)+hq*SQR(dq);
        }    /* end of plane search loop */
        /* x=x+p*dx */
        daxpy_(&n, &p, dx, &int1, x, &int1);
        /* z=z+q*dz */
        daxpy_(&m, &q, dz, &int1, z, &int1);
        /* u=A*x+b*/
        dcopy_(&m, b, &int1, u, &int1);
        dgemv_("N", &m, &n, &double1, A, &m, x, &int1, &double1, u, &int1);

        /* update gap (and deviation from centrality and store in hist) */
        if (out_mode == 2)
        {
            dgapdev(m,L,N,u,z,&gap,&dev);
            hist[2*iter_out]=gap;
            hist[2*iter_out+1]=dev;
        }
        else
        {
            gap=ddot_(&m,u,&int1,z,&int1);  /* gap = u'*z; */
            if (out_mode == 1)
                hist[iter_out]=gap;
        }
    } /* end of outer loop */
    *iter=iter_out;

    /* report reason for exit */
    *info = xabs + 2*xrel + 3*xtarget + 4*xiter;
    /* normal exit */
    return 0;
} /* end of socp() */


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    mwSize m, n, L, *N, iter, out_mode;
    double* f;
    double* A;
    double* b;
    double* x;
    double* z;
    double* x_opt;
    double* z_opt;
    double abs_tol, rel_tol, target, Nu;
    double* hist;
    mwSize mhist, nhist;
    mwSize ndbl, nptr, nint;
    double* dblwork;
    mwSize* intwork;
    mwSize i, j;
    mwSize info_time, info_socp, info_exit;
    double* info_ptr;
    static mwSize firstcall = 1;
    double* time;
#ifdef userusage
    struct rusage stats;
#else
    struct tms stats;
#endif
    mwSize int1 = 1;


    if (firstcall)
    {
        fprintf(stdout, "\nThis is the beta version of SOCP.,\n");
        fprintf(stdout, "COPYRIGHT (c) 1997, Miguel Lobo, Lieven Vandenberge, Stephen Boyd.\n\n");
        firstcall = 0;
    }

    /* check number of arguments */
    if (nrhs == 12)
    {
        if (nlhs != 5)
        mexErrMsgTxt("Five output arguments required.\n");
    }
    else
        mexErrMsgTxt("Twelve input arguments required.\n");

    /* get dimensions and ptr. to A */
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    A = mxGetPr(prhs[1]);

    /* check that dimensions of f, b, x and z match with A and get ptrs. */
    if (mxGetN(prhs[0]) != 1)
        mexErrMsgTxt("1st input argument must be a column vector.");
    if (mxGetM(prhs[0]) != n)
        mexErrMsgTxt("f and A do not agree in size.");
    f = mxGetPr(prhs[0]);

    if (mxGetN(prhs[2]) != 1)
        mexErrMsgTxt("3rd input argument must be a column vector.");
    if (mxGetM(prhs[2]) != m)
        mexErrMsgTxt("A and b do not agree in size.");
    b = mxGetPr(prhs[2]);

    if (mxGetN(prhs[4]) != 1)
        mexErrMsgTxt("5th input argument must be a column vector.");
    if (mxGetM(prhs[4]) != n)
        mexErrMsgTxt("A and x do not agree in size.");
    x = mxGetPr(prhs[4]);

    if (mxGetN(prhs[5]) != 1)
        mexErrMsgTxt("6th input argument must be a column vector.");
    if (mxGetM(prhs[5]) != m)
        mexErrMsgTxt("A and z do not agree in size.");
    z = mxGetPr(prhs[5]);

    /* get number of cones (=L) and size of each (=N[i]) */
    if (MIN(mxGetM(prhs[3]), mxGetN(prhs[3])) != 1)
        mexErrMsgTxt("4th input argument must be a vector.\n");
    L = MAX(mxGetM(prhs[3]), mxGetN(prhs[3]));
    N = mxCalloc(L, sizeof(mwSize)); /* new N of type int */
    for (i=0, j=0; i<L; ++i)
    {
        N[i] = (mwSize) mxGetPr(prhs[3])[i];
        j += N[i];
    }
    if (j != m)
        mexErrMsgTxt("Sum of elements of N must equal number of lines in A.\n");

    /* get stopping criteria */
    abs_tol = mxGetScalar(prhs[6]);
    rel_tol = mxGetScalar(prhs[7]);
    target = mxGetScalar(prhs[8]);
    iter = (mwSize) mxGetScalar(prhs[9]);

    /* get gap reduction vs. centering factor */
    Nu = mxGetScalar(prhs[10]);
    if (Nu < 0)
        mexErrMsgTxt("Parameter Nu must be non-negative");
    if (Nu <= 1)
        fprintf(stdout, "SOCP warning: Nu > 1 recommended.\n");

    /* output mode */
    out_mode = (mwSignedIndex) mxGetScalar(prhs[11]);
    if (out_mode<0 || out_mode>2)
        mexErrMsgTxt("Illegal value for out_mode; must be 0, 1 or 2.");

    /* allocate workspace for socp() */
    socp_getwork(L,N,n,iter,out_mode,&mhist,&nhist,&ndbl,&nint);

    dblwork = (double*) mxCalloc(ndbl, sizeof(double));
    intwork = (mwSize*) mxCalloc(nint, sizeof(mwSize));

    /* prepare output arguments */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    x_opt = mxGetPr(plhs[0]);
    dcopy_(&n, x, &int1, x_opt, &int1); /* copy x to x_opt */

    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    info_ptr = mxGetPr(plhs[1]);

    plhs[2] = mxCreateDoubleMatrix(m, 1, mxREAL);
    z_opt = mxGetPr(plhs[2]);
    dcopy_(&m, z, &int1, z_opt, &int1); /* copy z to z_opt */

    plhs[3] = mxCreateDoubleMatrix(mhist, nhist, mxREAL);
    hist = mxGetPr(plhs[3]);

    plhs[4] = mxCreateDoubleMatrix(1, 3, mxREAL); /* time stats: utime, stime, iters */
    time = mxGetPr(plhs[4]);


    /*
    * call socp
    */

    /* initialize user time and system time */
#ifdef userusage
    info_time = getrusage(RUSAGE_SELF,&stats);
    time[0] = - (double)stats.ru_utime.tv_sec
              - ((double)stats.ru_utime.tv_usec)/1e6;
    time[1] = - (double)stats.ru_stime.tv_sec
              - ((double)stats.ru_stime.tv_usec)/1e6;
#else
    info_time = times(&stats);
    time[0] = - (double)stats.tms_utime/CLK_TCK;
    time[1] = - (double)stats.tms_stime/CLK_TCK;
#endif

    /* call socp() */
    info_socp = socp(L, N, n, f, A, b,
	                 x_opt, z_opt,
	                 abs_tol,rel_tol,target,&iter,
		             Nu, &info_exit, out_mode, hist,
	                 dblwork, intwork);

    /* update user time and system time */
#ifdef userusage
    info_time = getrusage(RUSAGE_SELF,&stats);
    time[0] += (double)stats.ru_utime.tv_sec
            + ((double)stats.ru_utime.tv_usec)/1e6;
    time[1] += (double)stats.ru_stime.tv_sec
            + ((double)stats.ru_stime.tv_usec)/1e6;
#else
    info_time = times(&stats);
    time[0] += (double)stats.tms_utime/CLK_TCK;
    time[1] += (double)stats.tms_stime/CLK_TCK;
#endif

    time[2] = iter;

    *info_ptr = (double)info_exit;

    /* rescale hist */

    /* note: hist will padded with zeros after call to socp; the following */
    /* line truncates it */
    /* (if hist is used to output other info, this may need to be removed) */
    mxSetN(plhs[3],iter+1);

    /* free allocated memory */
    mxFree(N);
    mxFree(dblwork);
    mxFree(intwork);

    /* error handling */

    if (info_socp)
    {
        fprintf(stdout, "info from dgelss = %ld \n", info_socp);
        mexErrMsgTxt("Error in SOCP, call to LAPACK failed.\n");
    }
}
