/* hydrogen.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define ABS(a)	   (((a) < 0) ? -(a) : (a))
#define MIN(a,b)   (((a) < (b)) ? (a) : (b))
#define MAX(a,b)   (((a) > (b)) ? (a) : (b))

main()
{
  /* Subroutines */
    extern void do_mesh(int, double, double, double, double, 
		       double *, double *, double *);
    extern void init_pot (double, int, double *, double *);
    extern double solve_sheq(int, int, double, int, double,
			     double *, double *, double *, double *, 
			     double *);
    /* Variables */
    double *r, *y, *r2, *sqr, *vpot;
    int mesh;
    double zeta, dx, rmax, xmin, zmesh, e;
    int i, l, n;
    FILE *out;

/* initialize atomic charge (Z) */

    fprintf(stdout, " Atomic Charge > " );
    scanf("%lf",&zeta);
    fprintf(stdout, " Atomic Charge = %8.4f\n", zeta );

/* initialize logarithmic mesh */

    zmesh = zeta;
    rmax = 100.;
    xmin = -8.;
    dx = 0.01 ;
/* number of grid points */
    mesh = (int ) ((log(zmesh * rmax) - xmin) / dx);

    r = malloc( (mesh+1) * sizeof(double) );
    r2= malloc( (mesh+1) * sizeof(double) );
    sqr=malloc( (mesh+1) * sizeof(double) );

    do_mesh ( mesh, zmesh, xmin, dx, rmax, r, sqr, r2);

/* initialize the potential */

    vpot = malloc( (mesh+1) * sizeof(double) );
    init_pot(zeta, mesh, r, vpot);

/* open output file that will contain the wavefunctions */

    out = fopen("wfc.out", "w");

/* read number of nodes (stop if nodes < 0) */

L1:
    fprintf(stdout, " n,l >> " );
    scanf("%d %d",&n, &l);
    if (n < 1) {
      fclose(out);
      free(y); free(vpot); free(sqr); free(r2); free(r);
      exit (0);
    }
    if (n < l + 1) {
      fprintf(stderr, "error in main: n.lt.l -> wrong number of nodes\n");
      exit (1);
    }
    if (l < 0) {
      fprintf(stderr, "error in main: l < 0 unphysical\n");
      exit (1);
    }

/* solve the schroedinger equation in radial coordinates by Numerov method */

    y = malloc( (mesh+1) * sizeof(double) );
    e = solve_sheq(n, l, zeta, mesh, dx, r, sqr, r2, vpot, y);

/* write out the eigenvalue energy to be compared with the external potential */

    fprintf (stdout, "eigenvalue = %16.8e, eig*(n/zeta)^2 = %16.8e\n", 
	    e, e*(n*n/zeta/zeta) );

    fprintf(out, "#       r             R(r)          X(r)=rR(r)          e\n"); 
    for (i = 0; i <= mesh; ++i) {
      fprintf(out, "%14.6e %14.6e %14.6e %14.6e\n", 
	      r[i], y[i] / sqr[i], y[i] * sqr[i], e);
    }
    fprintf(out, "\n\n"); 
    goto L1;
} /* MAIN__ */

/* --------------------------------------------------------------------- */
/* Subroutine */ double solve_sheq(int n, int l, double zeta, int mesh, 
				   double dx, double *r, double *sqr, 
				   double *r2, double *vpot, double *y) 

{
    /* Local variables */
    static double eps=1e-10;
    static int maxiter=100;
    int i, j;
    double e, de, fac;
    int icl, kkk;
    double x2l2, elw, eup, ddx12, norm;
    int nodes;
    double sqlhf, ycusp, dfcusp;
    int ncross;
    double *f;

/* --------------------------------------------------------------------- */

/* solve the schroedinger equation in radial coordinates on a 
   logarithmic grid by Numerov method - atomic (Ry) units */

    ddx12 = dx * dx / 12.;
/* Computing 2nd power */
    sqlhf = (l + 0.5) * (l + 0.5);
    x2l2 = (double) (2*l+ 2);
/* set initial lower and upper bounds to the eigenvalue */

    eup = vpot[mesh];
    elw = eup;
    for (i = 0; i <= mesh; ++i) {
/*      if ( elw > sqlhf / r2[i] + vpot[i] ) 
	elw = sqlhf / r2[i] + vpot[i] ; */
        elw = MIN ( elw, sqlhf / r2[i] + vpot[i] );
    }
    if (eup - elw < eps) {
      fprintf (stderr, "%25.16e 25.16e\n", eup, elw);
      fprintf (stderr, "solve_sheq: lower and upper bounds are equal\n");
      exit(1);
    }
    e = (elw + eup) * .5;
    f = malloc( (mesh+1) * sizeof(double) );

/* start loop on energy */

    de= 1e+10; /* any number larger than eps */
    for ( kkk = 0; kkk < maxiter && ABS(de) > eps ; ++kkk ) {

/* set up the f-function and determine the position of its last */
/* change of sign */
/* f < 0 (approximately) means classically allowed   region */
/* f > 0         "         "        "      forbidden   " */

      icl = -1;
      f[0] = ddx12 * (sqlhf + r2[0] * (vpot[0] - e));
      for (i = 1; i <= mesh; ++i) {
	f[i] = ddx12 * (sqlhf + r2[i] * (vpot[i] - e));
/* beware: if f(i) is exactly zero the change of sign is not observed */
/* the following line is a trick to prevent missing a change of sign */
/* in this unlikely but not impossible case: */
	if (f[i] == 0.) {
	    f[i] = 1e-20;
	}
	if (f[i] != copysign(f[i], f[i - 1])) {
	    icl = i;
	}
      }
      if (icl < 0 || icl >= mesh - 2) {
        fprintf (stderr, "%4d %4d\n", icl, mesh);
        fprintf (stderr, "error in solve_sheq: last change of sign too far");
        exit(1);
      }

/* f function as required by numerov method */

      for (i = 0; i <= mesh; ++i) {
	f[i] = 1. - f[i];
	y[i] = 0.;
      }

/* determination of the wave-function in the first two points */

      nodes = n - l - 1;
      y[0] = pow (r[0], l+1) * (1. - zeta * 2. * r[0] / x2l2) / sqr[0];
      y[1] = pow (r[1], l+1) * (1. - zeta * 2. * r[1] / x2l2) / sqr[1];

/* outward integration, count number of crossings */

      ncross = 0;
      for (i = 1; i <= icl-1; ++i) {
 	y[i + 1] = ((12. - f[i] * 10.) * y[i] - f[i - 1] * y[i - 1])
		 / f[i + 1];
	if (y[i] != copysign(y[i],y[i+1]) ) {
	    ++ncross;
	}
      }
      fac = y[icl];

/* check number of crossings */

      if (ncross != nodes) {
        /* incorrect number of nodes: adjust energy bounds */
	if (ncross > nodes) {
	    eup = e;
	} else {
	    elw = e;
	}
	e = (eup + elw) * .5;
      } else {
        /* correct number of nodes: perform inward iteration */

/* determination of the wave-function in the last two points */
/* assuming y(mesh+1) = 0 and y(mesh) = dx */

        y[mesh] = dx;
        y[mesh - 1] = (12. - f[mesh] * 10.) * y[mesh] / f[mesh - 1];

/* inward integration */

        for (i = mesh - 1; i >= icl+1; --i) {
	  y[i - 1] = ((12. - f[i] * 10.) * y[i] - f[i + 1] * y[i + 1])
		 / f[i - 1];
	  if (y[i - 1] > 1e10) {
	    for (j = mesh; j >= i-1; --j) {
		y[j] /= y[i - 1];
	    }
	  }
        }

/* rescale function to match at the classical turning point (icl) */

        fac /= y[icl];
        for (i = icl; i <= mesh; ++i) {
 	  y[i] *= fac;
        }

/* normalize on the segment */

        norm = 0.;
        for (i = 1; i <= mesh; ++i) {
	  norm += y[i] * y[i] * r2[i] * dx;
        }
        norm = sqrt(norm);
        for (i = 0; i <= mesh; ++i) {
	  y[i] /= norm;
        }

/* find the value of the cusp at the matching point (icl) */

        i = icl;
        ycusp = (y[i - 1] * f[i - 1] + f[i + 1] * y[i + 1] + f[i] * 10. 
  	    * y[i]) / 12.;
        dfcusp = f[i] * (y[i] / ycusp - 1.);

/* eigenvalue update using perturbation theory */

        de = dfcusp / ddx12 * ycusp * ycusp * dx;
        if (de > 0.) {
	  elw = e;
        }
        if (de < 0.) {
	  eup = e;
        }

/* prevent e to go out of bounds, i.e. e > eup or e < elw */
/* (might happen far from convergence) */

        e = e + de;
        e = MIN (e,eup);
        e = MAX (e,elw);
        /* if ( e > eup ) e=eup;
           if ( e < elw ) e=elw; */
      }
    }
/* ---- convergence not achieved ----- */
    if ( ABS(de) > eps ) {
      if ( ncross != nodes ) {
         fprintf(stderr, "ncross=%4d nodes=%4d icl=%4d e=%16.8e elw%16.8e eup=%16.8e\n", 
		  ncross, nodes, icl,  e, elw, eup);
      } else {
         fprintf(stderr, "e=%16.8e  de= %16.8e\n", e, de);
      }
      fprintf(stderr, " solve_sheq not converged after %d iterations\n",maxiter);
      exit (1);
    }
/* ---- convergence has been achieved ----- */
    fprintf(stdout, "convergence achieved at iter # %4d, de = %16.8e\n",
	    kkk, de);
    free(f);
    return e;
} /* solve_sheq__ */

/* -------------------------------------------------------------------- */
/* Subroutine */ void do_mesh ( int mesh, double zmesh, double xmin, 
	double dx, double rmax, 
	double *r, double *sqr, double *r2)
{

    /* Builtin functions */
    double log(double), exp(double), sqrt(double);

    /* Local variables */
    int i;
    double x;

/* -------------------------------------------------------------------- */

/* initialize grid */

    for (i = 0; i <= mesh; ++i ) {
	x = xmin + dx * i;
	r[i] = exp(x) / zmesh;
	sqr[i] = sqrt(r[i]);
	r2[i] = r[i] * r[i];
    }
    fprintf(stdout, " radial grid information:\n");
    fprintf(stdout, " dx   = %12.6f", dx);
    fprintf(stdout, ", xmin = %12.6f", xmin);
    fprintf(stdout, ", zmesh =%12.6f\n", zmesh);
    fprintf(stdout, " mesh = %5d", mesh);
    fprintf(stdout, ", r(0) = %12.6f",  r[0]);
    fprintf(stdout, ", r(mesh) = %12.6f\n", r[mesh]);
    return; 
} /* do_mesh */

/* -------------------------------------------------------------------- */
/* Subroutine */ void init_pot(double zeta, int mesh, double *r, 
			      double *vpot)
{
    /* Local variables */
    int i;
    FILE *out;

/* -------------------------------------------------------------------- */

/* initialize potential */

    out = fopen("pot.out","w");
    fprintf(out, "#       r             V(r)\n");
    for (i = 0; i <= mesh; ++i) {
	vpot[i] = -2 * zeta / r[i];
	fprintf(out, "%16.8e %16.8e\n", r[i], vpot[i]);
    }
    fclose(out);
    return;
} /* init_pot */
