#include <common.cuh>
#include <config.cuh>

#ifdef PREVIOUS_CODE



//#include <general.h>

int i_dnnt(double x)
{
    return (int)(x >= 0. ? floor(x + 0.5) : -floor(0.5 - x));
}

int maxX(int a, int b) {
    return (a > b) ? a : b;
}

#define IPIX 1000
#define IPIX1 500

double absX(double a) {
    return (a > 0) ? a : -a;
}

double diff(double a, double b, double old) {
    double tmp = /* 2 * */ absX(a - b)/(a + b);
    return (tmp > old) ? tmp : old;
}


double H(double a) {
    return (a >= 0) ? 1.0f : 0.0f;
}

char fileread[64], filewrite[64], curve[64];
int sigma;

double gauss_(int, double , double, double *, int);
double power_(int, double , double, double *, int);
double ld_   (int, double , double, double *, int);
double ad_   (int, double , double, double *, int);

double cosalpha, sinalpha;
double pix_double__[1000*1000]	/* was [1000][1000] */, x[505], y[505];

int i0, i1, i2;
int ix, iy, ix1, ix2, iy1, iy2;
int pixmax;

int16_t pix[1000*1000]	/* was [1000][1000] */;
int16_t pix1[500*500]	/* was [500][500] */;

double xxx, yyy;
int ixxx, iyyy;
double alpha, x_end__, y_end__;
double value, slope;
double x_diff__, y_diff__;
double x_start__, y_start__;

int main(int argc, char** argv)
{
    double x0, x1, x2, y1, y2, y0;
    double pixL, R_1_2;

    int32_t p_size = (IPIX * IPIX + IPIX1*IPIX1) * sizeof(int16_t);

    if (argc == 1) {
        printf("Usage:\n\t%s iter pixlength sourcesize\n\t\titer - iteration number, e.g. 1\n\t\tpixlength - length of pixel in Einstein radii, e.g. 0.04\n\t\tsourcesize - size of the source (R_1/2) in pixlength, e.g. 5\n", argv[0]);
        exit(-1);
    }

    sprintf(fileread, "IRIS%s", argv[1]);
    sprintf(filewrite, "IRIS%s-track", argv[1]);

/*      reads   unformatted  (pixel)  data of 1000^2 pixels from file IRIS111 */
/*       draws a line through two points, determines lightcurve */
/* 	for circular source, radius sigma,  with flat profile */
/* 	   (Gaussian source with width sigma possible as well) */
/* 					Joachim Wambsganss, January 2008 */
/* 							jkw, 2/27/96 */

    ix1 = 0;
    iy1 = 500;
    ix2 = 1000;
    iy2 = 500;

/* source size in pixels (Gaussian width) */

    sigma = 3;

    x1 = (double) ix1;
    x2 = (double) ix2;
    y1 = (double) iy1;
    y2 = (double) iy2;


    sscanf(argv[2], "%lf", &pixL);
    sscanf(argv[3], "%d", &sigma);

    printf("Pixlength (pixL): %lf\n", pixL);
    printf("Sourcesize (R_1/2): %d\n", sigma);

    printf("fileread = %s\n", fileread);
    FILE * in_file = fopen(fileread, "r");
    if (!in_file)
    {
        printf("Unable to open file!\n");
        return 1;
    }
    int32_t f_size = 0;
    int res = 0;
    res = fread(&f_size, sizeof(int32_t), 1, in_file);
    if (f_size != p_size) {
        printf("Error: The header (first int32_t value) should be %d but is %d\n", p_size, f_size);
        return (-1);
    }
    res = fread(&pix, sizeof(int16_t), IPIX*IPIX, in_file);
    res = fread(&pix1, sizeof(int16_t), IPIX1*IPIX1, in_file);
    f_size = 0;
    res = fread(&f_size, sizeof(int32_t), 1, in_file);
    if (f_size != p_size) {
        printf("Error: The last int32_t value should be %d but is %d\n", p_size, f_size);
        return (-1);
    }
    fclose(in_file);

/* determine maximum of magnification pattern: */
    pixmax = 0;
    for (i1 = 1; i1 <= 1000; ++i1) {
        for (i2 = 1; i2 <= 1000; ++i2) {
            pixmax = maxX(pixmax, pix[i2 + i1 * 1000 - 1001]);
            pix_double__[i2 + i1 * 1000 - 1001] = pow(10.0, (double) ((double) (pix[i2 + i1 * 1000 - 1001] - 1024) / 256.f * .4f));
        }
    }
    printf("pixmax = %d\n", pixmax);
/* determine points of line: */

    slope = (y2 - y1) / (x2 - x1);
    printf("slope = %lf\n", slope);

    alpha = atan(slope);
    printf("alpha = %lf\n", alpha);

    sinalpha = sin(alpha);
    cosalpha = cos(alpha);

/* points for x = 1    and    y = 1: */

    y0 = y1 + (1. - x1) * slope;
    x0 = x1 + (1. - y1) / slope;

    printf("x0 = %lf\n", x0);
    printf("y0 = %lf\n", y0);

    if (x0 >= 1. && x0 <= 1e3) {
        x_start__ = x0;
        y_start__ = 1.0;
    } else if (y0 >= 1. && y0 <= 1e3) {
        x_start__ = 1.0;
        y_start__ = y0;
    }

    printf("x_start = %lf\n", x_start__);
    printf("y_start = %lf\n", y_start__);

/* points for x = ipix    and    y = ipix: */

    y0 = y1 + (1e3 - x1) * slope;
    x0 = x1 + (1e3 - y1) / slope;

    printf("x0 = %lf\n", x0);
    printf("y0 = %lf\n", y0);

    if (x0 >= 1. && x0 <= 1e3) {
        x_end__ = x0;
        y_end__ = y1 + (x_end__ - x1) * slope;
    } else if (y0 >= 1. && y0 <= 1e3) {
        y_end__ = y0;
        x_end__ = x1 + (y_end__ - y1) / slope;
    }

    printf("x_end = %lf\n", x_end__);
    printf("y_end = %lf\n", y_end__);

    x_diff__ = x_end__ - x_start__;
    y_diff__ = y_end__ - y_start__;
    printf("x_diff = %lf\n", x_diff__);
    printf("y_diff = %lf\n", y_diff__);


/* mark pixels along line: */

/* AND   determine light curve for Gaussian source: */

    sprintf(curve, "outline%s_gauss.dat", argv[1]);
    FILE *gauss_out = fopen(curve, "w");
    if (!gauss_out) {
        printf("Can't open out_line file for write!\n");
        exit(-2);
    }

    sprintf(curve, "outline%s_power.dat", argv[1]);
    FILE *power_out = fopen(curve, "w");
    if (!power_out) {
        printf("Can't open out_line file for write!\n");
        exit(-2);
    }

    sprintf(curve, "outline%s_ld.dat", argv[1]);
    FILE *ld_out = fopen(curve, "w");
    if (!ld_out) {
        printf("Can't open out_line file for write!\n");
        exit(-2);
    }

    sprintf(curve, "outline%s_ad.dat", argv[1]);
    FILE *ad_out = fopen(curve, "w");
    if (!ad_out) {
        printf("Can't open out_line file for write!\n");
        exit(-2);
    }

    sprintf(curve, "diff.dat", argv[1]);
    FILE *diff_out = fopen(curve, "a");
    if (!diff_out) {
        printf("Can't open out_line file for write!\n");
        exit(-2);
    }


    double value_gauss, value_power, value_ld, value_ad;

    double diff_gauss_power=0, diff_gauss_ld=0, diff_gauss_ad=0,
                               diff_power_ld=0, diff_power_ad=0,
                                                   diff_ld_ad=0;

/*        do i0 = 0,2*ipix */
    for (i0 = -2000; i0 <= 2000; ++i0) {
        xxx = x_start__ + i0 * cosalpha;
        yyy = y_start__ + i0 * sinalpha;
        ixxx = i_dnnt(xxx);
        iyyy = i_dnnt(yyy);
        if (ixxx >= sigma * 3 + 1 && ixxx <= 1000 - sigma * 3 && iyyy >= sigma * 3 + 1 && iyyy <= 1000 - sigma * 3) {
/* cc	         value = 10**(float(pix(ixxx,iyyy)-1024)/256.0) */

            value_gauss = gauss_(sigma, xxx, yyy, pix_double__, IPIX);
            fprintf(gauss_out, "%5d %15.4lf %15.4lf %10d %12.5lf %5d %5d\n", i0, xxx, yyy, pix[ixxx + iyyy * 1000 - 1001], value_gauss, ixxx, iyyy);

            value_power = power_(sigma, xxx, yyy, pix_double__, IPIX);
            fprintf(power_out, "%5d %15.4lf %15.4lf %10d %12.5lf %5d %5d\n", i0, xxx, yyy, pix[ixxx + iyyy * 1000 - 1001], value_power, ixxx, iyyy);

            value_ld = ld_(sigma, xxx, yyy, pix_double__, IPIX);
            fprintf(ld_out, "%5d %15.4lf %15.4lf %10d %12.5lf %5d %5d\n", i0, xxx, yyy, pix[ixxx + iyyy * 1000 - 1001], value_ld, ixxx, iyyy);

            value_ad = ad_(sigma, xxx, yyy, pix_double__, IPIX);
            fprintf(ad_out, "%5d %15.4lf %15.4lf %10d %12.5lf %5d %5d\n", i0, xxx, yyy, pix[ixxx + iyyy * 1000 - 1001], value_ad, ixxx, iyyy);

            diff_gauss_power = diff(value_gauss, value_power, diff_gauss_power);
            diff_gauss_ld    = diff(value_gauss, value_ld,    diff_gauss_ld   );
            diff_gauss_ad    = diff(value_gauss, value_ad,    diff_gauss_ad   );

            diff_power_ld    = diff(value_power, value_ld,    diff_power_ld   );
            diff_power_ad    = diff(value_power, value_ad,    diff_power_ad   );

            diff_ld_ad       = diff(value_ld,    value_ad,    diff_ld_ad      );

        }
    }

    fprintf(diff_out, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", diff_gauss_power, diff_gauss_ld, diff_gauss_ad, diff_power_ld, diff_power_ad, diff_ld_ad);


    for (i0 = -2000; i0 <= 2000; ++i0) {
        xxx = x_start__ + i0 * cosalpha;
        yyy = y_start__ + i0 * sinalpha;
        ixxx = i_dnnt(xxx);
        iyyy = i_dnnt(yyy);
        if (ixxx >= sigma * 3 + 1 && ixxx <= 1000 - sigma * 3 && iyyy >= sigma * 3 + 1 && iyyy <= 1000 - sigma * 3) {
            pix[ixxx + iyyy * 1000 - 1001] = 0;
            pix[ixxx + iyyy * 1000 - 1001] = 2000;
            pix[ixxx + iyyy * 1000 - 1001] = 5000;
        }
    }

/* 030711: indicate source in lower left corner: */
    for (i1 = sigma * -2; i1 <= sigma << 1; ++i1) {
	for (i2 = sigma * -2; i2 <= sigma << 1; ++i2) {
	    if (i1 * i1 + i2 * i2 <= (sigma << 2) * sigma) {
		pix[(sigma << 1) + 10 + i1 + ((sigma << 1) + 10 + i2) * 1000 
			- 1001] = (int16_t) ((pixmax << 1) / 3);
	    }
	}
    }
    for (i1 = -sigma; i1 <= sigma; ++i1) {
	for (i2 = -sigma; i2 <= sigma; ++i2) {
	    if (i1 * i1 + i2 * i2 <= sigma * sigma) {
		pix[(sigma << 1) + 10 + i1 + ((sigma << 1) + 10 + i2) * 1000 
			- 1001] = (int16_t) pixmax;
	    }
	}
    }


/* convert x and y into pixels: */
/* 		and add pixel values to magpat: */

    for (i1 = 1; i1 <= 505; ++i1) {
	ix = (x[i1 - 1] + 1.1f) * 1000 / 2.2f + 1;
	iy = (y[i1 - 1] + 1.1f) * 1000 / 2.2f + 1;
/* 	write(*,*)' i1,x(i1),ix,y(i1),iy: ', i1,x(i1),ix,y(i1),iy */
	if (ix >= 1 && ix <= 1000 && iy >= 1 && iy <= 1000) {
/* 	   pix(ix,iy) = 5000 */
/* 	   pix(ix,1000-iy) = 5000 */
	    pix[ix + iy * 1000 - 1001] = (int16_t) (pix[ix + iy * 1000 - 1001] + 500);
	    pix[ix + (1000 - iy) * 1000 - 1001] = (int16_t) (pix[ix + (1000 - iy) * 1000 - 1001] + 500);
	} else {
/* 	   write(*,*)' i1,x(i1),ix,y(i1),iy: ', i1,x(i1),ix,y(i1),iy */
	}
    }

    printf("filewrite = %s\n", filewrite);
    FILE * out_file = fopen(filewrite, "w");
    if (!in_file)
    {
        printf("Unable to open file for write!\n");
        return 1;
    }
    fwrite(&p_size, sizeof(int32_t), 1, out_file);
    fwrite(&pix, sizeof(int16_t), IPIX*IPIX, out_file);
    fwrite(&pix1, sizeof(int16_t), IPIX1*IPIX1, out_file);
    fwrite(&p_size, sizeof(int32_t), 1, out_file);
    fclose(out_file);

    return 0;
}

double gauss_(int isig, double xxx, double yyy, double* pix, int ipix)
{
    /* System generated locals */
    double value = 0;
    int pix_dim1, pix_offset;

    /* Local variables */
    double factorex;
    int i1, i2, ix, iy;
    double dum, delx, dely;
    int isig3;
    double dist2, sigsq2;
    int isig3sq;
    double normfac;

/* modified: disk of constant brightness */

    /* Parameter adjustments */
    pix_dim1 = ipix;
    pix_offset = 1 + pix_dim1;
    pix -= pix_offset;

    /* Function Body */
    ix = i_dnnt(xxx);
    iy = i_dnnt(yyy);
    delx = xxx - (double) ix;
    dely = yyy - (double) iy;


/* 030711: isig3 = RADIUS  for source!!! */

    isig3 = isig * 4;
/* 	isig3 =   isig */
    isig3sq = isig3 * isig3;
/* Computing 2nd power */
    sigsq2 = ((double)isig) * ((double)isig);
    value = 0.f;
/* 	 write(*,*)' in GAUSS' */
/* write(*,*)' isig,ix,iy,value = ',isig,ix,iy,value */
/* 	pause */

/* Gaussian source profiles with width  isig: */
    normfac = 0.f;
    int count = 0;
    for (i1 = -isig3; i1 <= isig3; ++i1) {
	for (i2 = -isig3; i2 <= isig3; ++i2) {

/* 	            dist2 = i2*i2 + i1*i1 */

	    dist2 = pow(i1 - delx ,2) + pow(i2 - dely ,2);
	    if (isig != 0) {
		if (dist2 <= (double) isig3sq) {
		    factorex = exp(-dist2 / sigsq2);
/* 	                   factorex = 1.0 */
		    dum = pix[ix + i1 + (iy + i2) * pix_dim1];
		    value += dum * factorex;
		    normfac += factorex;
		    count++;
		}
	    } else if (isig == 0) {
		factorex = 1.f;
		dum = pix[ix + i1 + (iy + i2) * pix_dim1];
		value += dum * factorex;
		normfac += factorex;
	    }
	}
    }
//    printf("normfac = %f\n", normfac);
    value /= normfac;
    return value;
}

double ld_(int isig, double xxx, double yyy, double* pix, int ipix)
{
    double p_ld = static_cast<double>(2);
    double R_1_2 = (double)isig * sqrt(log(2));
    double R = R_1_2/sqrt(1-pow(static_cast<double>(0.5),static_cast<double>(2)/(p_ld+2)));
    double R2 = R * R;

    /* System generated locals */
    double value = 0;
    int pix_dim1, pix_offset;

    /* Local variables */
    double factorex;
    int i1, i2, ix, iy;
    double dum, delx, dely;
    int isig3;
    double dist2, sigsq2, dist;
    int isig3sq;
    double normfac;

/* modified: disk of constant brightness */

    /* Parameter adjustments */
    pix_dim1 = ipix;
    pix_offset = 1 + pix_dim1;
    pix -= pix_offset;

    /* Function Body */
    ix = i_dnnt(xxx);
    iy = i_dnnt(yyy);
    delx = xxx - (double) ix;
    dely = yyy - (double) iy;


/* 030711: isig3 = RADIUS  for source!!! */

    isig3 = isig * 10;
/* 	isig3 =   isig */
    isig3sq = isig3 * isig3;
/* Computing 2nd power */
    sigsq2 = ((double)isig) * ((double)isig) * 2.f;
    value = 0.f;
/* 	 write(*,*)' in GAUSS' */
/* write(*,*)' isig,ix,iy,value = ',isig,ix,iy,value */
/* 	pause */

/* Gaussian source profiles with width  isig: */
    normfac = 0.f;
    int count = 0;
    for (i1 = -isig3; i1 <= isig3; ++i1) {
	for (i2 = -isig3; i2 <= isig3; ++i2) {

/* 	            dist2 = i2*i2 + i1*i1 */

	    dist2 = pow(i1 - delx ,2) + pow(i2 - dely ,2);
	    dist = sqrt(dist2);
	    if (isig != 0) {
		if (dist2 <= (double) isig3sq) {
//		    factorex = exp(-dist2 / sigsq2);
		    factorex = ((p_ld + 1)/(M_PI*R2)) * H(1 - dist2/R2) * pow(1 - dist2/R2, p_ld);
/* 	                   factorex = 1.0 */
		    dum = pix[ix + i1 + (iy + i2) * pix_dim1];
		    value += dum * factorex;
		    normfac += factorex;
		    count++;
		}
	    } else if (isig == 0) {
		factorex = 1.f;
		dum = pix[ix + i1 + (iy + i2) * pix_dim1];
		value += dum * factorex;
		normfac += factorex;
	    }
	}
    }
//    printf("normfac = %f\n", normfac);
    value /= normfac;
    return value;

}

double power_(int isig, double xxx, double yyy, double* pix, int ipix)
{
    double p = static_cast<double>(3)/static_cast<double>(2);
    double R_1_2 = (double)isig * sqrt(log(2));
    double R = R_1_2/sqrt((pow(static_cast<double>(2),static_cast<double>(1)/(p-1))-1)/log(2));
    double R2 = R * R;

    /* System generated locals */
    double value = 0;
    int pix_dim1, pix_offset;

    /* Local variables */
    double factorex;
    int i1, i2, ix, iy;
    double dum, delx, dely;
    int isig3;
    double dist2, sigsq2, dist;
    int isig3sq;
    double normfac;

/* modified: disk of constant brightness */

    /* Parameter adjustments */
    pix_dim1 = ipix;
    pix_offset = 1 + pix_dim1;
    pix -= pix_offset;

    /* Function Body */
    ix = i_dnnt(xxx);
    iy = i_dnnt(yyy);
    delx = xxx - (double) ix;
    dely = yyy - (double) iy;


/* 030711: isig3 = RADIUS  for source!!! */

    isig3 = isig * 10;
/* 	isig3 =   isig */
    isig3sq = isig3 * isig3;
/* Computing 2nd power */
    sigsq2 = ((double)isig) * ((double)isig) * 2.f;
    value = 0.f;
/* 	 write(*,*)' in GAUSS' */
/* write(*,*)' isig,ix,iy,value = ',isig,ix,iy,value */
/* 	pause */

/* Gaussian source profiles with width  isig: */
    normfac = 0.f;
    int count = 0;
    for (i1 = -isig3; i1 <= isig3; ++i1) {
	for (i2 = -isig3; i2 <= isig3; ++i2) {

/* 	            dist2 = i2*i2 + i1*i1 */

	    dist2 = pow(i1 - delx ,2) + pow(i2 - dely ,2);
	    dist = sqrt(dist2);
	    if (isig != 0) {
		if (dist2 <= (double) isig3sq) {
//		    factorex = exp(-dist2 / sigsq2);
		    factorex = ((p - 1)/(M_PI*R2))*(1/pow(1 + dist2/R2, p));
/* 	                   factorex = 1.0 */
		    dum = pix[ix + i1 + (iy + i2) * pix_dim1];
		    value += dum * factorex;
		    normfac += factorex;
		    count++;
		}
	    } else if (isig == 0) {
		factorex = 1.f;
		dum = pix[ix + i1 + (iy + i2) * pix_dim1];
		value += dum * factorex;
		normfac += factorex;
	    }
	}
    }
//    printf("normfac = %f\n", normfac);
    value /= normfac;
    return value;
}


double ad_(int isig, double xxx, double yyy, double* pix, int ipix)
{
    double R_1_2 = (double)isig * sqrt(log(2));
    double R = R_1_2/static_cast<double>(4);
    double R2 = R * R;

    /* System generated locals */
    double value = 0;
    int pix_dim1, pix_offset;

    /* Local variables */
    double factorex;
    int i1, i2, ix, iy;
    double dum, delx, dely;
    int isig3;
    double dist2, sigsq2, dist;
    int isig3sq;
    double normfac;

/* modified: disk of constant brightness */

    /* Parameter adjustments */
    pix_dim1 = ipix;
    pix_offset = 1 + pix_dim1;
    pix -= pix_offset;

    /* Function Body */
    ix = i_dnnt(xxx);
    iy = i_dnnt(yyy);
    delx = xxx - (double) ix;
    dely = yyy - (double) iy;


/* 030711: isig3 = RADIUS  for source!!! */

    isig3 = isig * 10;
/* 	isig3 =   isig */
    isig3sq = isig3 * isig3;
/* Computing 2nd power */
    sigsq2 = ((double)isig) * ((double)isig) * 2.f;
    value = 0.f;
/* 	 write(*,*)' in GAUSS' */
/* write(*,*)' isig,ix,iy,value = ',isig,ix,iy,value */
/* 	pause */

/* Gaussian source profiles with width  isig: */
    normfac = 0.f;
    int count = 0;
    for (i1 = -isig3; i1 <= isig3; ++i1) {
	for (i2 = -isig3; i2 <= isig3; ++i2) {

/* 	            dist2 = i2*i2 + i1*i1 */

	    dist2 = pow(i1 - delx ,2) + pow(i2 - dely ,2);
	    dist = sqrt(dist2);
	    if (isig != 0) {
		if (dist2 <= (double) isig3sq) {
//		    factorex = exp(-dist2 / sigsq2);
                    if (H(dist - R) > 0) factorex = (3 * R  / (2 * M_PI * pow(dist, static_cast<double>(3)))) * (1 - sqrt(R/dist));
/* 	                   factorex = 1.0 */
		    dum = pix[ix + i1 + (iy + i2) * pix_dim1];
		    value += dum * factorex;
		    normfac += factorex;
		    count++;
		}
	    } else if (isig == 0) {
		factorex = 1.f;
		dum = pix[ix + i1 + (iy + i2) * pix_dim1];
		value += dum * factorex;
		normfac += factorex;
	    }
	}
    }
//    printf("normfac = %f\n", normfac);
    value /= normfac;
    return value;
}

#endif