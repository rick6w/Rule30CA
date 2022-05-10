//Original C code by George Marsaglia, modified slightly for use with Rule 30 generator.

/* This file contains a sample main() for calling the
   three tough tests gcd,bday and gorilla, procs for which
   are included.   To run the tests on a RNG of your choice,
   you must include a proc for that generator, and tell
   the three test routines how to access it by means of
   a #define.   A sample RNG is included (one of the versions
   of the KISS generator). Just replace that proc with your own
   and change the #define IUNI accordingly.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/*  USER: IN THE FOLLOWING LINE, CHANGE Kiss() TO THE NAME OF THE (32-BIT) TEST RNG. */ 
#define IUNI rule30_wrap()

static double U[98];
/* The table used by the (97, 33) lagged-Fibonacci RNG
   This table is modified by rstart() and uni()       
   rstart() must be called prior to the calls of uni().
  It sets up the table needed by the lagged-Fibonacci RNG.
  The seeds, i, j, k, l must be in [1, 168].
*/

void rstart(int i,int j,int k,int l)
{
    int ii, m, n;
    double s, t;

    for (ii=1; ii<=97; ++ii)
    {   /*Fill entries of the table U[] */
        s = 0.0;
        t = 0.5;
        for (n=1; n<=52; ++n)
        {    /*Generate seeds bit by bit */
            m = (((i*j)%179)*k)%179;   /*First mini RNG */
            i = j;
            j = k;
            k = m;
            l = (53*l+1)%169; /*Second mini RNG */
            if (((l*m)%64)>=32)
                s = s+t;  /* Generate a bit */
            t = 0.5*t;
        }
        U[ii] = s;
    }
}  /* end rstart() */


//Rule 30 generator, wraps lsb and msb around.
unsigned int rule30_wrap()
{
  static unsigned int x = 0x00010000; //32-bit seed
  x = ((x>>1)^(x<<31)) ^ (x | ((x<<1)^(x>>31)));
  return x;
}


double chisq(double d, double z)
{
    double v,t;
    if (z<=0.) return(0.);
    t = (exp(log(z/d)/3.)-1.+2/(9*d))/sqrt((float)(2/(9*d)));
    if (t>5.) return (1.);
    if (t+5.<0.) return(0.);
    v = (t>0.) ? 1./(1.+.33267*t) : 1./(1.-.33267*t);
    v = exp(-.5*t*t)*(.1740120799+(.3739278012*v-.04793993633)*v)*v;
    if(t>0) v = 1.-v;
    return(v);
}

void gcd(unsigned int n)
{ //  gcd Test, n pairs of IUNI's
    double s,ch32,ch99,e,r=1.0e-10;
    unsigned int u, v, w, i, k, t[47], gcd[101];
    unsigned int kp[36]= 
    { 0,0,0,5506,29532,144541,590691,2065333,6277273,16797273,39964829,
      85160313,163520964,284315128,449367802,647663189,853365977,1030017639,
      1140689999,1160424679,1085307890,933657875,738971259,538010076,360199303,
      221583256,125137422,64787412,30697488,13285746,5238570,1876493,608568,
      177920,46632,13674
    };
    for (i=1; i<101; i++)
    {
        gcd[i] = 0;
        if (i<47) t[i] = 0;
    }
    for (i=1; i<=n; i++)
    {
        k = 0;
        do
        {u = IUNI;v = IUNI;}
        while (u==0 || v==0);
        do
        {
            w = u%v;
            u = v;
            v = w;
            k++;
        } while(v>0);
        if (u>100) u = 100;
        gcd[u]++;
        if (k<3) k = 3;
        if (k>35) k = 35;
        t[k]++;
    }
    ch32 = 0.0;
    printf("Euclid's algorithm:\n");
    for(i=3; i<=35; i++)
    {
        e = (r*kp[i])*n;
        s = (t[i]-e)*(t[i]-e)/e;
        ch32 += s;
    }

// User: you may want to insert printf("%4d %6.3f %7.2f\n",i,s,ch32);
// in the above loop to see the (O-E)^2/E values and their sums
    printf(    " p-value, steps to gcd:   %f\n",chisq(16.,.5*ch32));

    e = n*.61097691e-2;
    ch99 = (gcd[100]-e)*(gcd[100]-e)/e;
    for (i=1; i<100; i++)
    {
        e = n*0.6079271/(i*i);
        s = (gcd[i]-e)*(gcd[i]-e)/e;
        ch99 += s;
    }
// User: you may want to insert printf("%4d %6.3f %7.2f\n",i,s,ch99);
// in the above loop to see the (O-E)^2/E values and their sums

   printf(" p-value, dist. of gcd's: %f\n",chisq(49.5,.5*ch99));
   printf("\n");
}

// requires #define IUNI rngf() // with your rngf() function


int compi(const void *i, const void *j) // compare function for qsort
{
    if  ( *(unsigned int *)i <  *(unsigned int *)j )
        return -1;
    if  ( *(unsigned int *)i > *(unsigned int *)j  )
        return 1;
    return 0;
}


void bday(void)
{
    int i, j, k, m=4096;
    unsigned int t[4096],obs[11]={0,0,0,0,0,0,0,0,0,0,0};
    double w, x, ex[11], v=0.;
    printf("Birthday spacings test: 4096 birthdays, 2^32 days in year\n");
    for (k=1; k<=5000; k++)
    {//loop for sample of size 5000
        for (i=0; i<m; i++)
            t[i] = IUNI;    //   get 4096 IUNI's
        qsort(t,m,sizeof(int),compi); //  sort them
        for (i=m-1; i>0; i--)
            t[i] = t[i]-t[i-1]; // get spacings
        qsort(t,m,sizeof(int),compi);  // sort spacings
        j = 0;
        for (i=1; i<m; i++)
            j = j+(t[i]==t[i-1]) ;
        if(j>10) j = 10;
        obs[j]++;      // count duplicate spacings
    }// end  sample of 5000 loop
    printf("           Table of Expected vs. Observed counts:\n");
    ex[0] = 5000.*exp(-4.); ex[10]=5000.*.008132242;
    for(i=1; i<10; i++)
        ex[i] = 4.*ex[i-1]/i;
    printf("Duplicates   0");
    for (i=1; i<10; i++)
        printf("%6.0f",i+0.);
    printf("   >=10\n");
    printf("Expected ");
    for(i=0;i<10;i++)
        printf("%6.1f",ex[i]);
    printf("%6.1f\n",ex[10]);
    printf("%8s","Observed");
    for (i=0; i<10; i++)
        printf("%6.0f",obs[i]+0.0);
    printf("%6.0f\n ",obs[10]+0.0);
    printf("\b%9s","(O-E)^2/E");
    for (i=0; i<11; i++)
    {
        x = (obs[i]-ex[i])*(obs[i]-ex[i])/ex[i];
        v += x;
        printf("%6.1f",x);
    }
    x = .5*v;
    w = 24.+(24.+(12.+(4.+x)*x)*x)*x;
    w = 1.-exp(-x)*w/24.;
    printf("\n              Birthday Spacings: Sum(O-E)^2/E=%7.3f, p=%6.3f\n",v,w);
    printf("\n");
} // end bday

float ad32(float z)
{  // converts AD KS statistic into p-value, n=32
    float y;
    if (z<.5)
    {
        y = exp(-1.27*log(z));
        return( exp(.266-(.70025-(.009804-.000213*y)*y)*y));
    }
    if (z<1.)
        return ( (.53383+(1.65285-(1.988-.6634*z)*z)*z)*z-.21862 );
    if (z<2.)
        return ( .99987-.6616*exp(-1.0896*z)-.953*exp(-2.005*z)  );
    return( 1-.52686*exp(-1.05276*z)-.68053*exp(-1.62034*z) );
}

void gorilla(void)
{// begin gorilla
    unsigned int *t;
    float z,v,u[32];
    int s[256]={0}; // table for quick count-the-1's in a long stream
    unsigned int k,w,i,j,q,mask,m[32];
    int n=67108864;
    printf(" Gorilla test for 2^26 bits, positions 0 to 31:\n");
    printf(" Note: lengthy test---for example, ~20 minutes for 850MHz PC\n");
    for (i=1; i<256; i++)
    {
        q = i;
        s[i] = 1;
        while (q &= (q-1))
            s[i]++;
    } // set s table
    m[31] = 1;
    for (i=31; i>0; i--)
        m[i-1] = 2*m[i]; // set masks
    t = (unsigned int *) malloc(2097152*sizeof(unsigned int));
    if (t == NULL)
    {
        printf("Not enough memory for gorilla test\n");
        exit(1);
    }
    for (k=0; k<32; k++)
    { //k loops thru bits 0,1,..31
        mask = m[k];     // mask set to bit k
        for (i=0; i<2097152; i++)
             t[i] = 0;  // initialize t-array
        // fill the initial word with 25 bits
        w = 0;
        for (i=0; i<25; i++)
        {
            w = (w>>1);
            if (IUNI&mask)
                w = w+(1<<25);
        }
        for (i=0; i<n; i++)
        { //i loop, gets 2^26 w's, sets bit in t[] for each
            w = (w>>1);
            if (IUNI&mask)
                w = w+(1<<25);
            j = (w&2097151);
            t[j] = (t[j] | m[w>>21]);
        } //end i loop
        w=0;
        for (i=0; i<2097152; i++)
        {
            j = t[i]; // loop to count 1's in t[] array
            w += s[(j)&255]+s[(j>>8)&255]+s[(j>>16)&255]+s[(j>>24)&255];
        }        // end count loop
        z = ((n-w)-24687971.)/4170.;    // convert counts to normal z-score
        v = (z>0.) ? 1./(1.+.33267*z) : 1./(1.-.33267*z);
        v = exp(-.5*z*z)*(.1740120799+(.3739278012*v-.04793993633)*v)*v;
        if (z>0)
           v = 1.-v;
        u[k] = v;       // z-score converted to p-value
        if (k%8 == 0)
           printf(" Bits %2u to %2u--->",k,k+7);
        printf( "%6.3f", v);
        if (k%8 == 7) printf("\n");
    } //end k loop
    free(t);
    // Now do KS test on 32 u's.  first, sort u's.
    for (i=0; i<31; i++)
    {
        for (j=0; j<31-i; j++)
           if (u[j+1] < u[j])
           {
                v = u[j];
                u[j] = u[j+1];
                u[j+1] = v;
            }
    }
    // calculate Anderson-Darling statistic
    z = -1024.;
    for(i=0; i<32; i++)
    {
        v = u[i]*(1.-u[31-i]);
        if (v<1.e-30)
           v = 1.e-30;
        z = z-(i+i+1)*log((double)v);
    }
    printf(" KS test for the above 32 p values: %6.3f\n",ad32(z/32.));
}//end gorilla

int main(void)
{
    rstart(23,55,76,91);
    bday();
    gcd(10000000);
    gorilla();
    return 0;
}


