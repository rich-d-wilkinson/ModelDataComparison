/* Call this from R using 
   R CMD SHLIB Test.c
   Once compiled load the library in R using dyn.load("Test.so")
   Run as .C("hello",as.integer(n))
   Like this all the functions here are also in R */

#include <R.h> 		 /* required  */
#include <Rmath.h>  	 /* for distribution functions etc. */  

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define verbose 0

void hello(int *n)
{
	int i;
	for(i=0;i < *n; i++) {
	   Rprintf("Hello, world!\n");
	}
}



void kernel_smooth(double *x, int *n, double *xpts, int *nxpts, double *h, double *result)
{
	int i,j;
	double d, ksum;
	for(i=0; i<*nxpts;i++) {
		ksum =0;
		for(j=0;j<*n;j++) {
			d = xpts[i] - x[j];
			ksum += dnorm(d / *h,0, 1, 0);
		}
		result[i] = ksum / ((*n) * (*h));
	}
}

void TakahashiC(int *n, int *N, int *ii,int *jj, double *x, double *y)
{
   int j,i,t,t2,count,kcount,skjcount,jscount,jscount2;
   float Lii;
   float skj[*n],Lki[*n],Lki2[*n],k[*n], k2[*n], js[*n],js2[*n];




   count = *N;
   *(y+count-1) = 1 / pow(*(x+count-1),2);     
   count--;

   if (*n > 1) 
   {
	
         for (i=(*n-1); i >=1;i--) 
	  {
              
 	     /* Find row numbers with non-zero indices at this column */
	     kcount = 0;
	      Rprintf("i = %d\n",i);

	     /* While still not on the diagonal */
	     /*            k <- ii[ii > i & jj == i]  */
	     while((*(jj+count-1) == i)&(*(ii+count-1) > i)){
		   k[kcount] = *(ii + count - 1);
		   kcount++;
		   count--;
		  }
 		
	   
	    if (kcount == 0) {
 		  *(y+count-1) = 1 / pow(*(x+count-1),2); 
		   
	     } else {
  	    /*  Lii <- L[i, i]
                Lki <- L[k, i]  */
	        jscount = 0;
		Lii = *(x+count-1);
		if(verbose) Rprintf("i = %d: Lii = %lf \n",i,Lii);
		for(t=1; t <= kcount;t++) {
		  Lki[t-1]= *(x+count+kcount-t);
		 }	

		/* Reverse k and Lki */
		for(t=0; t < kcount;t++) {
		  Lki2[kcount-t-1] = Lki[t];
		  k2[kcount-t-1] = k[t];
		}
		for(t=0; t < kcount; t++) {
		 Lki[t] = Lki2[t];
		 k[t] = k2[t];
		 if(verbose) Rprintf("k = %d: Lki = %lf \n",k[t],Lki[t]);
		}
				
		/* js <- rev(jj[ii %in% k & jj >= i]) # going backwards */
	/*	for(t2 = 0; t2 <= (*N - count ); t2++) {
		for(t=1; t <= kcount;t++) {

		     if (ii[t2+count-1] == k[t-1])  {
			  js[jscount] = jj[t2+count-1];
			  // Rprintf("js = %d: \n",js[jscount]);
			  jscount++;
		      }
		  }
		  }

	       // Find unique js 
 		if (jscount > 1) {
			jscount2 = 1;
			js2[0] = js[0];
			for (t = 2; t <= jscount; t++) {
				if (js[t-1] != js[t-2]) {
					js2[jscount2] = js[t-1];
					jscount2++;
				}
   		        }
			jscount = jscount2;
			for  (t = 1; t <= jscount2; t++) {
				js[t-1] = js2[t-1];
			}
			jscount = jscount2;
		}
		
	*/

	for(t=0; t < kcount; t++) {
		js[jscount+1] = k[t];
		if(verbose) Rprintf("js = %d: \n",js[jscount]);
		jscount++;
		
	} 
	js[0] = i;
	jscount++;
	

	 /*skj <- apply(matrix(k),1,function(i){ S@x[ii == max(i,j) & jj == min(j,i)] } )*/
		for(j=jscount-1; j>=0; j--) {
		  if(verbose) Rprintf("Iterating j = %d\n",js[j]);
		  skjcount = 0;
		  for(t = count-1; t < *N; t++){
                  /*if ((*(ii+t) == max(js[j],k[t % kcount])) & (*(jj+t) == min(js[j],k[t % kcount]))) */
			for(t2 = 0; t2 < kcount; t2++) 
		        if ((*(ii+t) == max(js[j],k[t2])) & (*(jj+t) == min(js[j],k[t2])))
		{
			skj[skjcount] = *(y +t);
			if(verbose) Rprintf("skj = %lf \n",skj[skjcount]);
			skjcount++;
		
			}
		  }


	          double mysum =0;
		  for(t=1; t <= kcount;t++) {
			mysum = mysum + skj[t-1]*Lki[t-1];
		  } 
		  mysum = ((i==js[j])/Lii - mysum)/Lii;
		  *(y + count +j -1) = mysum;
		    if(verbose) Rprintf("y[%d] = %lf\n",count+j-1,mysum); 

/*		  skj <- S@x[ii == pmax(k, j) & jj == pmin(k, j)] 
                   S@x[ii == j & jj == i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii */
		  
		}

		
		

	     }
	     count--;
	    
		

/*	     if (*(ii+count-1) == *(jj+count-1))
	      {
		  *(y+count-1) = 1 / pow(*(x+count-1),2); 
	      }   
*/	  } 
  
    } 
/*            k <- ii[ii > i & jj == i]      # Find row numbers with non zero indices at this column
            if (length(k) == 0) 
	    {
                S@x[ii == i & jj == i] <- 1 / L[i, i]^2 
            }   else 
	    {
                Lii <- L[i, i]
                Lki <- L[k, i]
  
                js <- rev(jj[ii %in% k & jj >= i]) # going backwards
                for (j in js) 
	        {
                   skj <- S@x[ii == pmax(k, j) & jj == pmin(k, j)] # select from lower triangle
                   S@x[ii == j & jj == i] <- ((i==j) / Lii - sum(Lki * skj)) / Lii
                }
             }
          }
       } */

}
	



