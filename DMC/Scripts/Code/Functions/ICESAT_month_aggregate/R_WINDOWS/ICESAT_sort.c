#include <R.h> 		 /* required  */
#include <Rmath.h>  	 /* for distribution functions etc. */  

/* Averaging the ICESAT data on the 1km grid*/


void ICESAT_sort(int *x, int *y, double *z, int *n, int *x2, int *y2, double *z2, int *n2)

{
	int j, count,buffer;
	
	j = 0;
	count = 0;
	buffer = 1;
		
	while (j < *n) {
		if ((x[j] == x2[count]) & (y[j] == y2[count])) {
 			buffer++;
			z2[count] = z2[count] + z[j];
			j++;
		} else {
			if (buffer > 1) {
					z2[count] = z2[count]/buffer;
					buffer = 1;
					count++;
			}
			x2[count] = x[j];
			y2[count] = y[j];
			z2[count] = z[j];
			j++;
			
		}
	}
	*n2 = count;
}

