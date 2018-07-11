/*function to get the ith column of the modularity matrix for a multislice bipartite network

Version 1.0, January 27, 2012

Lucas Jeub


b=modf_bipord(i,m,n,T,A,k,d,mm,omega)

*/


#include "mex.h"
#include "matrix.h"
#include "math.h"

/*function*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	/*check for correct input (not exhaustive)*/
	
	/*#input arguments*/
	if (nrhs!= 9){
		mexErrMsgTxt("9 input arguments required");
	}
	
	/*#output arguments*/
	if (nlhs != 1){
		mexErrMsgTxt("1 output argument required");
	}
	
	/*check for sparse adjacency matrix*/
	if(!mxIsCell(prhs[4])){
		mexErrMsgTxt("Needs Cell adjacency matrix");
	}
		
	
	/*declare variables*/
	int i, ii,ir,s, N, m, n, T,nnz, count;
	double  *value, *valuemod, *modv, *kv, *dv, *mm, omega, v, kk, dd;
	mwIndex *indrow, *indcol, *indrowmod, *indcolmod, *kr,*kc,*dr,*dc;
	double *Acol;
	mxArray* As;
	
	/*initialise*/
	ir= mxGetPr(prhs[0])[0];
	
	m= mxGetPr(prhs[1])[0];
	n= mxGetPr(prhs[2])[0];
	T= mxGetPr(prhs[3])[0];
	
	kv=mxGetPr(prhs[5]);
	kr=mxGetIr(prhs[5]);
	kc=mxGetJc(prhs[5]);



	dv=mxGetPr(prhs[6]);
	dr=mxGetIr(prhs[6]);
	dc=mxGetJc(prhs[6]);

	mm=mxGetPr(prhs[7]);
	omega=mxGetPr(prhs[8])[0];
	
	N=m+n;
	s=ceil(((double) ir)/N)-1;
	ii=ir-(s)*N;
	
	As=mxGetCell(prhs[4],s);
	value=mxGetPr(As);
	indrow=mxGetIr(As);
	indcol=mxGetJc(As);
	nnz=mxGetNzmax(As);
	
			

		
	if (ii<=m){
		
		Acol=mxCalloc(n,sizeof(double)); 	/*array for the ith row of the adjacency matrix*/
		modv=mxCalloc(n,sizeof(double));	/*array for null-model contribution*/	
		
		plhs[0]=mxCreateSparse(N*T,1,n+2,0); /*create output matrix*/
	
		valuemod=mxGetPr(plhs[0]);	/*get array for values*/
		indrowmod=mxGetIr(plhs[0]);	/*get array for row indeces*/		
		indcolmod=mxGetJc(plhs[0]);	/*get array for column indeces*/
		
		/*k(ii,s)*/
		kk=0;
		for(i=kc[s];i<kc[s+1];i++){
			if(kr[i]==ii-1){
				kk=kv[i];
				break;	
			}
		}
		
		
		for(i=dc[s];i<dc[s+1];i++){
			modv[dr[i]]=kk*dv[i]/mm[s];
		}

		/*get values of the ith row of the adjacency matrix*/
	count=0;	
		for(i=0;i<nnz;i++){
			while(i==indcol[count+1]){
				count++;
				}
			if(indrow[i]==ii-1){
				Acol[count]=value[i];
			}
		}
						
				
	count=0; /*keep track of number of non-zero elements */
		
		/*interslice coupling*/
		if((-N+ir)>0){
		indrowmod[0]=-N+ir-1;
		valuemod[0]=omega;
		count++;
		}
		
		/*intraslice coupling*/
		for(i=0;i<n;i++){
			
			v=Acol[i]-modv[i]; /*calculate modularity*/
			
			if(v!=0){
			indrowmod[count]=m+i+(s)*N;
			valuemod[count]=v;
			count++;
			}
		}
		
		/*interslice coupling*/
		if((N+ir)<=N*T){
		indrowmod[count]=N+ir-1;
		valuemod[count]=omega;
		count++;
		}
		
		/*Asign Column Indeces (indcol[j] is the index of the 
		first non-zero element in the jth row (starting at zero) */

		indcolmod[0]=0;
		indcolmod[1]=count;
		
		
		
				
				
	}
	else{
		
		Acol=mxCalloc(m, sizeof(double)); /*array for ith column of adjacency matrix*/
		modv=mxCalloc(m, sizeof(double));
		
		plhs[0]=mxCreateSparse(N*T,1,m+2,0); /*create output matrix*/
		
		valuemod=mxGetPr(plhs[0]);	/*get array for values*/	
		indrowmod=mxGetIr(plhs[0]);	/*get array for row indeces*/
		indcolmod=mxGetJc(plhs[0]);	/*get array for column indeces*/

		dd=0;
		for(i=dc[s];i<dc[s+1];i++){
			if(dr[i]==ii-1-m){
				dd=dv[i];
				break;
			}
		}
		
		for(i=kc[s];i<kc[s+1];i++){
			modv[kr[i]]=dd*kv[i]/mm[s];
		}
		
		/*get values for the ith column of the adjacency matrix*/
		for(i=indcol[ii-m-1];i<indcol[ii-m];i++){
			Acol[indrow[i]]=value[i];
		}				

	
			
		
	count=0;	/*keep track of number of non-zero elements*/
	
		/*interslice coupling*/
		if((-N+ir)>0){
		indrowmod[0]=-N+ir-1;
		valuemod[0]=omega;
		count++;
		}
			
		/*intraslice coupling*/
		for (i=0; i<m; i++){
			v=Acol[i]-modv[i]; /*calculate modularity*/
			if(v!=0){
				valuemod[count]=v;
				indrowmod[count]=i+(s)*N;
				count++;
			}
		}
		
		/*interslice coupling */
		if((N+ir)<=N*T){
		indrowmod[count]=N+ir-1;
		valuemod[count]=omega;
		count++;
		}
		
		/*Asign Column Indeces (indcol[j] is the index of the 
		first non-zero element in the jth row (starting at zero) */
		
		indcolmod[0]=0;
		indcolmod[1]=count;		
		
				
	}
	
	
/* Free Memory*/	
mxFree(Acol);
mxFree(modv);


}


