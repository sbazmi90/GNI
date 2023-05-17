//Geometric integration of a point vortex system.(Symplectic Euler method)//

#include <stdlib.h>
#include <stdio.h>
//#include <conio.h>
#include <math.h>
#define N 9
#define h 1
#define G 0.07957747154

main()
{
  // files and parameters  
  FILE *af0,*af1,*af2,*af3,*af4,*af5,*af6,*af7,*af8,*af9,*af10;
  af0= fopen("data.txt", "r");
  af1=fopen("Object1.txt","w");
  af2=fopen("Object2.txt","w");
  af3=fopen("Object3.txt","w");
  af4=fopen("Object4.txt","w");
  af5=fopen("Object5.txt","w");
  af6=fopen("Object6.txt","w");
  af7=fopen("Object7.txt","w");
  af8=fopen("Object8.txt","w");
  af9=fopen("Object9.txt","w");
 
  af10=fopen("Energy.txt","w");

  int i,j,t,timestep; 
  double   x[N],y[N],kx1[N],ky1[N],gamma[N];     
  double  r,r2,px,py,lz,H,dx,dy;

  //Read initial states from data file//
  if (af0 == NULL)
    {
      printf("Error: cannot open file for reading...\n");
      exit (0);
    }
  else {
    i=1;
    while(fscanf(af0,  "%lf  %lf %lf", &gamma[i] , &x[i], &y[i]) != EOF)
      i++;				
  }

  for(i=1;i<=N;i++)
    {
      if(i==1) fprintf(af1,"0  %20.18f %20.18f \n",x[i],y[i]);
      if(i==2) fprintf(af2,"0  %20.18f %20.18f \n",x[i],y[i]);
      if(i==3) fprintf(af3,"0  %20.18f %20.18f \n",x[i],y[i]);
      if(i==4) fprintf(af4,"0  %20.18f %20.18f \n",x[i],y[i]);
      if(i==5) fprintf(af5,"0  %20.18f %20.18f \n",x[i],y[i]);
      if(i==6) fprintf(af6,"0  %20.18f %20.18f \n",x[i],y[i]);
      if(i==7) fprintf(af7,"0  %20.18f %20.18f \n",x[i],y[i]);
      if(i==8) fprintf(af8,"0  %20.18f %20.18f \n",x[i],y[i]);
      if(i==9) fprintf(af9,"0  %20.18f %20.18f \n",x[i],y[i]);
    }

  timestep = 15552000;

  //Start time evolution in the system//

  for(t=1; t<timestep ;t++)
    {
      px=0,py=0,lz=0,H=0;
  
      for(i=1;i<=N;i++)
	{
	  kx1[i]=0;
	  ky1[i]=0;      
	}
      
      for(i=1;i<=N;i++)
	{              
	  for(j=1;j<=N;j++)
	    {
	      if(i!=j)
		{
		  dx=x[i]-x[j]; 
		  dy=y[i]-y[j];
 
		  r=sqrt(dx*dx+dy*dy);
		  r2=r*r;
		  kx1[i]=kx1[i]-G*gamma[j]*(y[i]-y[j])/r2;
		  // ky1[i]=ky1[i]+G*(x[i]-x[j])/r2;	   
		  H=H-G*gamma[i]*gamma[j]*log(r);
		}//end of if
	    }//end of j
	  lz=lz-0.5*gamma[i]*(x[i]*x[i]+y[i]*y[i]);
	  px= px+gamma[i]*y[i];
	  py= py-gamma[i]*x[i];       
	}//end of i

      //////////////////////////////////////////////////////////

      for(i=1;i<=N;i++)
	{
	  x[i]=x[i]+h*kx1[i];
	}
      
      for(i=1;i<=N;i++)
	{              
	  for(j=1;j<=N;j++)
	    {
	      if(i!=j)
		{
		  dx=x[i]-x[j]; 
		  dy=y[i]-y[j];
 
		  r=sqrt(dx*dx+dy*dy);
		  r2=r*r;
		  ky1[i]=ky1[i]+G*gamma[j]*(x[i]-x[j])/r2;
	   
		}//end of if
	    }//end of j       
	}//end of i
      ///////////////////////////////////////////////
      for(i=1;i<=N;i++)
	{
	  y[i]=y[i]+h*ky1[i];
	}
      if(t%7200==0)
	{
	  for(i=1;i<=N;i++)
	    {  
	      if(i==1) fprintf(af1,"%d  %20.18f %20.18f \n",t/7200,x[i],y[i]);
	      if(i==2) fprintf(af2,"%d  %20.18f %20.18f \n",t/7200,x[i],y[i]);
	      if(i==3) fprintf(af3,"%d  %20.18f %20.18f \n",t/7200,x[i],y[i]);
	      if(i==4) fprintf(af4,"%d  %20.18f %20.18f \n",t/7200,x[i],y[i]);
	      if(i==5) fprintf(af5,"%d  %20.18f %20.18f \n",t/7200,x[i],y[i]);
	      if(i==6) fprintf(af6,"%d  %20.18f %20.18f \n",t/7200,x[i],y[i]);
	      if(i==7) fprintf(af7,"%d  %20.18f %20.18f \n",t/7200,x[i],y[i]);
	      if(i==8) fprintf(af8,"%d  %20.18f %20.18f \n",t/7200,x[i],y[i]);
	      if(i==9) fprintf(af9,"%d  %20.18f %20.18f \n",t/7200,x[i],y[i]);
	    }

	  fprintf(af10,"%d  %20.18f %20.18f %20.18f %20.18f \n",t/7200, px, py, lz, H);
	}
        
    }// end of time step
 
  fclose(af0);
  fclose (af1);
  fclose (af2);
  fclose (af3);
  fclose (af4);
  fclose (af5);
  fclose (af6);
  fclose (af7);
  fclose (af8);
  fclose (af9);
  fclose (af10);
  
 
  return 0; 
}
