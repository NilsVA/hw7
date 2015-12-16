#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void function(double* kx,double* ky,double* x,double* y,double my)
{       
        double r = sqrt(pow((x[0]+my),2)+pow(y[0],2));
        double s = sqrt(pow(x[0]-1-my,2)+pow(y[0],2));
	kx[0]=x[1];
        kx[1]=x[0]+2*y[1]-(1-my)*(x[0]+my)/pow(r,3)-my*(x[0]-1+my)/pow(s,3);
        
        ky[0]=y[1];
        ky[1]=y[0]-2*x[1]-(1-my)*y[0]/pow(r,3)-my*y[0]*pow(s,3);
}


void rk45(double* x,double* y,double* xrk4,double* xrk5,double* yrk4,double* yrk5,double dt, double my);

int main()
{
	
	const double tmax=100.0;
        double x[2];
        double y[2];
        double xrk4[2];xrk4[0]=1;xrk4[1]=1;
        double yrk4[2];yrk4[0]=1;yrk4[1]=1;
        double xrk5[2];xrk5[0]=0;xrk5[1]=0;
        double yrk5[2];yrk5[0]=0;yrk5[1]=0;
        double kx[2];
        double ky[2];
        double tol = 1e-5;
	double dt=1e-3;
	int n = 0;
	
	
	x[0]=0.994;x[1]=0;
        y[0]=0;y[1]=(-1)*(2.00158510637908);
        double my=0.012277471;
        ofstream out("data.txt");
        while (abs(sqrt(pow(xrk4[0]+xrk5[0],2)+pow(xrk4[1]+xrk5[1],2)))>tol && n < 1000)
        {
            rk45(x,y,xrk4,xrk5,yrk4,yrk5,dt,my);
            
            dt= 0.9*dt*tol/abs(sqrt(pow(xrk4[0]+xrk5[0],2)+pow(xrk4[1]+xrk5[1],2)));
            
            int i=0;
            
            out << i*dt << "\t" << dt << endl;
            i++;
            n++;
        }
        out.close();
        
	
}

void rk45(double* x,double* y,double* xrk4,double* xrk5,double* yrk4,double* yrk5,double dt, double my)
{       int dim = 2;
	double kx1[dim];double ky1[dim];double kx2[dim];double ky2[dim];
        double kx3[dim];double ky3[dim];double kx4[dim];double ky4[dim];
        double kx5[dim];double ky5[dim];double kx6[dim];double ky6[dim];
        double kx7[dim];double ky7[dim];
	double ftempx[dim];
	double ftempy[dim];
        
	
	function(kx1,ky1,x,y,my);
        
        for(int i=0;i<dim;i++){
        ftempx[i]=x[i]+(dt/5.)*kx1[i]/5.;        
        ftempy[i]=y[i]+(dt/5.)*ky1[i]/5.; 
        }
        
        function(kx2,ky2,ftempx,ftempy,my);
        
        for(int i=0;i<dim;i++){
        ftempx[i]=x[i]+(3.*dt/10.)*(3.*kx1[i]/40. + 9./40.*kx2[i]); 
        ftempy[i]=y[i]+(3.*dt/10.)*(3.*ky1[i]/40. + 9./40.*ky2[i]); 
        }
        
        function(kx3,ky3,ftempx,ftempy,my);  
        
        for(int i=0;i<dim;i++){
        ftempx[i]=x[i]+(4.*dt/5.)*(44.*kx1[i]/45. - 56./15. * kx2[i] + 32./9. * kx3[i]); 
        ftempy[0]=y[i]+(4.*dt/5.)*(44.*ky1[i]/45. - 56./15. * ky2[i] + 32./9. * ky3[i]); 
        }
        
        function(kx4,kx5,ftempx,ftempy,my);
        
        for(int i=0;i<dim;i++){
        ftempx[i]=x[i]+(8.*dt/9.)*(19372.*kx1[i]/6561. - 25360./2187. * kx2[i] + 64448./6561. * kx3[i] - 212./729. * kx4[i]); 
	ftempy[i]=y[i]+(8.*dt/9.)*(19372.*ky1[i]/6561. - 25360./2187. * ky2[i] + 64448./6561. * ky3[i] - 212./729. * ky4[i]);
        }
        
        function(kx5,ky5,ftempx,ftempy,my);
        
        for(int i=0;i<dim;i++){
        ftempx[i]=x[i]+(dt)*(9017.*kx1[i]/3168. - 355./33. * kx2[i] + 46732./5247. * kx3[i] - 49./176. * kx4[i] - 5103./18656. * kx5[i]); 
	ftempy[i]=y[i]+(dt)*(9017.*ky1[i]/3168. - 355./33. * ky2[i] + 46732./5247. * ky3[i] - 49./176. * ky4[i] - 5103./18656. * ky5[i]); 
        }
        
        function(kx6,ky6,ftempx,ftempy,my);
        
        for(int i=0;i<dim;i++){
        ftempx[i]=x[i]+(dt)*(35.*kx1[i]/384. + 500./1113. * kx3[i] + 125./192. * kx4[i] - 2187./6784. * kx5[i] + 11./84. * kx6[i]); 
	ftempy[i]=y[i]+(dt)*(35.*ky1[i]/384. + 500./1113. * ky3[i] + 125./192. * ky4[i] - 2187./6784. * ky5[i] + 11./84. * ky6[i]); 
        }
        
        function(kx7,ky7,ftempx,ftempy,my);
        
        for(int i=0;i<dim;i++){
	xrk4[i]=x[i]+(dt)*(kx1[i]* 5179./57600. + 7571./16695. * kx3[i] + 393./640.* kx4[i] - 92097/339200* kx5[i] + 187./2100. * kx6[i] + 1./40.* kx7[i]);
        yrk4[i]=y[i]+(dt)*(ky1[i]* 5179./57600. + 7571./16695. * ky3[i] + 393./640.* ky4[i] - 92097/339200* ky5[i] + 187./2100. * ky6[i] + 1./40.* ky7[i]);
        }
	
	for(int i=0;i<dim;i++){
        xrk5[i]=x[i]+(dt)*(kx1[i]* 35./384. + 500./1113. * kx3[i] + 125./192.* kx4[i] - 2187/6784* kx5[i] + 11./84. * kx6[i]);
	yrk5[i]=y[i]+(dt)*(ky1[i]* 35./384. + 500./1113. * ky3[i] + 125./192.* ky4[i] - 2187/6784* ky5[i] + 11./84. * ky6[i]);
        }
	
	
	
}
