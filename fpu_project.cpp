//
// 
//  
//
//  Created by Preetha Saha on 4/26/17.
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>

using namespace std;
const double Pi=3.14159;

//fourier coefficients of the position of the masses

vector<double> Fourier_transform(const vector<double>& x){
    vector<double> a(x.size());
    const int N=x.size() - 1;
    for(int k=1;k<N;k++){
        for(int i=1;i<N;i++){
            a[k] += x[i]*sin((i*k*Pi)/double(N));
        }
    }
    return a;
}

//kinetic energy
double Energy_kin(const vector<double>& x, const vector<double>& dxdt, const double alpha){
    const int N=x.size() - 1;
    double E_kin=0;
    for(int i=0;i<N;i++){
        E_kin += (dxdt[i]*dxdt[i])/2.0;
    }
    
    return E_kin;
}

////Potential energy

double Energy_pot(const vector<double>& x, const vector<double>& dxdt,const double alpha){
    
    double E_pot = 0;
    const int N=x.size() - 1;
    for (int i = 1; i <= N; ++i) {
        double dx = x[i] - x[i-1];
        E_pot += pow(dx, 2.0) / 2.0 + alpha * pow(dx, 3.0) / 3.0;
    }
    return E_pot;
    
}

// Local energy distribution amonst the masses, function as given in question part 3
double Energy_local(const int i,const vector<double>& x, const vector<double>& dxdt,const double alpha){
    
    double E_local = 0;
    const int N=x.size() - 1;
    
    double dx_plus = x[i+1] - x[i];
    double dx_minus = x[i] - x[i-1];
    
    E_local = .25*(pow(dx_plus, 2.0) + pow(dx_minus, 2.0)) ;
    
    return E_local;
    
}

//acceleration at every mass point
void cal_accelerations(const vector<double>& x,vector<double>& d2xdt2,const double alpha)
{
    int N = x.size() - 1;
    for (int i = 1; i < N; ++i) {
        double dx_plus = x[i+1] - x[i], dx_minus = x[i] - x[i-1];
        d2xdt2[i] = dx_plus - dx_minus +
        alpha * (pow(dx_plus, 2.0) - pow(dx_minus, 2.0));
    }
}




//starting condition with the purturbation amplitude A=0.2 and sine wave
double initial_condition(vector<double>& x,vector<double>& dxdt,vector<double>& d2xdt2, const double alpha){
    double N=x.size()-1;
    for(int i=0;i<=N;i++){
        x[i]=0.2*sqrt(2.0/N)*sin((i*Pi)/N);// purturbation amplitude 0.2
        
        dxdt[i] = d2xdt2[i] = 0.0;//initial velocities and acceleration is 0
        
    }
    double E_0=Energy_pot(x,dxdt,alpha);
    cal_accelerations(x, d2xdt2, alpha);
    
    return E_0;
    
}


//starting condtion for question part 3

void initial_condition_2(vector<double>& x,vector<double>& dxdt,vector<double>& d2xdt2, const double alpha2){
    double N=x.size()-1;
    const double L=N*1;
    for(int i=0;i<=N;i++){
        x[i]=0;
        
        dxdt[i]=0.8*exp(-(pow(i-L/2,2))/100);
        d2xdt2[i] = 0.0;
        
    }
    cal_accelerations(x, d2xdt2, alpha2);
    
}


//velocity verlet integration algorithm for calculating successive velocities in time

void velocity_Verlet(const double dt,vector<double>& x,vector<double>& dxdt,vector<double>& d2xdt2,const double alpha)
{
    int N = x.size() - 1;
    for (int i = 1; i < N; ++i) {
        x[i] += dxdt[i] * dt + 0.5 * d2xdt2[i] * dt * dt;
        dxdt[i] += 0.5 * d2xdt2[i] * dt;
    }
    cal_accelerations(x, d2xdt2, alpha);
    for (int i = 1; i < N; ++i)
        dxdt[i] += 0.5 * d2xdt2[i] * dt;
}
int main(){
    FILE *f = fopen("E18.txt", "w+");
    
    
    int N = 36;
    double alpha = 0;
    double dt = .18;// velocity verlet time step dt
    
    
    vector<double> x = vector<double>(N+1);
    vector<double> dxdt = x, d2xdt2 = x, a = x;
    
    double E_0=initial_condition(x, dxdt, d2xdt2, alpha);// starting potential energy at time 0
    
    
    double t=0;
    
    while (t<=150) {
       
            
            {
                double Kin_energy = Energy_kin(x, dxdt, alpha);
                double Pot_energy = Energy_pot(x, dxdt, alpha);
                double E_err=(E_0-Kin_energy-Pot_energy)/E_0;// error in total energy found numerically
                
                fprintf(f,"%.10f\t%.10f\t%.10f\t%.10f\n",t,Kin_energy,Pot_energy,E_err);
            }
            
            
            
        
        
        velocity_Verlet(dt, x, dxdt, d2xdt2, alpha);//updating velocity at every time step
        t+=dt;
    }
    
    return 0;
}












