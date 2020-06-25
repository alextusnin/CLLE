#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <time.h>
#include <random>
#include <cstdio>
#include <fstream>

#include "./../NR/NR_C301/code/nr3.h"
#include "./../NR/NR_C301/code/stepper.h"
#include "./../NR/NR_C301/code/stepperdopr853.h"
#include "./../NR/NR_C301/code/odeint.h"


#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
struct rhs_lle_pseud_spectral{
    Int Nphi;
    Doub det, f, d2, dphi;
    double* Dint;
    Complex i=1i;
    fftw_plan plan_direct_2_spectrum;
    fftw_plan plan_spectrum_2_direct;
    fftw_plan plan_disp_spectrum_2_direct;
    fftw_complex *buf_direct;
    fftw_complex *buf_spectrum;
    fftw_complex *buf_disp_direct;
    fftw_complex *buf_disp_spectrum;

    rhs_lle_pseud_spectral(Int Nphii, const double* Dinti, Doub deti, Doub fi, Doub d2i, Doub dphii)
    {
        Nphi = Nphii;
        buf_direct = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nphi);
        buf_spectrum = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nphi);
        buf_disp_direct = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nphi);
        buf_disp_spectrum = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nphi);
        plan_direct_2_spectrum = fftw_plan_dft_1d(Nphi, buf_direct,buf_spectrum, FFTW_BACKWARD, FFTW_ESTIMATE);
        plan_spectrum_2_direct = fftw_plan_dft_1d(Nphi, buf_spectrum,buf_direct, FFTW_FORWARD, FFTW_ESTIMATE);
        plan_disp_spectrum_2_direct = fftw_plan_dft_1d(Nphi, buf_disp_spectrum,buf_disp_direct, FFTW_FORWARD, FFTW_ESTIMATE);
        det = deti;
        f = fi;
        d2 = d2i;
        dphi = dphii;
        Dint = new (std::nothrow) double[Nphi];
        for (int i_phi = 0; i_phi<Nphi; i_phi++){
            Dint[i_phi] = Dinti[i_phi];
        }
    }
    ~rhs_lle_pseud_spectral()
    {
        delete [] Dint;
        free (buf_direct);
        free (buf_spectrum);
        free (buf_disp_direct);
        free (buf_disp_spectrum);
        fftw_destroy_plan(plan_direct_2_spectrum);
        fftw_destroy_plan(plan_spectrum_2_direct);
        fftw_destroy_plan(plan_disp_spectrum_2_direct);
    }
    void operator() (const Doub x, VecDoub &y, VecDoub &dydx) {
        for (int i_phi = 0; i_phi<Nphi; i_phi++){
            buf_direct[i_phi][0] = y[i_phi];
            buf_direct[i_phi][1] = y[i_phi+Nphi];
        }

        fftw_execute(plan_direct_2_spectrum);
        
        for (int i_phi = 0; i_phi<Nphi; i_phi++){
            buf_disp_spectrum[i_phi][0] = -Dint[i_phi]*buf_spectrum[i_phi][0];
            buf_disp_spectrum[i_phi][1] = -Dint[i_phi]*buf_spectrum[i_phi][1];
        }
        fftw_execute(plan_disp_spectrum_2_direct);
        
        
        for (int i_phi = 0; i_phi<Nphi; i_phi++){

            dydx[i_phi] = -y[i_phi] + det*y[i_phi+Nphi]  - buf_disp_direct[i_phi][1]/Nphi - (y[i_phi]*y[i_phi]+y[i_phi+Nphi]*y[i_phi+Nphi])*y[i_phi+Nphi] + f;
            dydx[i_phi+Nphi] = -y[i_phi+Nphi] - det*y[i_phi]  + buf_disp_direct[i_phi][0]/Nphi + (y[i_phi]*y[i_phi]+y[i_phi+Nphi]*y[i_phi+Nphi])*y[i_phi];

        }
    }


};

struct rhs_lle{
    Int Nphi;
    Doub det, f, d2, dphi;
    double* Dint;
    double* Disp;
    Complex i=1i;

    rhs_lle(Int Nphii, const double* Dinti, Doub deti, Doub fi, Doub d2i, Doub dphii)
    {
        Nphi = Nphii;
        det = deti;
        f = fi;
        d2 = d2i;
        dphi = dphii;
        Dint = new (std::nothrow) double[Nphi];
        Disp = new (std::nothrow) double[2*Nphi];
        for (int i_phi = 0; i_phi<Nphi; i_phi++){
            Dint[i_phi] = Dinti[i_phi];
        }
    }
    ~rhs_lle()
    {
        delete [] Dint;
        delete [] Disp;
    }
    void operator() (const Doub x, VecDoub &y, VecDoub &dydx) {
        
        Disp[0] = d2*(y[1] - 2*y[0]+ y[Nphi-1])/dphi/dphi;
        Disp[Nphi-1] = d2*(y[0] - 2*y[Nphi-1]+ y[Nphi-2])/dphi/dphi;
        
        Disp[Nphi] = d2*(y[Nphi+1] - 2*y[Nphi]+ y[2*Nphi-1])/dphi/dphi;
        Disp[2*Nphi-1] = d2*(y[Nphi] - 2*y[2*Nphi-1]+ y[2*Nphi-2])/dphi/dphi;


        for (int i_phi = 1; i_phi<Nphi-1; i_phi++){
            Disp[i_phi] = d2*(y[i_phi+1] - 2*y[i_phi]+ y[i_phi-1])/dphi/dphi;
            Disp[i_phi+Nphi] = d2*(y[i_phi+Nphi+1] - 2*y[i_phi+Nphi]+ y[i_phi+Nphi-1])/dphi/dphi;
        }
        
        for (int i_phi = 0; i_phi<Nphi; i_phi++){

            dydx[i_phi] = -y[i_phi] + det*y[i_phi+Nphi]  - Disp[i_phi+Nphi] - (y[i_phi]*y[i_phi]+y[i_phi+Nphi]*y[i_phi+Nphi])*y[i_phi+Nphi] + f;
            dydx[i_phi+Nphi] = -y[i_phi+Nphi] - det*y[i_phi]  + Disp[i_phi]+ (y[i_phi]*y[i_phi]+y[i_phi+Nphi]*y[i_phi+Nphi])*y[i_phi];

        }
    }


};

void printProgress (double percentage)
{
    int val = (int) (percentage*100 );
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}

std::complex<double>* InitialValue( const double f, const double detuning, const int Nphi)
{
    std::cout<<"Initializing linear resonance\n";
    std::complex<double> i=1i;
    std::complex<double>* res = new (std::nothrow) std::complex<double>[Nphi];
    for (int i_phi=0; i_phi<Nphi; i_phi++){
      res[i_phi] = f/(1.+i*detuning);
    }
    return res;
}

std::complex<double>* WhiteNoise(const double amp, const int Nphi)
{
    
    std::complex<double>* noise_spectrum = new (std::nothrow) std::complex<double>[Nphi];//contains white noise in spectal domain
    std::complex<double>* res = new (std::nothrow) std::complex<double>[Nphi];//contains white noise in spectal domain
    fftw_complex noise_direct[Nphi];
    fftw_plan p;
    
    p = fftw_plan_dft_1d(Nphi, reinterpret_cast<fftw_complex*>(noise_spectrum), noise_direct, FFTW_BACKWARD, FFTW_ESTIMATE);
    double phase;
    double noise_amp;
    const std::complex<double> i(0, 1);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    for (int j=0; j<Nphi; j++){
       phase = distribution(generator) *2*M_PI-M_PI;
       noise_amp  = distribution(generator) *amp;
       noise_spectrum[j] = noise_amp *std::exp(i*phase)/sqrt(Nphi);
    }


    fftw_execute(p);
    for (int j=0; j<Nphi; j++){
        res[j].real(noise_direct[j][0]);
        res[j].imag(noise_direct[j][1]);
    }
    fftw_destroy_plan(p);
    delete [] noise_spectrum;
    return res;
}

std::complex<double>** PropagateDOPRI(const double f,  const double *detuning, const double J, const double *phi, const double* Dint, const int Ndet, const int Nt, const double dt,  int Nphi, double noise_amp){
    std::complex<double> **res = new (std::nothrow) std::complex<double>*[Ndet];
    VecDoub res_buf(2*Nphi);
    const Doub atol = 1e-9, rtol=atol, dtmin=0.0, t0=0.,t1=dt*(Nt-1);

    for (int i=0; i<Ndet; i++){
        res[i] = new (std::nothrow) std::complex<double>[Nphi];
    }
    std::complex<double>* noise = new (std::nothrow) std::complex<double>[Nphi];
    VecDoub noise_buf(2*Nphi);
    res[0] = InitialValue(f, detuning[0], Nphi);
    noise=WhiteNoise(noise_amp,Nphi);
    for (int i_phi=0; i_phi<Nphi; i_phi++){
        res[0][i_phi]+=noise[i_phi];
        res_buf[i_phi] = res[0][i_phi].real(); 
        res_buf[i_phi+Nphi] = res[0][i_phi].imag();

    }
    Output out; 
    //rhs_lle lle(Nphi, Dint, detuning[0],f,Dint[1],std::abs(phi[1]-phi[0]));
    rhs_lle_pseud_spectral lle(Nphi, Dint, detuning[0],f,Dint[1],std::abs(phi[1]-phi[0]));
    for (int i_det=0; i_det<Ndet; i_det++){
        lle.det = detuning[i_det];
        noise=WhiteNoise(noise_amp,Nphi);
        //Odeint<StepperDopr853<rhs_lle> > ode(res_buf,t0,t1,atol,rtol,dt,dtmin,out,lle);
        Odeint<StepperDopr853<rhs_lle_pseud_spectral> > ode(res_buf,t0,t1,atol,rtol,dt,dtmin,out,lle);
        ode.integrate();
        for (int i_phi=0; i_phi<Nphi; i_phi++){
            res[i_det][i_phi].real(res_buf[i_phi]);
            res[i_det][i_phi].imag(res_buf[i_phi+Nphi]);
            res[i_det][i_phi]+=noise[i_phi];
            res_buf[i_phi] = res[i_det][i_phi].real(); 
            res_buf[i_phi+Nphi] = res[i_det][i_phi].imag();

        }
        printProgress((i_det+1.)/Ndet);
        
    }
    delete [] noise;
    return res;
}



std::complex<double>** PropagateSS(const double f,  const double *detuning, const double J, const double *phi, const double* Dint, const int Ndet, const int Nt, const double dt, const int Nphi, double noise_amp)
{
    
    std::complex<double> i = 1i;
    int check;
    
    std::complex<double> **res = new (std::nothrow) std::complex<double>*[Ndet];
    for (int i=0; i<Ndet; i++){
        res[i] = new (std::nothrow) std::complex<double>[Nphi];
    }
    std::complex<double>* noise = new (std::nothrow) std::complex<double>[Nphi];
    
    res[0] = InitialValue(f, detuning[0], Nphi);
    noise=WhiteNoise(noise_amp,Nphi);
    
    fftw_plan plan_direct_2_spectrum;
    fftw_plan plan_spectrum_2_direct;
    
    std::complex<double> buf;
    fftw_complex buf_direct[Nphi], buf_spectrum[Nphi];
    plan_direct_2_spectrum = fftw_plan_dft_1d(Nphi, buf_direct,buf_spectrum, FFTW_BACKWARD, FFTW_PATIENT);
    plan_spectrum_2_direct = fftw_plan_dft_1d(Nphi, buf_spectrum,buf_direct, FFTW_FORWARD, FFTW_PATIENT);
    
    for (int i_phi=0; i_phi<Nphi; i_phi++){
          buf_direct[i_phi][0] = res[0][i_phi].real() + noise[i_phi].real();
          buf_direct[i_phi][1] = res[0][i_phi].imag() + noise[i_phi].imag();
    }
    fftw_execute(plan_direct_2_spectrum);
    std::cout<<"Split step is running\n";

    for (int i_det=0; i_det<Ndet; i_det++){
        noise=WhiteNoise(noise_amp,Nphi);
        for (int i_t=0; i_t<Nt; i_t++){
            for (int i_phi=0; i_phi<Nphi; i_phi++){
                buf.real( buf_direct[i_phi][0] );
                buf.imag( buf_direct[i_phi][1]);
                buf+=(noise[i_phi]);
                buf*= std::exp(dt *(i*buf*std::conj(buf)  +i*J*(std::cos(phi[i_phi])+0.*std::sin(2*phi[i_phi]))  ) );
                buf_direct[i_phi][0] = buf.real();
                buf_direct[i_phi][1] = buf.imag();
            }
            fftw_execute(plan_direct_2_spectrum);//First step terminated
            for (int i_phi=0; i_phi<Nphi; i_phi++){
                buf.real( buf_spectrum[i_phi][0]);
                buf.imag( buf_spectrum[i_phi][1]);
                buf = std::exp(dt * (-1. - i*detuning[i_det] - i*Dint[i_phi]))*buf + f*Nphi/(-1. - i*detuning[i_det] - i*Dint[i_phi]) *(std::exp(dt * (-1. - i*detuning[i_det] - i*Dint[i_phi] ) ) - 1.) *((i_phi==0)? 1.0 : 0.0) ; 
                buf_spectrum[i_phi][0] = buf.real()/Nphi;
                buf_spectrum[i_phi][1] = buf.imag()/Nphi;
            }
            fftw_execute(plan_spectrum_2_direct);
            //Second step terminated
        }
        for (int i_phi=0; i_phi<Nphi; i_phi++){
            res[i_det][i_phi].real(buf_direct[i_phi][0]);
            res[i_det][i_phi].imag(buf_direct[i_phi][1]);
        }
        //std::cout<<(i_det+1.)/Ndet*100.<<"% is done\n";
        printProgress((i_det+1.)/Ndet);
    }
    fftw_destroy_plan(plan_direct_2_spectrum);
    fftw_destroy_plan(plan_spectrum_2_direct);
    std::cout<<"Split step is finished\n";
    return res;
}

void SaveData( std::complex<double> **A, const double *detuning, const double *phi, const int Ndet, const int Nphi)
{

    std::ofstream outFile;
    outFile.open("Field.bin", std::ios::binary);
    for (int i =0; i<Ndet; i++){
        for (int j=0; j<Nphi; j++){
            outFile.write(reinterpret_cast<const char*>(&A[i][j]),sizeof(std::complex<double>));
        }
    }
    outFile.close();
    outFile.open("Detuning.bin", std::ios::binary);
    for (int i =0; i<Ndet; i++){
        outFile.write(reinterpret_cast<const char*>(&detuning[i]),sizeof(double));
    }
    outFile.close();
    outFile.open("Phi.bin", std::ios::binary);
    for (int i =0; i<Nphi; i++){
        outFile.write(reinterpret_cast<const char*>(&phi[i]),sizeof(double));
    }
    outFile.close();
}
int main(int argc, char* argv[])
{
    double det_min = -8.0;
    double det_max = 20.0;
    const int Ndet = 400;
    double delta_det = (det_max-det_min)/(Ndet-1);
    double f = sqrt(14.1);
    double *detuning = new (std::nothrow) double[Ndet];
    for (int i = 0; i<Ndet; i++) detuning[i] = det_min+i*delta_det;
    double Tmax = 1.;
    double dt=0.5e-3;
    int Nt = int(Tmax/dt)+1;
    const int Nphi = pow(2,9);
    double *phi = new (std::nothrow) double[Nphi];
    double dphi = 2*M_PI/(Nphi-1);
    double noise_amp = 1e-9;
    
    double *fftw_freq = new (std::nothrow) double[Nphi];
    double *Dint = new (std::nothrow) double[Nphi];
    double D2 = 2*M_PI*0.48; //MHz
    double kappa = 2*M_PI*50; //MHz
    double d2 = D2/kappa;
    double J = -9.*0;
    if (Nphi%2 == 0){
        for (int i = 0; i<Nphi/2; i++){
            fftw_freq[i] = (i);
            Dint[i] = fftw_freq[i]*fftw_freq[i]*d2;
        }
        for (int i = Nphi/2; i<Nphi; i++){
            fftw_freq[i] = (i-Nphi);
            Dint[i] = fftw_freq[i]*fftw_freq[i]*d2;
        }
        
    }
    else{
        for (int i = 0; i<(Nphi+1)/2; i++){
            fftw_freq[i] = (i);
            Dint[i] = fftw_freq[i]*fftw_freq[i]*d2;
        }
        for (int i = (Nphi+1)/2; i<Nphi; i++){
            fftw_freq[i] = (i-Nphi);
            Dint[i] = fftw_freq[i]*fftw_freq[i]*d2;
        }
        
    }
    for (int i = 0; i<Nphi; i++){
        phi[i] = i*dphi;
    }
    std::complex<double> **A;
    //A= PropagateSS(f, detuning,J, phi,   Dint ,Ndet, Nt, dt, Nphi,noise_amp);
    A= PropagateDOPRI(f, detuning,J, phi,   Dint ,Ndet, Nt, dt, Nphi,noise_amp);
    
    SaveData(A, detuning, phi, Ndet, Nphi);
    for (int i=0; i<Ndet; ++i){
        delete[] A[i];
    }
    delete [] A;
    delete [] fftw_freq;
    delete [] Dint;
    return 0;
}
