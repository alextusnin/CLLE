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
std::complex<double>* InitialValue( const double f, const double detuning, const int Nphi)
{
    std::cout<<"Initializing linear resonance\n";
    std::complex<double> i = 1i;
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
//        noise=WhiteNoise(noise_amp,Nphi);
        for (int i_t=0; i_t<Nt; i_t++){
            for (int i_phi=0; i_phi<Nphi; i_phi++){
                buf.real( buf_direct[i_phi][0] );
                buf.imag( buf_direct[i_phi][1]);
                buf+=(noise[i_phi]);
                buf*= std::exp(dt *(i*buf*std::conj(buf) + f/buf +i*J*(std::cos(phi[i_phi])+0.45*std::sin(2*phi[i_phi]))  ) );
                buf_direct[i_phi][0] = buf.real();
                buf_direct[i_phi][1] = buf.imag();
            }
            fftw_execute(plan_direct_2_spectrum);//First step terminated
            for (int i_phi=0; i_phi<Nphi; i_phi++){
                buf.real( buf_spectrum[i_phi][0]);
                buf.imag( buf_spectrum[i_phi][1]);
                buf *= std::exp(dt * (-1. - i*detuning[i_det] - i*Dint[i_phi] )  ); 
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
    double det_min = -30.;
    double det_max = 30.;
    const int Ndet = 750;
    double delta_det = (det_max-det_min)/(Ndet-1);
    double f = sqrt(10.);
    double *detuning = new (std::nothrow) double[Ndet];
    for (int i = 0; i<Ndet; i++) detuning[i] = det_min+i*delta_det;
    double Tmax = 10.;
    double dt=1e-4;
    int Nt = int(Tmax/dt)+1;
    const int Nphi = pow(2,9);
    double *phi = new (std::nothrow) double[Nphi];
    double dphi = 2*M_PI/(Nphi-1);
    double noise_amp = 1e-12;
    
    double *fftw_freq = new (std::nothrow) double[Nphi];
    double *Dint = new (std::nothrow) double[Nphi];
    double D2 = 2*M_PI*0.48; //MHz
    double kappa = 2*M_PI*50; //MHz
    double d2 = D2/kappa;
    double J = -9.;
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
    A= PropagateSS(f, detuning,J, phi,   Dint ,Ndet, Nt, dt, Nphi,noise_amp);
    
    SaveData(A, detuning, phi, Ndet, Nphi);
    for (int i=0; i<Ndet; ++i){
        delete[] A[i];
    }
    delete [] A;
    delete [] fftw_freq;
    delete [] Dint;
    return 0;
}
