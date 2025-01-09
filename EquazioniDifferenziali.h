#ifndef __EQ_DIFF__
#define __EQ_DIFF__

#include "VectorOperations.h"
#include <cmath>

using namespace std;

// ===============================================================
// Classe astratta, restituisce la derivata nel punto x al tempo t
// ===============================================================

class FunzioneVettorialeBase {

public:
  virtual vector<double> Eval(double t, const vector<double> & x) const = 0;
};

// ===============================================================
// Caso fisico concreto
// ===============================================================

class OscillatoreArmonico : public FunzioneVettorialeBase {

public:

  OscillatoreArmonico(double omega0) { m_omega0 = omega0; };
	//OscillatoreArmonico(double freq){ m_omega0 = 2*M_PI*freq};

	
  vector<double> Eval(double t, const vector<double> &x) const override { //valutazione derivata

		vector<double> v{0.,0.};
		v[0] = x[1]; //valutazione della velocità
		v[1] = x[0]*(-pow(m_omega0, 2)); //valutazione accelerazione

		return v;
	}


	
	void SetOmega(double omega){
		m_omega0 = omega;
	}

	void SetOmega2(double freq){
		m_omega0 = 2*M_PI*freq;	
	}

	double GetFreq() const{
		return m_omega0/(2*M_PI);
	}

	double GetOmega() const{
		return m_omega0;
	}

	

private:
  double m_omega0;
};


// ===============================================================
// Classe astratta per un integratore di equazioni differenziali
// ===============================================================

class EquazioneDifferenzialeBase {
public:
  virtual vector<double> Passo(double t, const vector<double> &x, double h, const FunzioneVettorialeBase &f) const = 0;
};

// ===============================================================
// Integratore concreto, metodo di Eulero
// ===============================================================

class Eulero : public EquazioneDifferenzialeBase {
public:

  vector<double> Passo(double t, const vector<double> &x, double h, const FunzioneVettorialeBase &f) const override {

		return x + (f.Eval(t,x))*h;

  }

  // =====================================================================================
  // Stima dell'errore dato un passo, vettore di dati iniziali e l'equazione differenziale
  // =====================================================================================
  vector<double> stima_err(double tmax,double h,const vector<double>& x, const FunzioneVettorialeBase &f){

    vector<double> x_N{x};
    vector<double> x_2N{x};
    unsigned int nstep{static_cast<unsigned int>(tmax/h+0.5)};
    double t{0.};
  
    for(uint i{};i<2*nstep;i++){
      if(i%2==0) x_N = Passo(t,x_N,h,f);
      x_2N = Passo(t,x_2N,h/2., f);
      t += h/2.;
    }

    return fabs(x_N-x_2N)*2.;//errore di ordine 1
  }

  vector<double>  stima_h(const vector<double>& err_req,double tmax,double h,const vector<double>& x, const FunzioneVettorialeBase &f){
    
    vector<double> dt = stima_err(tmax,h,x, f);
    for(uint i=0;i< dt.size();i++){
      dt[i] = (err_req[i]/dt[i]) * h;
    }
    //dereferenzio l'iteratore all'elemento minimo
    // return *std::min_element(dt.begin(),dt.end());
    return dt;
  }
};

class RK: public EquazioneDifferenzialeBase{
public:
  vector<double> Passo(double t, const vector<double> &x, double h, const FunzioneVettorialeBase &f) const override {

		vector<double> k1 = f.Eval(t,x);
		vector<double> k2 = f.Eval(t + h*0.5, x + (h*0.5*k1));
		vector<double> k3 = f.Eval(t + h*0.5, x + (h*0.5*k2));		
		vector<double> k4 = f.Eval(t + h, x + (h*k3));

		return x + (h/6.)*(k1 + k2*(2.) + (2.)*k3 + k4);
		
  }

  // =====================================================================================
  // Stima dell'errore dato un passo, vettore di dati iniziali e l'equazione differenziale
  // =====================================================================================
  vector<double> stima_err(double tmax,double h,const vector<double>& x, const FunzioneVettorialeBase &f){

    vector<double> x_N{x};
    vector<double> x_2N{x};
    unsigned int nstep{static_cast<unsigned int>(tmax/h+0.5)};
    double t{0.};
  
    for(uint i{};i<2*nstep;i++){
      if(i%2==0) x_N = Passo(t,x_N,h, f);
      x_2N = Passo(t,x_2N,h/2.);
      t += h/2.;
    }

    return fabs(x_N-x_2N)*16./15.;
  }

  vector<double> stima_h(const vector<double>& err_req,double tmax,double h,const vector<double>& x, const FunzioneVettorialeBase &f){
    
    vector<double> dt = stima_err(tmax,h,x, f);
    for(uint i=0;i< dt.size();i++){
      dt[i] = pow(err_req[i] / dt[i] ,0.25)* h;
    }
    //dereferenzio l'iteratore all'elemento minimo
    //return *std::min_element(dt.begin(),dt.end());
    return dt;
  }
};

#endif // __EQ_DIFF__
