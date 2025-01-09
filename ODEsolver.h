#pragma once

#include "VectorOperations.h"

#include <array>
#include <cmath>
#include <functional>


// ====================================================================
// auto f = [](double t, const array<double> &x) -> array<double>
// ====================================================================

// ===============================================================
// Classe astratta, restituisce la derivata nel punto x al tempo t
// ===============================================================
class FunzioneVettorialeBase {

public:
  virtual ~FunzioneVettorialeBase() {}

  virtual vector<double> Eval(double t, const vector<double> &x) const = 0;

  vector<double> operator()(double t, const vector<double> &x) const {
    return Eval(t, x);
  } //per usare Eval posso fare: f(x) ( = f.Eval(x))
};

// ===============================================================
// Oscillatore Armonico semplice, smorzato o forzato
// ===============================================================
class OscillatoreArmonico : public FunzioneVettorialeBase {

public:
  OscillatoreArmonico() : m_omega0{1.}, m_gamma{0.}, m_omega{0.} {}

  OscillatoreArmonico(double omega0)
    : m_omega0{omega0}, m_gamma{0.}, m_omega{0.} {}

  OscillatoreArmonico(double omega0, double gamma)
    : m_omega0{1.}, m_gamma{gamma}, m_omega{0.} {}

  OscillatoreArmonico(double omega0, double gamma, double omega)
    : m_omega0{omega0}, m_gamma{gamma}, m_omega{omega} {}

  ~OscillatoreArmonico() override {}

  void setOmega(double omega) { m_omega = omega; }
  double getOmega() const { return m_omega; }

  vector<double> Eval(double t, const vector<double> &x) const override {
    return vector<double>{ x[1], -pow(m_omega0,2) * x[0] - m_gamma * x[1] + sin(m_omega * t)};
  }

private:
  double m_omega0, m_gamma, m_omega;
};

// ===============================================================
// Pendolo
// ===============================================================
class Pendolo : public FunzioneVettorialeBase {

public:
  Pendolo() : m_l{1.} {}
  Pendolo(double length) : m_l{length} {}

  ~Pendolo() override {}

  vector<double> Eval(double t,const vector<double> &x) const override {
    return vector<double>{x[1], -g * sin(x[0]) / m_l};
  }

private:
  double m_l;
  const double g{9.8};
};

// ===============================================================
// Gravitazione
// ===============================================================
class Gravitazione : public FunzioneVettorialeBase {
public:
  Gravitazione(double M) : m_M{M} {}
  ~Gravitazione() override {}

  double getM() const { return m_M; }
  void setM(double M) { m_M = M; }
  double const G{6.6742e-11};

  virtual vector<double>
  Eval(double t, const vector<double> &x) const override {
    double k = -G * m_M / pow(x[0] * x[0] + x[1] * x[1], 1.5);
    return vector<double>{x[2], x[3], k * x[0], k * x[1]};
  }

private:
  double m_M{};
};

// ===============================================================
// Elettromagnetismo
// ===============================================================
class Carica : public FunzioneVettorialeBase {
public:
  Carica(double M, double q) : m_M{M}, m_q{q} {}
  Carica(double M, double q, vector<double> E, vector<double> B)
    : m_M{M}, m_q{q}, m_E{std::move(E)}, m_B{std::move(B)} {}

  ~Carica() override {}

  double getM() const { return m_M; }
  double getC() const { return m_q; }

  double getEx() const { return m_E[0]; }
  double getEy() const { return m_E[1]; }
  double getEz() const { return m_E[2]; }

  double getBx() const { return m_B[0]; }
  double getBy() const { return m_B[1]; }
  double getBz() const { return m_B[2]; }

  void setM(double M) { m_M = M; }
  void setq(double q) { m_q = q; }

  void setEx(double Ex) { m_E[0] = Ex; }
  void setEy(double Ey) { m_E[1] = Ey; }
  void setEz(double Ez) { m_E[2] = Ez; }

  void setBx(double Bx) { m_B[0] = Bx; }
  void setBy(double By) { m_B[1] = By; }
  void setBz(double Bz) { m_B[2] = Bz; }

  void printB() const { std::cout << "B = " << m_B << "\n"; }
  void printE() const { std::cout << "E = " << m_E << "\n"; }

  vector<double>
  Eval(double t, const vector<double> &x) const override {

    vector<double> result;
    double k{m_q / m_M};
    result[0] = x[3]; // dx/dt
    result[1] = x[4]; // dy/dt
    result[2] = x[5]; // dz/dt

    result[3] = k * (x[4] * getBz() - x[5] * getBy() + getEx()); // yz-zy
    result[4] = k * (x[5] * getBx() - x[3] * getBz() + getEy()); // zx-xz
    result[5] = k * (x[3] * getBy() - x[4] * getBx() + getEz()); // xy-yx

    return result;
  }

private:
  double m_M{}, m_q{};
  vector<double> m_E{0., 0., 0.}, m_B{0., 0., 0.};
};

// ===============================================================
// Classe astratta per un integratore di equazioni differenziali
// ===============================================================
// template <std::uint16_t N> class EquazioneDifferenzialeBase {
// public:
// virtual vector<double> Passo(double t, const vector<double> &x, double h, 
//   const std::function< vector<double,N> (double,vector<double,N>)> &f) const = 0;

// };

// ===============================================================
// Integratore concreto, metodo di Eulero
// ===============================================================
template <std::uint16_t N> class Eulero{
public:
  Eulero() = default;
  Eulero(std::function< vector<double,N> (double,vector<double,N>)> function): f(function) {}

  void SetFunction(std::function< vector<double,N> (double,vector<double,N>)> function) { f = function; }

  vector<double> Passo(double t, const vector<double> &x, double h) const {
    
    return x + f(t, x) * h;
  }

   // =====================================================================================
  // Stima dell'errore dato un passo, vettore di dati iniziali e l'equazione differenziale
  // =====================================================================================
  vector<double,N> stima_err(double tmax,double h,const vector<double>& x){

    vector<double> x_N{x};
    vector<double> x_2N{x};
    unsigned int nstep{static_cast<unsigned int>(tmax/h+0.5)};
    double t{0.};
  
    for(uint i{};i<2*nstep;i++){
      if(i%2==0) x_N = Passo(t,x_N,h);
      x_2N = Passo(t,x_2N,h/2.);
      t += h/2.;
    }

    return fabs(x_N-x_2N)*2.;//errore di ordine 1
  }

  vector<double,4>  stima_h(const vector<double>& err_req,double tmax,double h,const vector<double>& x){
    
    vector<double,N> dt = stima_err(tmax,h,x);
    for(uint i=0;i< N;i++){
      dt[i] = (err_req[i]/dt[i]) * h;
    }
    //dereferenzio l'iteratore all'elemento minimo
    // return *std::min_element(dt.begin(),dt.end());
    return dt;
  }

  private:
  std::function<vector<double,N> (double,vector<double,N>)> f;
};

// ===============================================================
// Integratore concreto, metodo di Runge-Kutta ordine 4
// ===============================================================
template <std::uint16_t N>
class RungeKutta{
public:

  RungeKutta() = default;
  RungeKutta(std::function< vector<double,N> (double,vector<double,N>)> function): f{function} {}

  void SetFunction(std::function< vector<double,N> (double,vector<double,N>)> function) { f = function; }

  vector<double> Passo(double t, const vector<double> &x, double h) const{

    vector<double> k1{move(f(t, x))};
    vector<double> k2{move(f(t + h / 2, x + 0.5 * k1 * h))};
    vector<double> k3{move(f(t + h / 2, x + 0.5 * k2 * h))};
    vector<double> k4{move(f(t + h, x + k3 * h))};
    return x + (k1 + 2. * k2 + 2. * k3 + k4) * h / 6.;
  }

  // =====================================================================================
  // Stima dell'errore dato un passo, vettore di dati iniziali e l'equazione differenziale
  // =====================================================================================
  vector<double,N> stima_err(double tmax,double h,const vector<double>& x){

    vector<double> x_N{x};
    vector<double> x_2N{x};
    unsigned int nstep{static_cast<unsigned int>(tmax/h+0.5)};
    double t{0.};
  
    for(uint i{};i<2*nstep;i++){
      if(i%2==0) x_N = Passo(t,x_N,h);
      x_2N = Passo(t,x_2N,h/2.);
      t += h/2.;
    }

    return fabs(x_N-x_2N)*16./15.;
  }

  vector<double,4> stima_h(const vector<double>& err_req,double tmax,double h,const vector<double>& x){
    
    vector<double,N> dt = stima_err(tmax,h,x);
    for(uint i=0;i< N;i++){
      dt[i] = pow(err_req[i] / dt[i] ,0.25)* h;
    }
    //dereferenzio l'iteratore all'elemento minimo
    //return *std::min_element(dt.begin(),dt.end());
    return dt;
  }

  private:
  std::function<vector<double,N> (double,vector<double,N>)> f;
};  
