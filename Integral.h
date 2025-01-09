#pragma once
#include "FunzioneBase.h"
#include "RandomGen.h"
#include "VectorOperations.h"
#include <cmath>
#include <iostream>

using namespace std;


class Integral {
  
protected:
  double m_a, m_b;
  int m_sign;
  unsigned int m_step;
  double m_prec;
  double m_err;

public:
  Integral() : m_a{}, m_b{}, m_sign{}, m_step(), m_prec{}, m_err{} {}

  double integrate(double a, double b, FunzioneBase &f,
                  double value) {
    checkInterval(a, b);
    m_err = 0.;
		if(value > 1) {
			m_step = static_cast<unsigned int> (value);
			m_prec= 0.;
			return calculate(m_step, f); // Fixed_step
		}
		else{
			m_step = 1u;
			m_prec = value;
			return calculate(f); // Fixed_Precision
		}
  }

  double getA() const { return m_a; }
  double getB() const { return m_b; }
  int getSign() const { return m_sign; }
  unsigned int getPartition() const { return m_step; }
  double getError() const { return m_err; }
  double getPrecision() const { return m_prec; }

  void setPrecision(double p) { m_prec = p; }

private:
  virtual double calculate(unsigned int nstep,
                           FunzioneBase &f) const = 0;

  virtual double calculate(FunzioneBase &f) = 0;

  void checkInterval(double a, double b) {
    m_a = std::min(a, b);
    m_b = std::max(a, b);
    m_sign = (a < b) ? 1 : -1;
  }
};

class MidPoint : public Integral {
public:
  MidPoint() : Integral() {}

private:
  // ================================================
  // Fixed Partition
  // ================================================

  double calculate(unsigned int nstep,
                   FunzioneBase &f) const override {
    double h{(m_b - m_a) / nstep};
    double sum{};

    for (unsigned int i{}; i < nstep; i++) {
      sum += f.Eval(m_a + (i + 0.5) * h);
    }

    return m_sign * sum * h;
  }
  // ================================================
  // Fixed Precision
  // ================================================

  virtual double calculate(FunzioneBase &f) override {
    double old_sum{};
    double sum{};

    do {
      m_step *= 2;
      old_sum = sum;                       // salvo l'integrale con N intervalli
      sum = calculate(m_step, f);          // calcolo con 2N intervalli
      m_err =  fabs(sum - old_sum) / 3.;   // Errore al secondo ordine
    } while (m_err > m_prec);

    return sum; // Formula al secondo ordine
  }
};

class Simpson : public Integral {
public:

  Simpson() : Integral() {}

private:
  // ================================================
  // Fixed Partition
  // ================================================

  double calculate(unsigned int nstep,
                   FunzioneBase &f) const override {
    unsigned int true_n = (nstep % 2 == 0) ? nstep : (nstep - 1);
    double h{(m_b - m_a) / true_n};
    double sum{(f.Eval(m_a) + f.Eval(m_b)) / 3.};

    for (unsigned int i{1}; i < true_n; i++) {
      sum += 2. * (1 + i % 2) * f.Eval(m_a + i * h) / 3.;
    }
    return m_sign * sum * h;
  }

  // ================================================
  // Fixed Precision
  // ================================================

  virtual double calculate(FunzioneBase &f) override {
    double old_sum{};
    double sum{};

    do {
      old_sum = sum;
      sum = calculate(m_step, f);
      m_err = ( fabs(sum - old_sum) / 15); // Errore al quarto ordine
      m_step *= 2;
    } while (m_err > m_prec);

    return sum; // Formula al quarto ordine
  }
};

class Trapezoide : public Integral {
public:
  Trapezoide() : Integral() {}

private:
  // ================================================
  // Fixed Partition
  // ================================================

  double calculate(unsigned int nstep,
                   FunzioneBase &f) const override {
    double sum{(f.Eval(m_a) + f.Eval(m_b)) / 2};
    double h{(m_b - m_a) / nstep}; // Passo h

    for (unsigned int i{1}; i < nstep; i++) {
      sum += f.Eval(m_a + i * h);
    }
    return m_sign * sum * h;
  }

  // ================================================
  // Fixed Precision
  // ================================================

  // Versione ottimizzata che tiene conto del calcolo precedente

  double calculate(FunzioneBase &f) override {
		double sum{(f.Eval(m_a) + f.Eval(m_b)) / 2}; //passo m_step = 1
    double h{(m_b - m_a)}; //passo m_step = 1
    double Int{}; // Valore dell'integrale
		double old_Int{}; // variabile di appoggio
		
    do {
      m_step *= 2; // raddoppio il numero di intervalli
			h /= 2.; // Dimezzo il passo h
      old_Int = Int;   // Salvo la somma precedente
      for (unsigned int i{1}; i < m_step; i += 2) { // Sommo sui dispari in modo che i/m_step siano coprimi
        sum += f.Eval(m_a + i * h);
      }
      Int = sum * h;
      m_err = ( fabs(Int - old_Int) / 3); // Aggiorno il nuovo errore
    } while (m_err > m_prec);

    return m_sign*Int; // Formula al secondo ordine
  }
};


class IntegralMC{
	
	public:
	
		IntegralMC(unsigned int seed){m_myrand = new RandomGen(seed);}
		~IntegralMC(){;}

		double IntegraMed(double a,double b, FunzioneBase & f, int punti){
			double sum = 0;
			
			for(int k = 0; k < punti; k++)
				sum += f.Eval(m_myrand->Unif(a,b));
			
			return static_cast<double>(b-a) * sum / static_cast<double> (punti);
		}

		double IntegraHOM(double a, double b, FunzioneBase & f, double max, int Ntot){
			
			int Nhit = 0;
			double x,y;
			for(int k = 0; k < Ntot;k++){
				if(f.Eval(m_myrand->Unif(a,b))>m_myrand->Unif(0,max))
					Nhit++;
			}
			return (b-a)*Nhit*max/static_cast<double> (Ntot);
		}

		double Ave_err(double a,double b,const FunzioneBase &f,unsigned int N){
			double media{},media2{},x{},appo{};
		  
			for(unsigned int i=0;i<N;i++){//calcolo media e media dei quadrati
		    x = f.Eval(m_myrand->Unif(a,b));
		    appo = static_cast<double>(i)/static_cast<double>(i+1);
		    media = appo*media + x/static_cast<double>(i+1);
		    media2 = appo*media2 + x/static_cast<double>(i+1)*x;
		  }
		  //la varianza è media dei quadrati meno quadrato della media (media2-media*media)
		  //L'errore della Media e' stimato con la deviazione standard della media delle valutazioni di f
			return (b-a)*sqrt((media2-media*media)/static_cast<double>(N-1));
		}

	double Average(const vector<double>& inf,const vector<double>& sup,
	const FunzioneBase& f,int punti){
	        
	  vector<double> a{m_myrand->Unif(inf[0],sup[0]),m_myrand->Unif(inf[1],sup[1])};
	  double accu{f.Eval2(a[0], a[1])};
	
	  for(int i{1};i<punti;i++){
	    a = vector<double>{m_myrand->Unif(inf[0],sup[0]),m_myrand->Unif(inf[1],sup[1])};
			//cout << a[0] << "\t"<< a[1] << endl;
	    accu+= f.Eval2(a[0], a[1]);
	  }
	  return (sup[0]-inf[0])*(sup[1]-inf[1])*accu/static_cast<double>(punti);
	}


	protected:

		RandomGen * m_myrand;
};
