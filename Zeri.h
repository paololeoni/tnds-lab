#ifndef __bisezione__
#define __bisezione__

#include <float.h> // Serve per includere DBL_MAX

#include "Solutore.h"
#include <iostream>


using namespace std;

// funzione sign utile negli algoritmi di ricerca zeri. Potete aggiungerla come
// funzione che eredita da FunzionaBase 

int sign(double x){
  if(x<0)
    return -1;
	if(x==0){
		return 0;
	}
  else
    return 1;
}

class Bisezione : public Solutore {

public:

// costruttori 

	Bisezione () : Solutore() {;} ;
  Bisezione (double a, double b, double prec) : Solutore (a,b,prec) { ; } ;

  ~Bisezione () {;};

  double CercaZeri (double xmin, 
										double xmax, 
										const FunzioneBase& f,
										double prec = 0.001,
										unsigned int nmax = 100 ) override{ //se non specifico prec e nmax quando chiamo cercazeri questi sono i valori di default

	if(xmin>xmax){
		cout << "Errore, l'estremo superiore è minore dell'estremo inferiore, scambio i valori..." << endl;
		m_a = xmax;
		m_b = xmin;
	}

	else{
		m_a = xmin;
		m_b = xmax;
	}

	m_prec = prec;
	m_nmax = nmax;
	m_niterations = 0;

	if (sign(f.Eval(m_a)) * sign(f.Eval(m_b)) > 0) {
    cout << "non vale teorema degli zeri"
          << endl; // check per controllare che valga il Teo zeri
    exit(-1);
  }
		
	if (sign(f.Eval(m_a)) * sign(f.Eval(m_b)) == 0) {
    cout << "lo zero si trova in uno dei due estremi"
          << endl; // check per controllare che valga il Teo zeri
    exit(-1);
  }	
		
	while((abs(m_b - m_a)) > m_prec){

		double c = (m_b + m_a)*0.5;
		double fa = f.Eval(m_a);
		double fb = f.Eval(m_b);

		double fc = f.Eval(c);

		if(m_niterations > m_nmax){
			cout << "Numero massimo di iterzioni raggiunto.\n";
			break;
		}

		m_niterations++;
		
		if(sign(fa)*sign(fc)<0){
			m_b = c;
			fb = fc;
		}

		else if(sign(fb)*sign(fc)<0){
			m_a = c;
			fa = fc;	
		}

		else if(sign(fb)*sign(fc)==0 || sign(fa)*sign(fc)==0){
			return c; 	
		}

		else{
			cout << "Non ci sono zeri nell'intervallo proposto" << "\n";
		}

	//m_prec = abs(m_a-m_b);	 
	//double prec1 = abs(m_a-m_b);
	
	}

	m_prec = abs(m_a-m_b); //nota bene quanta differenza fa se lo metto qua o dentro il while (qua non rispetta la condizione del while)
	cout << "\n\t" << "iterazioni: " << m_niterations << endl;
	//m_prec = prec1;
	return (m_a+m_b)*0.5;

	
};
 

int sign1(double x){
  if(x<0)
    return -1;
  else
    return 1;
}

class Secante : public Solutore {

public:

// costruttori 

Secante () : Solutore() {;} ;
Secante (double xmin, double xmax, double prec) : Solutore (xmin, xmax,prec) { ; } ;

~Secante () {;};

double CercaZeri (double xmin, 
										double xmax, 
										const FunzioneBase& f,
										double prec = 0.001,
										unsigned int nmax = 100 ) override {


	if(xmin>xmax){
		cout << "Errore, l'estremo superiore è minore dell'estremo inferiore, scambio i valori..." << endl;
		m_a = xmax;
		m_b = xmin;
	}

	else{
		m_a = xmin;
		m_b = xmax;
	}

	m_prec = prec;
	m_nmax = nmax;
	m_niterations = 0;

	if (sign1(f.Eval(m_a)) * sign1(f.Eval(m_b)) > 0) {
    cout << "non vale teorema degli zeri"
          << endl; // check per controllare che valga il Teo zeri
    exit(-1);
  }

	double fa=f.Eval(m_a); ;
	double fb= f.Eval(m_b) ;
	double c = 0;

	
	do{
		
		c = (m_b - (fb *( m_b - m_a))/(fb - fa));

		fb = f.Eval(m_b);
		fa= f.Eval(m_a);
			
		if( sign1(f.Eval(c)) == sign1(fb)) m_b = c ;
		else m_a = c ;
		}while( abs((fb *( m_b - m_a))/(fb - fa)) > prec );

	return c;
	
}



};
};
#endif
