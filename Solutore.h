#ifndef __solutore__
#define __solutore__

#include "FunzioneBase.h"

class Solutore {

public:

	// Costruttori 

	Solutore() { m_a=0; m_b=0;} ;
	Solutore(double a, double b, double prec){m_a=a; m_b=b; m_prec=prec;};

	// altri ?

	// Distruttore virtuale

	virtual ~Solutore() {;}; 

  virtual double CercaZeri (double xmin, 
														double xmax, 
														const FunzioneBase& f,
														double prec = 0.001,
														unsigned int nmax = 100 ) = 0 ;

	// ======================================
	// 	metodi SET
	// ======================================

  void SetPrecisione (double epsilon) { m_prec=epsilon; }

  void SetNMaxiterations (unsigned int n) { m_nmax = n; }


	// ======================================
	// 	metodi GET
	// ======================================

  unsigned int GetNMaxiterations () {return m_nmax; }

  unsigned int GetNiterations () {return m_niterations; }

 	double GetPrecisione () { return m_prec; }


protected:

  double m_a, m_b;
  double m_prec;
  unsigned int m_nmax, m_niterations;
	
};

#endif // __solutore__
