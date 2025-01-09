#ifndef __FunzioneBase__
#define __FunzioneBase__

#include "Integral.h"


class FunzioneBase {

  public:

  virtual double Eval (double x) const = 0;
	//irtual  double Eval2 (	MidPoint myInt) const = 0;
  virtual ~FunzioneBase() {;};

};

#endif // __FunzioneBase__



#ifndef __f__
#define __f__


//============================================
// VARIE FIGLIE
// ===========================================
#include "FunzioneBase.h"
#include "Integral.h"
#include <cmath>

class fx : public FunzioneBase {

public:

	fx() {}
  //fx(double &n) : m_n{n} {};
  ~fx() {}

  double getd() { return d; }
	void setn(double &n) { m_n = n; };
	double getn() const { return m_n;}

  double Eval(double x) const override {
    return (1. / d) * cos(static_cast<double>((2. * M_PI / lambda)) * (sqrt(L * L + (m_n - x) * (m_n - x)) - sqrt(L * L + m_n * m_n)));
  }

	double Eval2(MidPoint myInt) const{
    return ;
		
	}

private:
  double lambda = 589E-9;
  double d = 10E-6;
  double L = 1.;
  double m_n;
};

#endif // __Sinx__
