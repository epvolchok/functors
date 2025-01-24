#include <iostream>
#include <cmath>

//simple differentiation
template<typename F, typename T>
T inline fin_diff(F f, const T& x, const T& h)
{
	return (f(x+h)-f(x))/h;
}

inline double finn_diff(double f(double),  const double& x,const double& h)
{
	return (f(x+h)-f(x))/h;
}
//just some function
double sin_cos(double x)
{
	return sin(x) + cos(x);
}
//functor with the function and arbitrary parameter alpha
class psc_f
{
	private: 
	double alpha;
	
	public:
	psc_f() : alpha{1.0} {}
	psc_f(double alpha) : alpha(alpha) {}
	
	double operator () (double x) const
	{
		return alpha*(sin(x) + cos(x));
	}
	
};
//first derivative
//it can get for entry its own results
template <typename F, typename T>
class derivative
{
	private:
	const F& f;
	T h;
	public:
	derivative(const F& f, const T& h) : f(f), h(h) {}
	
	T operator () (const T& x) const
	{
		return (f(x+h) - f(x))/h;
	}
};
//second derivative
//without calculating the first one manually
template <typename F, typename T>
class second_derivative
{
	private:
	T h;
	derivative<F, T> fp;
	
	
	public:
	second_derivative(const F& f, const T& h) : h(h), fp(f, h) {}
	
	T operator () (const T& x) const
	{
		return (fp(x+h) - fp(x))/h;
	}
};
//n-th order derivative by recursion
template <typename F, typename T, unsigned N>
class nth_derivative
{
	private:
	using prev_derivative = nth_derivative<F, T, N-1>;
	T h;
	prev_derivative fp;
	
	public:
	nth_derivative(const F& f, const T& h) : h(h), fp(f, h) {}
	
	T operator () (const T& x) const
	{
		return (fp(x+h) - fp(x))/h;
	}
	 
};
/*template <typename F, typename T> //specification of the template to stop recursion
class nth_derivative<F, T, 1>
{
	private:
	const F& f;
	T h;
	public:
	nth_derivative(const F& f, const T& h) : f(f), h(h) {}
	
	T operator () (const T& x) const
	{
		return (f(x+h) - f(x))/h;
	}
};
*/
//or by inheritance
template <typename F, typename T> 
class nth_derivative<F, T, 1>
	: public derivative<F, T>
{
	using derivative<F, T>::derivative; //inheritates the basic class constructor
};

template <unsigned N, typename F, typename T>
nth_derivative<F, T, N> make_nth_derivative(const F& f, const T& h)
{
	return nth_derivative<F, T, N>(f, h);
}

int main()
{
	psc_f psc_o; //with function
	using d_psc_f = derivative<psc_f, double>; //first derivative
	d_psc_f d_psc_o(psc_o, 0.001);
	using dd_psc_f = derivative<d_psc_f, double>; // second derivative with substitution
	dd_psc_f dd_psc_o(d_psc_o, 0.001);
	using d2_psc_f = second_derivative<psc_f, double>; // second derivative
	d2_psc_f d2_psc_o(psc_o, 0.001);
	nth_derivative<psc_f, double, 3> d3_psc_o(psc_f(1.0), 0.001); // highest order derivatives
	auto d7_psc_o = make_nth_derivative<7>(psc_o, 0.0001);
	
	std::cout<<"ordinary func = "<<finn_diff(sin_cos, 1.0, 0.001)<<std::endl;
	std::cout<<"dif functor alpha=1 = "<<fin_diff(psc_o, 1.0, 0.001)<<std::endl; 
	std::cout<<"dif functor alpha=2 = "<<fin_diff(psc_f(2.0), 1.0, 0.001)<<std::endl;
	std::cout<<"1 derivative "<<d_psc_o(1.0)<<std::endl;
	std::cout<<"2 derivative manually "<<dd_psc_o(1.0)<<std::endl;
	std::cout<<"2nd derivative "<<d2_psc_o(1.0)<<std::endl;
	std::cout<<"3nd derivative "<<d3_psc_o(1.0)<<std::endl;
	std::cout<<"7th derivative "<<d7_psc_o(1.0)<<std::endl;
	
	//lambda expressions
	auto d7_psc_l = make_nth_derivative<7>([](double x){return sin(x) + cos(x); }, 0.0001);
	
	double alpha = 2.0, beta = 1.5;
	auto sc_l = [alpha, beta](double x)->double {return sin(alpha*x) + cos(beta*x); }; //parameters are copying but they're const
	
	auto sc_l2 = [&alpha, &beta](double x)->double {return sin(alpha*x) + cos(beta*x); };
	alpha = 2.5;
	beta = 1.0;
	double a = fin_diff(sc_l2, 1.0, 0.001); //can be changed here, even inside l-expression
	
	
	std::cout<<"lambda functor 1st derivative "<<fin_diff(sc_l, 1.0, 0.001)<<std::endl;
	std::cout<<"lambda functor 7th derivative "<<d7_psc_l(1.0)<<std::endl;
	
	
	
	
	return 0;
}
