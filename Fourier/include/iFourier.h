#pragma once

#include "Header.h"
#include "Complex.h"
#include "Renderer.h"

#define CircleSmoothness 400

#ifdef DarkTheme
	#define CircleColor Color(150,150,150)
	#define LineColor	Color::White
#else
	#define CircleColor Color(175,175,175)
	#define LineColor	Color::Black
#endif

class iFourier {
private:
	std::vector<Complex> Points;
	std::vector<Complex> diffPoints;
	std::vector<Complex> Coef;
	std::vector<int>	 Order;

	std::vector<Complex> Plot;

	float time;

	float Norm2;
	float dNorm2;
	float ddNorm2;
	float Error = 1;

	std::vector<int> Lf;

public:
	iFourier() {}
	iFourier(std::vector<Complex>& points);

	Complex setPoints(std::vector<Complex>& points);
	void setPoints(int number);
	void setOrder(std::vector<int> order);
	void setTime(float t);
	
	std::vector<Complex>&	getPoints();
	std::vector<Complex>&	getCoef();
	std::vector<int>&		getOrder();
	std::vector<Complex>&	getPlot();
	float					getTime();
	float					getNorm2();
	float					getdNorm2();
	float					getddNorm2();
	Complex					getCircleCenter(int i);
	float					getCircleRadius(int i);

	static float norm2(std::vector<Complex>& points);
	static float dnorm2(std::vector<Complex>& points);
	static float ddnorm2(std::vector<Complex>& points);

	float norm2();
	float dnorm2();
	float ddnorm2();
	float discreteError();

	void computeDiffPoints();
	void preCompute(int apriori);
	void computeLf(float error = 0.f);
	void OrderLf();

	void fullComputation(float error = 0.f, int preCompute = 0);

	void generateCoef();
	void generatePlot(int s);
	Complex eval(float t, unsigned int d = 0);

	void renderPlot(Renderer& renderer);
	void renderPartialPlot(Renderer& renderer);
	void renderCircles(Renderer& renderer, bool opaque = true);
	void renderPoints(Renderer& renderer);

};
