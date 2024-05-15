#include "iFourier.h"

iFourier::iFourier(std::vector<Complex>& points)
	: Points{ points }
{
}

Complex iFourier::setPoints(std::vector<Complex>& points)
{
	Complex average = 0;
	for (unsigned int i = 0; i < points.size(); i++)
		average += points[i];
	average = average / (int)points.size();
	for (unsigned int i = 0; i < points.size(); i++)
		points[i] += -average;
	Points = points;
	return -average;
}

void iFourier::setPoints(int number)
{
	Points.clear();
	for (int i = 0; i < number; i++) {
		Points.push_back(eval(i * 2 * Pi / number));
	}
	Points.push_back(Points[0]);
}

void iFourier::setOrder(std::vector<int> order)
{
	Order = order;
}

void iFourier::setTime(float t)
{
	time = t;
	while (time > 1)
		time--;
}

std::vector<Complex>& iFourier::getPoints()
{
	return Points;
}

std::vector<Complex>& iFourier::getCoef()
{
	return Coef;
}

std::vector<int>& iFourier::getOrder()
{
	return Order;
}

std::vector<Complex>& iFourier::getPlot()
{
	return Plot;
}

float iFourier::getTime()
{
	return time;
}

float iFourier::getNorm2()
{
	return Norm2;
}

float iFourier::getdNorm2()
{
	return dNorm2;
}

float iFourier::getddNorm2()
{
	return ddNorm2;
}

Complex iFourier::getCircleCenter(int i)
{
	return eval(2 * Pi * time, i);
}

float iFourier::getCircleRadius(int i)
{
	return Coef[i].abs();
}

float iFourier::norm2(std::vector<Complex>& points)
{
	float x = 0;
	int m = points.size();
	for (int i = 0; i < m; i++)
		x += points[i].a * points[i].a + points[i].b * points[i].b;
	return x / m;
}

float iFourier::dnorm2(std::vector<Complex>& points)
{
	std::vector<Complex> dpoints;
	int m = points.size();
	for (int i = 1; i < m; i++)
		dpoints.push_back((points[i] - points[i - 1]) * m / (2 * Pi));
	return norm2(dpoints);
}

float iFourier::ddnorm2(std::vector<Complex>& points)
{
	std::vector<Complex> dpoints;
	int m = points.size();
	for (int i = 1; i < m; i++)
		dpoints.push_back((points[i] - points[i - 1]) * m / (2 * Pi));
	std::vector<Complex> ddpoints;
	m = dpoints.size();
	for (int i = 1; i < m; i++)
		ddpoints.push_back((dpoints[i] - dpoints[i - 1]) * m / (2 * Pi));
	return norm2(ddpoints);
}

float iFourier::norm2()
{
	float x = 0;
	int m = Points.size();
	for (int i = 0; i < m; i++)
		x += Points[i].a * Points[i].a + Points[i].b * Points[i].b;
	Norm2 = x / m;
	return Norm2;
}

float iFourier::dnorm2()
{
	std::vector<Complex> dpoints;
	int m = Points.size();
	for (int i = 1; i < m; i++)
		dpoints.push_back((Points[i] - Points[i - 1]) * m / (2 * Pi));
	dNorm2 = norm2(dpoints);
	return dNorm2;
}

float iFourier::ddnorm2()
{
	std::vector<Complex> dpoints;
	int m = Points.size();
	for (int i = 1; i < m; i++)
		dpoints.push_back((Points[i] - Points[i - 1]) * m / (2 * Pi));
	std::vector<Complex> ddpoints;
	m = dpoints.size();
	for (int i = 1; i < m; i++)
		ddpoints.push_back((dpoints[i] - dpoints[i - 1]) * m / (2 * Pi));
	ddNorm2 = norm2(ddpoints);
	return ddNorm2;
}

float iFourier::discreteError()
{
	int size = Points.size();
	std::vector<Complex> error;
	for (int i = 0; i < size; i++)
		error.push_back(Points[i] - eval(2 * i * Pi / size));
	return sqrt(norm2(error) / Norm2);
}

void iFourier::computeDiffPoints()
{
	diffPoints.clear();
	int m = Points.size();

	for (int p = 0; p < m; p++)
		diffPoints.push_back(Points[p] - eval(2 * Pi * p / m));
}

void iFourier::preCompute(int apriori)
{
	Lf.clear();
	for (int i = 1; i <= apriori; i++)
	{
		Lf.push_back(i);
		Lf.push_back(-i);
	}
	OrderLf();
	generateCoef();
	computeDiffPoints();
}

void iFourier::computeLf(float error)
{
	if (error)
		Error = error;

	float N = norm2(diffPoints);
	float dN = dnorm2(diffPoints);
	float ddN = ddnorm2(diffPoints);
	float moderror = error * sqrtf(Norm2 / N);
	float delta = 1.f / moderror * sqrt(ddN * N - dN * dN);
	float l0 = (dN - delta) / N;
	float l1 = (dN + delta) / N;
	int l = 1;
	if(Lf.size())
		l = -Lf[Lf.size() - 1] + 1;
	while (l * l < l1) {
		if (l * l > l0)
		{
			Lf.push_back(l);
			Lf.push_back(-l);
		}
		l++;
	}
}

void iFourier::OrderLf()
{
	Order = Lf;
}

void iFourier::fullComputation(float error, int apriori)
{
	norm2();
	dnorm2();
	ddnorm2();
	preCompute(apriori);
	computeLf(error);
	OrderLf();
	generateCoef();
}

void iFourier::generateCoef()
{
	Coef.clear();
	int m = Points.size();
	for (unsigned int i = 0; i < Order.size(); i++) {
		Coef.push_back(Complex(0));
		for (unsigned int j = 0; j < m; j++)
			Coef[i] += Points[j] * Complex::exp(-I * float(2 * Pi * j * Order[i] / m));
		Coef[i] = Coef[i] / float(m);
	}
}

void iFourier::generatePlot(int s)
{
	Plot.clear();
	for (int i = 0; i < s; i++) {
		Plot.push_back(eval(2 * Pi * i / s));
	}
	Plot.push_back(Plot[0]);
}

Complex iFourier::eval(float t, unsigned int d)
{
	if (!d)
		d = Order.size();
	Complex a(0);
	for (unsigned int i = 0; i < d; i++)
		a += Coef[i] * Complex::exp(I * Order[i] * t);
	return a;
}

void iFourier::renderPlot(Renderer& renderer)
{
	renderer.RenderPlot(Plot);
}

void iFourier::renderPartialPlot(Renderer& renderer)
{
	int i = int(Plot.size() * time);
	std::vector<Complex> plot(Plot.begin(), Plot.begin() + i);
	plot.push_back(eval(time * 2 * Pi));
	renderer.RenderPlot(plot);
}

void iFourier::renderCircles(Renderer& renderer, bool opaque)
{
	renderer.RenderPoint(0);
	std::vector<Complex> Circle;
	Complex center, nextcenter(0);
	for (unsigned int i = 0; i < Order.size(); i++) {
		int points = int(Coef[i].abs() * renderer.getScale() * 2 * Pi);
		if (!points)
			continue;
		if (points > CircleSmoothness)
			points = CircleSmoothness;
		center = nextcenter;
		nextcenter = eval(2 * Pi * time, i + 1);
		Circle.clear();
		
		for (int j = 0; j <= points; j++)
			Circle.push_back(Coef[i] * Complex::exp(I * 2 * Pi * j / points) + center);
		if (!opaque){
			renderer.RenderPlot(Circle, CircleColor * Color(255,255,255,100));
			renderer.RenderLine(center, nextcenter, LineColor * Color(255, 255, 255, 100));
			renderer.RenderPoint(nextcenter, Color(100,100,100,100));
		}
		else {
			renderer.RenderPlot(Circle, CircleColor);
			renderer.RenderLine(center, nextcenter, LineColor);
			renderer.RenderPoint(nextcenter);
		}
		
	}
}

void iFourier::renderPoints(Renderer& renderer)
{
	renderer.RenderPoints(Points);
}
