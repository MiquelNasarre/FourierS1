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
	float x = 0.f;
	int N = points.size();
	for (int i = 0; i < N; i++)
	{
		Complex z = points[i] + points[(i + 1) % N];
		x += z.a * z.a + z.b * z.b;
	}
	return 0.25f * x / N;
}

float iFourier::dnorm2(std::vector<Complex>& points)
{
	float x = 0.f;
	int N = points.size();
	for (int i = 0; i < N; i++)
	{
		Complex z = points[(i + 1) % N] - points[i];
		x += z.a * z.a + z.b * z.b;
	}
	return N * x;
}

float iFourier::ddnorm2(std::vector<Complex>& points)
{
	float x = 0.f;
	int N = points.size();
	for (int i = 0; i < N; i++)
	{
		Complex z = points[(i + 1) % N] + points[(i - 1 + N) % N] - 2 * points[i];
		x += z.a * z.a + z.b * z.b;
	}
	return 2 * N * N * N * x;
}

float iFourier::norm2()
{
	float x = 0.f;
	int N = Points.size();
	for (int i = 0; i < N; i++)
	{
		Complex z = Points[i] + Points[(i + 1) % N];
		x += z.a * z.a + z.b * z.b;
	}
	Norm2 = 0.25f * x / N;
	return Norm2;
}

float iFourier::dnorm2()
{
	float x = 0.f;
	int N = Points.size();
	for (int i = 0; i < N; i++)
	{
		Complex z = Points[(i + 1) % N] - Points[i];
		x += z.a * z.a + z.b * z.b;
	}
	dNorm2 = N * x;
	return dNorm2;
}

float iFourier::ddnorm2()
{
	float x = 0.f;
	int N = Points.size();
	for (int i = 0; i < N; i++)
	{
		Complex z = Points[(i + 1) % N] + Points[(i - 1 + N) % N] - 2 * Points[i];
		x += z.a * z.a + z.b * z.b;
	}
	ddNorm2 = 2 * N * N * N * x;
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
	float delta = 1.f / moderror * sqrt(abs(ddN * N - dN * dN));
	float l0 = (dN - delta) / N;
	float l1 = (dN + delta) / N;
	int l = 1;
	if(Lf.size())
		l = -Lf[Lf.size() - 1] + 1;
	while (4 * Pi * Pi * l * l < l1) {
		if (4 * Pi * Pi * l * l > l0)
		{
			Lf.push_back(l);
			Lf.push_back(-l);
		}
		l++;
	}
	printf("\n||f|| = %f\n||df|| = %f\n||Df|| = %f\nDesviacio = %f\nLOWER BOUND: %f\nUPPER BOUND: %f\n", N, dN, ddN, ddN * N - dN * dN, l0, l1);
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
	int N = Points.size();
	for (unsigned int i = 0; i < Order.size(); i++) {
		Coef.push_back(Complex(0));
		for (int j = 0; j < N; j++)
			Coef[i] += (Points[j] + Points[j + 1]) / 2.f * Complex::exp(-I * 2.f * Pi * (j + 0.5) * float(Order[i]) / float(N));
		Coef[i] = Coef[i] / float(N);
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
