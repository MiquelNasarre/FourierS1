#pragma once

#include "Header.h"
#include "Button.h"
#include "DecimalBar.h"

class Settings {
public:
	static const unsigned int DiscreteErrorText = 3;
	static const unsigned int CoefficientsText  = 5;

	Button Play;
	Button Pause;
	Button IncreaseSpeed;
	Button DecreaseSpeed;
	Button DataWindow;
	Button Draw;
	Button Save;
	Button Load;
	Button Show;
	Button Hide;
	DecimalBar ErrorBar;
	DecimalBar NumberPoints;
	DecimalBar PreComputedFreq;

	std::vector<Text> texts;
public:
	Texture TexHoder;
	Sprite Holder;

	bool playing = false;
	float error = 0.05f;
	int points = 200;
	int apriori = 10;
	Font font;
	bool pointsVisible = false;

	Settings();

	void render(RenderWindow& window);
};