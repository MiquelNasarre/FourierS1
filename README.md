# Fourier Series for $L^2(\mathbb{S}^1)$ Functions

This project is the first implementation of the formulas found on my final project. 
This one centers on functions from $\mathbb{S}^1$ to the complex plain.

It consists in an application which allows you to draw any closed curve. Then it takes those 
points and using the formulas shown in the *"Computation Aproximations"* section of the 
paper, it computes a number of Fourier coefficients according to the error you have chosen. 

It then displays the plot of the Fourier series in an interactive manner inspired in the 3blue1brown 
Youtube video [But what is a Fourier series?](https://www.youtube.com/watch?v=r6sGWTCMz2k&t=376s). 

https://github.com/MiquelNasarre/FourierS1/assets/124403865/fafab4c6-e3ca-4b7c-a5cb-ce5e8db3b79b

This video shows the opening example of the program.

## Requirements
- [Visual Studio](https://visualstudio.com) (For the program to be able to compile all the diffent C++ files correctly)
- [SFML 2.6.0](https://www.sfml-dev.org/download.php) (The files needed are already included in the repository)

If you are not planning to build and run the program yourself and you do not need the code, you can also find the 
program executable with its associated files in the following [link](). Please note that the program will only work on 
Windows 7+ operating systems.

## How to Use it
First, clone the repository into your computer and open the solution file with Visual Studio, then compile the only
project inside and the App should be up and running, you should compile it in 32bit (x86).
Note: if you don't operate in Windows you might have some compatibility problems with the SFML renderer.

## How it Works
Since this program is my first take on creating Windows desktop applications, it mostly relies on its dependencies for 
all the technical side. The program uses the SFML tools for all the window creation, interaction and rendering.

The programs structure is poorly organized, but after some effort is easy to follow its workings. The 
main file just calls the creation of a *Visual* class and runs the main Loop function, so all of the 
program revolves around such loop. Every frame it does an event check and renders everything accordingly, 
exiting the loop once the window is closed.

All of the widgets found on the menu bar are self designed, you can find them in the resources. Most 
of them have their own class and interactions. This is the only program with such features since the 
3D program uses the *imGui* functions for all the widgets.

For the Fourier series it has a class named *iFourier* which is used to handle most of the computations 
regarding the coefficients, as well as the plotting and rendering of the functions.

Also I would like to remark the use of the *Complex* struct, that later would be implemented in the 
[Math Library](https://github.com/MiquelNasarre/Math) used to do all calculations on the 3D program.

## Features and Math Involved
As mentioned before the main goal of this program is to showcase how the Fourier series work and implement 
the error functions found on the paper. For this reason I have implemented features that tackle all the 
math seen in the paper for this particular case. Lets go through its functionalities.

### Drawing
The drawing tool is mostly how you are supposed to design functions with the program. It is pretty straight 
forward and easy to use. You just press the *Draw* button and draw any shape you want in the complex plane 
using the mouse. Once you're done, press enter and the program will compute everything and display it as a 
function.

To see the points you have drawn just click on the *Show* button, then the points stored on the *Fourier* class will be displayed. 
You can also modify the number of points you want using the Points input on the top of the menu bar.

Now, how does it work? Once you press the *Draw* button, every time you are pressing the left button of the mouse 
and the mouse position changes it records a new location in an array. Once all the locations are recorded and you press 
*Enter*, the program considers every location is evenly spaced in time, and it distributes the function points accordingly. 

This assumption greatly simplifies the math in the paper, since all of the points are evenly spaced. Every time we compute 
$t_{i+1} -t_i$ it will just be $\frac{2\pi}{N}$ where $N$ is the total number of points. So the formulas we will be using 
will be simplified using this equality.

For computing the actual coefficients, after choosing which ones to compute $-$ which will be discussed in the following 
section $-$ and centering the figure in order that the first coefficient is $0$, it uses the modified formula from the paper

$$
f_n \approx \frac{1}{N} \sum_{i=0}^{N-1} \frac{z_i+z_{i+1}}{2} \overline{\psi_n} \left( 2\pi \frac{i + 1/2}{N} \right)
$$

where

$$
\psi_n (x) = e^{inx}
$$

and the $z_i$'s are the complex points drawn on the plain. 

Then to display it, it calculates the Fourier series for each and 
every coefficient in order to draw the next coefficient circle in the correct position, and for the last one it just displays 
the function up to a given value of $t$ to give the illusion that the circles are drawing the function. For reference de Fourier 
series looks as follows.

$$
f^L (x) = \sum_{\ell\in L} f_\ell \psi_\ell (x)
$$


### Error Computations
For knowing how many coefficients to compute we will be using the formulas from the paper. You can specify your desired error 
in the *Error* input just below the *Points* input. You can also specify how many coefficients will be computed prior to the estimation 
in order to make the estimation more precise. For reference the formula for the upper and lower bound of the coefficients to compute 
is

$$
L_f^\pm(\epsilon) = \frac{||df||^2 \pm \epsilon^{-1}
\sqrt{||\Delta f||^2||f||^2 - ||df||^4}}{||f||^2}
$$

Where $L_f^\pm$ determines the upper and lower bounds for the coefficients so that $\ell^2 \in \left\[L_f^- , L_f^+ \right]$.

In case we have already some coefficients computed $I\subset\Lambda$ we may use the formula for the new set of frequencies

$$
L_f(\epsilon,I) = I \cup L_{f^{\Lambda\setminus I}}
\left(\frac{\epsilon||f||}{||f^{\Lambda\setminus I}||}\right)
$$

To calculate all the norms we will be using the modified formulas from the paper

$$
||f||^2 \approx \frac{1}{4N}\sum_{i=0}^{N-1} |z_i + z_{i+1}|^2
$$

$$
||df||^2 \approx N \sum_{i=0}^{N-1} |z_{i+1} - z_i|^2
$$

$$
||\Delta f||^2 \approx 
2N^3\sum_{i=0}^{N-1} |z_{i+1} + z_{i-1} - 2z_i |^2
$$

Using all this as explained in the paper we obtain a subset of frequencies that will guarantee that the error is 
below our threshold. You can see the actual error obtained from the coefficients, as well as how many were calculated 
in the menu bar, where it says *Discrete error* and *Coefficients* respectively.

### Coefficients View
As an interactive and intuitive way to view the coefficients I created a small popup window in order to display them 
as a discrete graph. Just click the *Coefficients* button at the beggining of the menu bar and you will have a window 
like the following one.

![Coefficients](https://github.com/MiquelNasarre/FourierS1/assets/124403865/4f3d4d9e-b7b8-4ab3-8586-d8ba743114a8)

That way you can visually see which coefficients have a bigger importance on the series, looking like a Fourier transform.

### File Handling
The program also allows for saving and loading of figures, saved as sets of points, the files are stored on the 
[saveFiles](https://github.com/MiquelNasarre/FourierS1/tree/master/Fourier/saveFiles) folder. The format is pretty 
straight forward and you can use this feature in order to load files from other programs you might have or to save 
any of your drawings.

## Issues
Since this is not the main program of the final project, and it is my first attempt at graphics programming, it has 
a bunch of problems which will not be adressed in the coming future.

It does not allow for resizing or full screen. Which as a regular user of desktop applications can get really frustrating.

It has a strong dependence on external code, so it does not allow for great control of the features in your program.

It does not allow to save the coefficients it calculates, which would be interesting if you want to plot your Fourier Series somewhere else.

It does not allow for multiple plots at the same time or interpolation between them, which is a nice feature to implement in such programs.

Overall the program is not complete, but still it is able to do what it is supposed to, and it made the grounding for 
the $\mathbb{S}^2$ version of the program, which has much more care and effort put into it and I greatly reccomment you take a look at [it]().
