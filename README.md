# Fourier Series for $L^2(\mathbb{S}^1)$ Functions

This project is the first implementation of the formulas found on my final project. 
This one centers on functions from $\mathbb{S}^1$ to the complex plane.

It consisists in an application which allows you to draw any closed curve. Then it takes those 
points and using the formulas shown in the *"Computation Aproximations"* section of the 
paper, it computes a number of Fourier coefficients according to the error you have chosen. 

It then displays the plot of the Fourier series in an interective manner inspired in the 3blue1brown 
Youtube video [But what is a Fourier series?](https://www.youtube.com/watch?v=r6sGWTCMz2k&t=376s). 

https://github.com/MiquelNasarre/FourierS1/assets/124403865/fafab4c6-e3ca-4b7c-a5cb-ce5e8db3b79b

This video shows the opening example of the program.

## Requirements
- [Visual Studio](https://visualstudio.com) (Is mandatory for the program to be able to unite all the diffent C++ files)
- [SFML 2.6.0](https://www.sfml-dev.org/download.php) (The files needed are already included in the repository)

If you are not planning to build and run the program yourself and you do not need the code, you can also find the 
program executable with its associated files in the following [link](). Please note that the program will only work on 
Windows 7+ operating systems.

## How to Use it
First we clone the repository into your computer and open the solution file with Visual Studio, then compile the only
project inside and the App should be up and running, you should compile it in 32bit (x86).
Note: if you don't operate in Windows you might have some compatibility problems with the SFML renderer.

## How it Works
Since this program is my first take on creating Windows desktop applications, it mostly relies on its dependencies for 
all the techincal side. The program uses the SFML tools for all the window creation, interaction and rendering.

The programs structure is poorly organized, but after some effort is easy to follow its workings. The 
main file just calls the creation of a *Visual* class and runs the main Loop function, so all of the 
program revolves around such loop. Every frame it does an event check and renders everything accordingly, 
exiting the loop once the window is closed.

All of the widgets found on the menu bar are self designed, you can find them in the resources. Most 
of them have their own class and interections. This is the only program with such features since the 
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
using the mouse, once you're done press enter and the program will compute everything and display it as a 
function.

Now how does it work? Once you press the draw button, every time you are pressing the left button of the mouse 
and the mouse position changes it records a new location in an array. Once all the locations are recorded and you press 
enter the program considers every location is evenly spaced in time, and it distributes the function points accordingly. 
This assumption greatly simplifies the math in the paper, since all of the points are evenly spaced, every time we compute 
$t_{i+1} -t_i$ it will just be $\frac{2\pi}{N}$ where $N$ is the total number of points.

To see the points you have drawn just click on the shom button, then the points stored on the *Fourier* class will be displayed. 
You can also modify the number of points you want using the Points input on the top of the menu bar.

### Error Computations

### Coefficients View

### File Handling


