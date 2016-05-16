using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ParabelFitTest {
	class Program {
		static double focusPosition = 49;		//The position of the focus peak
		static double minPos = 0;						//Minimum position reachable
		static double maxPos = 52;					//Maximum position reachable
		static double sigma = 1000;					//Spread of the gaussian distribution
		static double amplitude = 50;				//Amplitude of the intensity signal
		static double offset = 50;					//Offset of the intensity signal
		static double noise = 0;						//Added noise to the signal

		static double startPos = 25;				//Start position of the successive search algorithm
		static int iterLimit = 5;						//Limit of initial iterations

		static double blackLevel = 50;			//The level of least contrast
		static double maxIntensity = 255;		//Maximum intensity
		static int rep = 1;									//Number of repetitions for statistical analysis

		static void Main(string[] args) {
			//Initialise variables for statistical analysis
			double[] pos = new double[rep]; 
			double mean, std;
			//Do the repeatet measurments for multiple positions of the focus peak and rep repetitons of the search
			for(int i = 30; i < 52; i++)
			{
				mean = 0.0;
				std = 0.0;
				for(int a = 0; a < rep; a++)
				{
					pos[a] = calculateMaxForFocus(i);
					mean += pos[a];
				}
				mean /= rep;
				for(int a = 0; a < rep; a++)
				{
					std += (pos[a] - mean) * (pos[a] - mean);
				}
				std /= rep;
				Console.WriteLine("Position: {0}, mean: {1}, varianz: {2}, noise {3}", i, mean, std, noise);
				//pos = findMaxDt(i);
				//Console.WriteLine("Position: {0}, calcualted dT: {1}, difference: {2}\n", i, pos, i - pos);
			}
			Console.ReadKey();

		}

		/// <summary>
		/// Calculates the position of the maximum for the given focus position fpos
		/// </summary>
		/// <param name="fPos">The position of the fcous peak</param>
		/// <returns>The calcualted posotion of the focus peak</returns>
		static double calculateMaxForFocus(double fPos)
		{
			//Initialise variables
			focusPosition = fPos;
			double x1 = 0.0;
			double x2 = 0.0;
			double x3 = 0.0;
			double x4 = 0.0;
			double xy = 0.0;
			double x2y = 0.0;
			double y = 0.0;
			
			int N = 0;
			double pos = startPos;
			double stepSize = 0.5;
			double guess = 0.0;
			double lastGuess = 0.0;
			int iter = 0;
			double [,] vals = new double[10,2];
			double A = 0.0, B = 0.0, C = 0.0, R = 0.0;
			//Calculate inital values
			//Repeat this step with differences in the start position and step size if the inital guess is to far off
			do
			{
				x1 = 0.0;
				x2 = 0.0;
				x3 = 0.0;
				x4 = 0.0;
				xy = 0.0;
				x2y = 0.0;
				y = 0.0;
				R = 0.0;
				stepSize = 0.7;// * (iter*2 + 1);
				if(guess > maxPos)
				{
					pos = startPos + iter * 7;
					if((maxPos - pos) / stepSize < 10)
						pos = maxPos - stepSize * 10;
				}
				else
				{
					pos = startPos;
				}
				N = 0;
				for(int i = 0; i < 10; i++)
				{
					vals[i, 0] = pos;
					vals[i, 1] =Math.Log(intensity(pos));
					x1 += pos;
					x2 += pos*pos;
					x3 += pos*pos*pos;
					x4 += pos*pos*pos*pos;
					xy += pos*Math.Log(intensity(pos));
					x2y += pos * pos * Math.Log(intensity(pos));
					y += Math.Log(intensity(pos));
					pos+=stepSize;
					N++;
				}
				x1 /= N;
				x2 /= N;
				x3 /= N;
				x4 /= N;
				xy /= N;
				x2y /= N;
				y /= N;
				if(fPos >= 40)
					focusPosition = fPos;
				guess = calcMax(x1, x2, x3, x4, xy, x2y, y);
				A = calcA(x1, x2, x3, x4, xy, x2y, y);
				B = calcB(x1, x2, x3, x4, xy, x2y, y);
				C = calcC(x1, x2, x3, x4, xy, x2y, y);
				//Calculate some sort of residual error
				for(int i = 0; i < 10; i++)
				{
					R += (vals[i, 1] - (A * vals[i, 0] * vals[i, 0] + B * vals[i, 0] + C)) * (vals[i, 1] - (A * vals[i, 0] * vals[i, 0] + B * vals[i, 0] + C));
				}
				R /= 10;
				//Console.WriteLine("Pos {0}, R: {1}", fPos, R);
				iter++;
			}
			while((guess < minPos || guess > maxPos)&& iter < iterLimit);
			if(iter >= iterLimit)
			{
				Console.WriteLine("Focus peek detection failed");
				return 0.0;
			}
			//Console.WriteLine("fPos: {0}, Pos {1}, iter {2}, guess {3}", fPos, pos, iter, guess);
			//Limit the guess
			if(guess < pos)
				guess = pos + 10 * stepSize;
			else if(guess > maxPos)
			{
				guess = maxPos;
			}
			lastGuess = guess;
			stepSize = (guess - pos) / 500;
			pos += stepSize;
			double diff = 1;
			//No iterate for some time to find a good guess of the focus peak
			while(Math.Abs(diff) > 0.0001 && N < 500)
			{
				x1 = x1 * N/(N+1) + pos / (N+1);
				x2 = x2 * N / (N + 1) + pos * pos / (N + 1);
				x3 = x3 * N / (N + 1) + pos * pos * pos / (N + 1);
				x4 = x4 * N / (N + 1) + pos * pos * pos * pos / (N + 1);
				xy = xy * N / (N + 1) + pos * Math.Log(intensity(pos)) / (N + 1);
				x2y = x2y * N / (N + 1) + pos * pos * Math.Log(intensity(pos)) / (N + 1);
				y = y * N / (N + 1) + Math.Log(intensity(pos)) / (N + 1);
				N++;
				guess = calcMax(x1, x2, x3, x4, xy, x2y, y);;
				diff = guess - lastGuess;
				lastGuess = guess;
				pos += stepSize;
			}
			//Console.WriteLine("Finished with diff = {0} and N = {1}", diff, N);
			return guess;
		}

		/// <summary>
		/// Alternative algorithm using the slope of the derivative, performs even worse
		/// </summary>
		/// <param name="fPos">The position of the focus peak</param>
		/// <returns>The calculated position of the focus peak</returns>
		static double findMaxDt(double fPos)
		{
			double pos = startPos;
			double stepSize = 0.5;
			int N = 40;

			double[] dY = new double[N];
			double a = 0.0, b=0.0, x, y;
			double x1 = 0.0;
			double x2 = 0.0;
			double xy = 0.0;
			double y1 = 0.0;
			for(int i = 0; i < N; i++)
			{
				dY[i] = (Math.Log(intensity(pos + stepSize) / intensity(pos))) / stepSize;
				pos += stepSize;
			}
			for(int i = 0; i < N; i++)
			{
				x = startPos + i*stepSize;
				y = dY[i];
				x1 += x;
				x2 += x*x;
				xy += x*y;
				y1 += y;
			}
			x1 /= N;
			x2 /= N;
			xy /= N;
			y1 /= N;

			a = (N*xy-N*x1*y1)/(N*x2-N*x1*x1);
			b = y1-a*x1;
			//a = (x2*y1-x1*xy)/(N*x2-x1*x1);
			//b = (N*xy-x1*y1)/(N*x2-x1*x1);
			return -b/a;
		}

		/// <summary>
		/// Function to return the intensity for a given position depending on offset, noise, amplitude, focusPosition and sigma
		/// </summary>
		/// <param name="position">The position to calcualte the intensity for</param>
		/// <returns>The intensity</returns>
		static double intensity(double position)
		{
			Random r = new Random();
			return (offset + (r.NextDouble()-0.5) * noise + amplitude * Math.Exp(-((position - focusPosition) * (position - focusPosition)) / sigma)) * 2;
		}

		/// <summary>
		/// Calculate A from the quadratic formular Ax^2+Bx+c
		/// </summary>
		/// <param name="x1">1/N sum x</param>
		/// <param name="x2">1/N sum x^2</param>
		/// <param name="x3">1/N sum x^3</param>
		/// <param name="x4">1/N sum x^4</param>
		/// <param name="xy">1/N sum xy</param>
		/// <param name="x2y">1/N sum x^2y</param>
		/// <param name="y">1/N sum y</param>
		/// <returns>A</returns>
		static double calcA(double x1, double x2, double x3, double x4, double xy, double x2y, double y) {
			return (x2y * (x2 - x1 * x1) + xy * (x1 * x2 - x3) + y * (x1 * x3 - x2 * x2)) / (x2 * x4 - x1 * x1 * x4 + 2 * x1 * x2 * x3 - x2 * x2 * x2 - x3 * x3);
		}

		/// <summary>
		/// Calculate B from the quadratic formular Ax^2+Bx+c
		/// </summary>
		/// <param name="x1">1/N sum x</param>
		/// <param name="x2">1/N sum x^2</param>
		/// <param name="x3">1/N sum x^3</param>
		/// <param name="x4">1/N sum x^4</param>
		/// <param name="xy">1/N sum xy</param>
		/// <param name="x2y">1/N sum x^2y</param>
		/// <param name="y">1/N sum y</param>
		/// <returns>B</returns>
		static double calcB(double x1, double x2, double x3, double x4, double xy, double x2y, double y) {
			return (x2y * (x3 - x1 * x2) + xy * (x2 * x2 - x4) + y * (x1 * x4 - x2 * x3)) / (x2 * x4 - x1 * x1 * x4 + 2 * x1 * x2 * x3 - x2 * x2 * x2 - x3 * x3);
		}

		/// <summary>
		/// Calculate C from the quadratic formular Ax^2+Bx+c
		/// </summary>
		/// <param name="x1">1/N sum x</param>
		/// <param name="x2">1/N sum x^2</param>
		/// <param name="x3">1/N sum x^3</param>
		/// <param name="x4">1/N sum x^4</param>
		/// <param name="xy">1/N sum xy</param>
		/// <param name="x2y">1/N sum x^2y</param>
		/// <param name="y">1/N sum y</param>
		/// <returns>C</returns>
		static double calcC(double x1, double x2, double x3, double x4, double xy, double x2y, double y) {
			return (x2y * (x1*x3 - x2 * x2) + xy * (x2 * x3 - x1*x4) + y * (x2 * x4 - x3 * x3)) / (x2 * x4 - x1 * x1 * x4 + 2 * x1 * x2 * x3 - x2 * x2 * x2 - x3 * x3);
		}

		/// <summary>
		/// Calculate the maximum of the parabola from the fit data
		/// </summary>
		/// <param name="x1">1/N sum x</param>
		/// <param name="x2">1/N sum x^2</param>
		/// <param name="x3">1/N sum x^3</param>
		/// <param name="x4">1/N sum x^4</param>
		/// <param name="xy">1/N sum xy</param>
		/// <param name="x2y">1/N sum x^2y</param>
		/// <param name="y">1/N sum y</param>
		/// <returns>Position of the maximum</returns>
		static double calcMax(double x1, double x2, double x3, double x4, double xy, double x2y, double y){
      return 0.5 * (x2y*(x3-x1*x2)+xy*(x2*x2-x4)+y*(x1*x4-x2*x3))/(x2y*(x2-x1*x1)-xy*(x3-x1*x2)+y*(x1*x3-x2*x2));
    }
	}
}
