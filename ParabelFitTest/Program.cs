using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ParabelFitTest {
	class Program {
		static double focusPosition = 49;
		static double minPos = 0;
		static double maxPos = 52;
		static double sigma = 1000;
		static double amplitude = 50;
		static double offset = 50;
		static double noise = 0.5;

		static double startPos = 25;
		static int iterLimit = 5;

		static double blackLevel = 50;	//The level of least contrast
		static double maxIntensity = 255;
		static int rep = 1;

		static void Main(string[] args) {
			double[] pos = new double[rep]; 
			double mean, std;
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

		static double calculateMaxForFocus(double fPos)
		{
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
				for(int i = 0; i < 10; i++)
				{
					R += (vals[i, 1] - (A * vals[i, 0] * vals[i, 0] + B * vals[i, 0] + C)) * (vals[i, 1] - (A * vals[i, 0] * vals[i, 0] + B * vals[i, 0] + C));
				}
				R /= 10;
				Console.WriteLine("Pos {0}, R: {1}", fPos, R);
				iter++;
			}
			while((guess < minPos || guess > maxPos)&& iter < iterLimit);
			if(iter >= iterLimit)
			{
				Console.WriteLine("Focus peek detection failed");
				return 0.0;
			}
			//Console.WriteLine("fPos: {0}, Pos {1}, iter {2}, guess {3}", fPos, pos, iter, guess);
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

		static double intensity(double position)
		{
			Random r = new Random();
			return (offset + (r.NextDouble()-0.5) * noise + amplitude * Math.Exp(-((position - focusPosition) * (position - focusPosition)) / sigma)) * 2;
		}

		static double calcA(double x1, double x2, double x3, double x4, double xy, double x2y, double y) {
			return (x2y * (x2 - x1 * x1) + xy * (x1 * x2 - x3) + y * (x1 * x3 - x2 * x2)) / (x2 * x4 - x1 * x1 * x4 + 2 * x1 * x2 * x3 - x2 * x2 * x2 - x3 * x3);
		}

		static double calcB(double x1, double x2, double x3, double x4, double xy, double x2y, double y) {
			return (x2y * (x3 - x1 * x2) + xy * (x2 * x2 - x4) + y * (x1 * x4 - x2 * x3)) / (x2 * x4 - x1 * x1 * x4 + 2 * x1 * x2 * x3 - x2 * x2 * x2 - x3 * x3);
		}

		static double calcC(double x1, double x2, double x3, double x4, double xy, double x2y, double y) {
			return (x2y * (x1*x3 - x2 * x2) + xy * (x2 * x3 - x1*x4) + y * (x2 * x4 - x3 * x3)) / (x2 * x4 - x1 * x1 * x4 + 2 * x1 * x2 * x3 - x2 * x2 * x2 - x3 * x3);
		}

		static double calcMax(double x1, double x2, double x3, double x4, double xy, double x2y, double y){
      return 0.5 * (x2y*(x3-x1*x2)+xy*(x2*x2-x4)+y*(x1*x4-x2*x3))/(x2y*(x2-x1*x1)-xy*(x3-x1*x2)+y*(x1*x3-x2*x2));
    }
	}
}
