// Threaded two-dimensional Discrete FFT transform
// Tom Wells
// ECE8893 Project 2


#include <iostream>
#include <string>
#include <math.h>

#include "Complex.h"
#include "InputImage.h"

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.
Complex* imageData;//image array
Complex* W;
int nThreads=16; //number of worker threads
int imageWidth; //length of a row
int N; //number of points in whole image for the provided bit reversal function

pthread_mutex_t activeMutex;
pthread_mutex_t targetMutex;
int active;
pthread_barrier_t barrier;
using namespace std;

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = imageWidth; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}
void swap(Complex* h, int a, int b)
{
	if(a!=b)
	{
		Complex temp = h[a];
		h[a] = h[b];
		h[b] = temp;
	}
}
void Transpose(Complex*input, int n)
{
	for(int i = 0; i<n; i++)
	{
		for(int j = i+1; j<n ; j++)
			swap(input, n*i +j, n*j +i);//inplace matrix transpose 
	}
}
 void flipBit(Complex* input,bool debug)
{
	for(int i=0; i<imageWidth;i++)
	{
		int d = ReverseBits(i);
		if(i<d)
			swap(input, i, d);
		if(debug)
		cout<<i <<"goes to"<<d<<"\n";
	}
}
                  
void Transform1D(Complex* h, int O)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
	flipBit(h,false);
  for(int NP = 2; NP<=O; NP*=2)
	{
		for(int i = 0; i < O/NP; i++)
		{
			for(int j = 0; j < NP/2; j++)
			{	
				int index = i*(NP)+j;
				//int index = 0;
				int indexW = (j*O)/((double)(NP));//j*(O/(NP*2));//was j*NP/2

				Complex temp = h[index];
//				double arguement = (2*M_PI*j*O)/((double)O);
//				Complex w = Complex(cos(arguement), -sin(arguement));
				h[index] = h[index] + W[indexW]*h[index+(NP/2)];
				h[index+NP/2] = temp - W[indexW]*h[index+(NP/2)];
			
			
	
			}
		}
	}

//testing function
//	for(int i = 0; i<N; i++)
//	{
//		h[i] = h[(N-1)-i];
//	}
}
void* Transform2DThread(void *v)
{ // This is the thread startign point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete
	unsigned long rank = (unsigned long)v;
 	int targetDataSize = N/nThreads;
	int rowDataStart = rank *targetDataSize;
//	Complex rawArray[imageWidth];	
	for(int i = 0; i < targetDataSize/imageWidth; i++)
	{
		Transform1D(&imageData[rowDataStart+(imageWidth*i)], imageWidth);

	}
	pthread_mutex_lock(&activeMutex);
	active--;
	pthread_mutex_unlock(&activeMutex);
	pthread_barrier_wait(&barrier);
	for(int i = 0; i < targetDataSize/imageWidth; i++)
	{
		Transform1D(&imageData[rowDataStart+(imageWidth*i)], imageWidth);
	}
	pthread_mutex_lock(&activeMutex);
	active--;
	pthread_mutex_unlock(&activeMutex);
  return 0;
}

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.

	InputImage image(inputFN);  // Create the helper object for reading the image
	imageWidth = image.GetWidth();
	N = imageWidth *imageWidth; //total size  
	imageData = image.GetImageData();
	active = nThreads;
	pthread_mutex_init(&activeMutex,0);
//	pthread_mutex_init(&targetMutex,0);
	pthread_barrier_init(&barrier,0,nThreads+1);
	Complex Weights[imageWidth/2];
	W = Weights;
	for(int i = 0; i<imageWidth/2; i++)
	{
		double arguement = i*M_PI*2/((double) imageWidth);
			W[i] =  Complex(cos(arguement), -sin(arguement));
	
	}
	cout<< "Weights successful" << "\n";
 	 // Create 16 threads
	for(int i = 0; i<nThreads; i++)
	{
		pthread_t thread;
		pthread_create(&thread, 0, Transform2DThread,(void*)i);
	}
	
	while(active != 0)
	{}//spin
	image.SaveImageData("MyAfter1d.txt",imageData,imageWidth, image.GetHeight());
	Transpose(imageData,imageWidth);
	//reset active threads
	active = nThreads;
	
	//release barrier
	pthread_barrier_wait(&barrier);
	while(active !=0)
	{} //spin
	Transpose(imageData, imageWidth);
	image.SaveImageData("Tower-DFT2D.txt",imageData, imageWidth, image.GetHeight());
	
  
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  
  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  
