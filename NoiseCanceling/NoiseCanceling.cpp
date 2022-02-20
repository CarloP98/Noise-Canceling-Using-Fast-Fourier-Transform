#include <iostream>
#include <fstream>
#include <complex>
#include <valarray>
#include<ctime>
#include <stdlib.h>
#include <math.h> 
#include "filter.h"
#pragma warning(disable:4996)

const double PI = 3.141592653589793238460;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
typedef struct Header
{
	char            ChunkID[4];
	unsigned int    ChunkSize;
	char            Format[4];
	char            Subchunk1ID[4];
	unsigned int    Subchunk1Size;
	unsigned short  AudioFormat;
	unsigned short  NumChannels;
	unsigned int    SampleRate;
	unsigned int    ByteRate;
	unsigned short  BlockAlign;
	unsigned short  BitsPerSample;
	char            Subchunk2ID[4];
	unsigned int    Subchunk2Size;
}Header;

void fft(CArray& x)
{
	const size_t N = x.size();
	if (N <= 1) return;
	// divide
	CArray even = x[std::slice(0, N / 2, 2)];
	CArray  odd = x[std::slice(1, N / 2, 2)];
	// conquer
	fft(even);
	fft(odd);
	// combine
	for (size_t k = 0; k < N / 2; ++k)
	{
		Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
		x[k] = even[k] + t;
		x[k + N / 2] = even[k] - t;
	}
}

double* hamming(int size) {
	int m; int n;
	double* out = new double[size] { 0 };
	if (size % 2 == 0)
		m = n = size / 2;
	else {
		m = ceil(size / 2);
		n = floor(size / 2);
	}
	for (int x = 0; x < size; x++)
		out[x] = 0.54 - 0.46 * cos(2 * PI * (x) / ((double)size - 1));
	return out;
}

double* bandreject(int size, int fc1, int fc2, int fs) {
	double Fc1 = fc1 / double(fs);
	double Fc2 = fc2 / double(fs);
	double* A = new double[size] { 0 };
	double* B = new double[size] { 0 };
	double* output = new double[size] { 0 };
	double* window = hamming(size);
	for (int i = 0; i < size; i++) {
		if (2 * i == size)
			A[i] = 2 * PI * Fc1;
		else
			A[i] = sin(2 * PI * Fc1 * (i - ((double)size / 2))) / (i - ((double)size / 2));
		A[i] *= window[i];
	}
	for (int i = 0; i < size; i++) {
		if (2 * i == size)
			B[i] = 2 * PI * Fc2;
		else
			B[i] = sin(2 * PI * Fc2 * (i - ((double)size / 2))) / (i - ((double)size / 2));
		B[i] *= window[i];
	}
	double Asum = 0;
	double Bsum = 0;
	for (int x = 0; x < size; x++) {
		Asum += A[x];
		Bsum += B[x];
	}
	for (int x = 0; x < size; x++) {
		B[x] /= Asum;
		A[x] /= Bsum;
	}
	for (int x = 0; x < size; x++) {
		B[x] *= -1;
	}
	B[int(size / 2)] += 1;
	for (int x = 0; x < size; x++) {
		output[x] = (A[x] + B[x]) * .7;
	}
	return output;
}

int main()
{
	FILE* inFile = fopen("original.wav", "rb");
	if (!inFile) { std::cout << "cannot find file" << std::endl; return 0; }
	FILE* outFile = fopen("modified.wav", "wb");
	std::cout << "reading file..." << std::endl;
	Header FileHeader;
	Header OutputHeader;
	fread(&FileHeader, sizeof(FileHeader), 1, inFile);
	OutputHeader = FileHeader;
	fwrite(&OutputHeader, sizeof(FileHeader), 1, outFile);

	short* u = new short[FileHeader.Subchunk2Size / 2];
	int MAX = 5000;
	Complex* DATA = new Complex[MAX];

	for (int x = 0; x < FileHeader.Subchunk2Size / 2; x++)
	{
		fread(&u[x], sizeof(short), 1, inFile);
		if (x < MAX)
			DATA[x] = (double)u[x];
	}

	//Dominant Spectral Component
	std::cout << "finding dominant spectral component..." << std::endl;
	CArray data(DATA, MAX);
	fft(data);
	double maxSpec = (data[0].real()) * (data[0].real()) + (data[0].imag()) * (data[0].imag());
	double tmp = 0;
	int maxIndex = 0;

	for (int j = 1; j < MAX; j++)
	{
		tmp = (data[j].real()) * (data[j].real()) + (data[j].imag()) * (data[j].imag());
		if (tmp > maxSpec)
		{
			maxSpec = tmp;
			maxIndex = j;
		}
	}
	int DSC = maxIndex * FileHeader.SampleRate * FileHeader.NumChannels / MAX;

	std::cout << "component fount at " << DSC << "Hz" << std::endl;

	int ss = 1025;
	short* w = new short[(FileHeader.Subchunk2Size / 2) + ss - 1];

	std::cout << "Generating bandstop filter coefficients..." << std::endl;
	double* test = bandreject(ss, DSC-100, DSC+100, FileHeader.SampleRate*2);

	std::cout << "Applying filter..." << std::endl;

	for (int i = 0; i < ((FileHeader.Subchunk2Size / 2) + ss - 1); i++)
		for (int j = 0; j < ss; j++)
			if (i - j >= 0 && i - j < FileHeader.Subchunk2Size / 2)
				w[i] += u[i - j] * test[j];

	std::cout << "Writing output to new wav file..." << std::endl;
	for (int x = 0; x <= FileHeader.Subchunk2Size / 2; x++)
		fwrite(&w[x], sizeof(short), 1, outFile);

	fclose(inFile);
	fclose(outFile);
	std::cout << "***DONE***" << std::endl;

	return 0;
}