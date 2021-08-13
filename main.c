//Volkan Can Bacaksız
//365259
//1. Öğretim
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define INPUT "PianoC6.wav"
#define OUTPUT "output.wav"
#define MAXCOEFFSIZE 22 		// Chebyshev max coeff size=>20 + need for function=>2 = 22
#define PI 3.141592				// PI number
#define TRUE 1
#define FALSE 0
#define MAXTXTFILENAMELENGTH 40 // Maximum name length for txt file

struct HEADER
{
	unsigned char riff[4];			  // RIFF string
	unsigned int overallSize;		  // overall size of file in bytes
	unsigned char wave[4];			  // WAVE string
	unsigned char fmtChunkMarker[4];  // fmt string with trailing null char
	unsigned int lengthOfFmt;		  // length of the format data
	unsigned int formatType;		  // format type. 1-PCM, 3- IEEE float, 6 - 8bit A law, 7 - 8bit mu law
	unsigned int channels;			  // num of channels
	unsigned int sampleRate;		  // sampling rate (blocks per second)
	unsigned int byterate;			  // SampleRate * NumChannels * BitsPerSample/8
	unsigned int blockAlign;		  // NumChannels * BitsPerSample/8
	unsigned int bitsPerSample;		  // bits per sample, 8- 8bits, 16- 16 bits etc
	unsigned char dataChunkHeader[4]; // DATA string or FLLR string
	unsigned int dataSize;			  // NumSamples * NumChannels * BitsPerSample/8 - size of the next chunk that will be read
};

int userRequestControl(float *aCoeffs, float *bCoeffs, int coeffSize);

float *getCoeffsFromTxt(const char *txtFilename, int coeffSize, float *coeffsArray);

int userGainChoice(float *aCoeffs, float *bCoeffs, int coeffSize);

int findSignalSize(const char *filename);

void readWavFile(const char *filename, float *signal, const int signalSize);

int calculateCoeffs(float *aCoeffs, float *bCoeffs, int coeffSize, int gainRequest);

float *needForCalculateCoefficient(float cutOffFrequency, int passType, int percentRipple, int numOfPole, int numOfStep, float *needElements);

void filter(const float *signalIn, float *signalOut, const float *aCoeffs,
			const float *bCoeffs, const int signalSize, int coeffSize);

void writeWavFile(const char *filename, float *signal, const int signalSize);

struct HEADER header;
struct HEADER tempHeader;

int main()
{
	const char *readFileName = INPUT;
	int signalSize = findSignalSize(readFileName);

	float *signalIn = (float *)malloc(sizeof(float) * signalSize);
	float *signalOut = (float *)malloc(sizeof(float) * signalSize);

	float *aCoeffs = (float *)malloc(sizeof(float) * MAXCOEFFSIZE);
	float *bCoeffs = (float *)malloc(sizeof(float) * MAXCOEFFSIZE);
	int coeffSize = 0;
	coeffSize = userRequestControl(aCoeffs, bCoeffs, coeffSize);
	if (coeffSize != 0)
	{
		readWavFile(readFileName, signalIn, signalSize);

		filter(signalIn, signalOut, aCoeffs, bCoeffs, signalSize, coeffSize);

		const char *writeFileName = OUTPUT;
		writeWavFile(writeFileName, signalOut, signalSize);

		printf("\nWav File Name --> %s\n\n", writeFileName);
		free(signalIn);
		free(signalOut);
		free(aCoeffs);
		free(bCoeffs);
	}
	else
	{
		printf("\nApplication Quit...\n\n");
	}
	system("pause");
	return 0;
}

//get input from user for how the program works
int userRequestControl(float *aCoeffs, float *bCoeffs, int coeffSize)
{
	char request;
	char prevRequest;
	int status = TRUE;
	while (status != FALSE)
	{
		if (prevRequest != '\n')
		{
			printf("\n");
			printf("--------------------------------- \n");
			printf("If U Want To Give The Coeffs Please Enter 1\n");
			printf("Create Coeffs For Filter Please Enter 2\n");
			printf("Quit The Program Please Enter 0\n");
		}
		scanf("%c", &request);
		switch (request)
		{
		case '1':
			printf("\nEnter Text Name For A Coeffs\n");
			printf("Example ---> aCoeffs.txt\n\n");
			char aCoeffsFile[MAXTXTFILENAMELENGTH];
			scanf("%s", &aCoeffsFile);
			printf("\nEnter Text Name For B Coeffs\n");
			printf("Example ---> bCoeffs.txt\n\n");
			char bCoeffsFile[MAXTXTFILENAMELENGTH];
			scanf("%s", &bCoeffsFile);
			printf("\nEnter Coeff Size\n\n");
			scanf("%d", &coeffSize);
			aCoeffs = getCoeffsFromTxt(aCoeffsFile, coeffSize, aCoeffs);
			bCoeffs = getCoeffsFromTxt(bCoeffsFile, coeffSize, bCoeffs);
			status = FALSE;
			break;
		case '2':
			coeffSize = userGainChoice(aCoeffs, bCoeffs, coeffSize);
			status = FALSE;
			break;
		case '0':
			return 0;
			break;
		default:
			break;
		}
		prevRequest = request;
	}
	return coeffSize;
}

//write for user want to give input from txt
float *getCoeffsFromTxt(const char *txtFilename, int coeffSize, float *coeffsArray)
{
	FILE *coeffsTxt;
	coeffsTxt = fopen(txtFilename, "rb");
	float coeff;

	for (int i = 0; i < coeffSize; i++)
	{
		fscanf(coeffsTxt, "%f", &coeff);
		coeffsArray[i] = coeff;
		printf("%f\n", coeffsArray[i]);
	}
	fclose(coeffsTxt);

	return coeffsArray;
}

//get input from user for calculate coeffs with gain or without gain
int userGainChoice(float *aCoeffs, float *bCoeffs, int coeffSize)
{
	char request;
	char prevRequest;
	int status = TRUE;
	while (status != FALSE)
	{
		if (prevRequest != '\n')
		{
			printf("\n");
			printf("--------------------------------- \n");
			printf("Create Coeffs With Gain Please Enter 1\n");
			printf("Create Coeffs Without Gain Please Enter 2\n");
			printf("Quit The Program Please Enter 0\n");
		}
		scanf("%c", &request);
		switch (request)
		{
		case '1':
			coeffSize = calculateCoeffs(aCoeffs, bCoeffs, coeffSize, TRUE);
			status = FALSE;
			break;
		case '2':
			coeffSize = calculateCoeffs(aCoeffs, bCoeffs, coeffSize, FALSE);
			status = FALSE;
			break;
		case '0':
			status = FALSE;
			return 0;
		default:
			break;
		}
		prevRequest = request;
	}
	return coeffSize;
}

//find signal size for transactions
int findSignalSize(const char *filename)
{
	unsigned int buffer4;
	FILE *inFile;
	unsigned int xChannels = 2;
	long numSamples;

	inFile = fopen(filename, "rb");

	if (inFile == NULL)
	{
		printf("Error opening input wav file\n");
		system("pause");
		exit(1);
	}

	// read header parts
	fread(tempHeader.riff, sizeof(tempHeader.riff), 1, inFile);
	fread((unsigned char *)&tempHeader.overallSize, 4, 1, inFile);
	printf("Overall Size of input WAV file is %ld\n", tempHeader.overallSize);
	fread(tempHeader.wave, sizeof(tempHeader.wave), 1, inFile);
	fread(tempHeader.fmtChunkMarker, sizeof(tempHeader.fmtChunkMarker), 1, inFile);
	fread((unsigned char *)&tempHeader.lengthOfFmt, sizeof(buffer4), 1, inFile);
	fread((unsigned char *)&tempHeader.formatType, 2, 1, inFile);
	fread((unsigned char *)&tempHeader.channels, 2, 1, inFile);
	printf("Number of Channels in input Wav file is %d\n", tempHeader.channels);
	fread((unsigned char *)&tempHeader.sampleRate, 4, 1, inFile);
	printf("Sample Rate used in input WAV file is %d Hz\n", tempHeader.sampleRate);
	fread((unsigned char *)&tempHeader.byterate, 4, 1, inFile);
	fread((unsigned char *)&tempHeader.blockAlign, 2, 1, inFile);
	fread((unsigned char *)&tempHeader.bitsPerSample, 2, 1, inFile);
	fread(tempHeader.dataChunkHeader, sizeof(tempHeader.dataChunkHeader), 1, inFile);
	fread((unsigned char *)&tempHeader.dataSize, 4, 1, inFile);
	printf("Data Size of input WAV file is %ld\n", tempHeader.dataSize);

	numSamples = (8 * tempHeader.dataSize) / (tempHeader.channels * tempHeader.bitsPerSample);

	printf("Number of Samples in input WAV file: %ld bytes\n\n", numSamples);
	fclose(inFile);
	return numSamples;
}

//read wav file
void readWavFile(const char *filename, float *signal, const int signalSize)
{
	unsigned int buffer4;
	unsigned short int buffer2;
	FILE *inFile;
	unsigned int xChannels = 2;
	long numSamples;
	long int i = 0;
	short dataBuffer;

	inFile = fopen(filename, "rb");

	if (inFile == NULL)
	{
		printf("Error opening input wav file\n");
		system("pause");
		exit(1);
	}

	//read wav file header
	fread(header.riff, sizeof(header.riff), 1, inFile);
	fread((unsigned char *)&header.overallSize, 4, 1, inFile);
	printf("\nOverall Size of input WAV file is %ld\n", header.overallSize);
	fread(header.wave, sizeof(header.wave), 1, inFile);
	fread(header.fmtChunkMarker, sizeof(header.fmtChunkMarker), 1, inFile);
	fread((unsigned char *)&header.lengthOfFmt, sizeof(buffer4), 1, inFile);
	fread((unsigned char *)&header.formatType, 2, 1, inFile);
	fread((unsigned char *)&header.channels, 2, 1, inFile);
	printf("Number of Channels in input Wav file is %d\n", header.channels);
	fread((unsigned char *)&header.sampleRate, 4, 1, inFile);
	printf("Sample Rate used in input WAV file is %d Hz\n", header.sampleRate);
	fread((unsigned char *)&header.byterate, 4, 1, inFile);
	fread((unsigned char *)&header.blockAlign, 2, 1, inFile);
	fread((unsigned char *)&header.bitsPerSample, 2, 1, inFile);
	fread(header.dataChunkHeader, sizeof(header.dataChunkHeader), 1, inFile);
	fread((unsigned char *)&header.dataSize, 4, 1, inFile);
	printf("Data Size of input WAV file is %ld\n", header.dataSize);

	numSamples = (8 * header.dataSize) / (header.channels * header.bitsPerSample);

	printf("Number of Samples in input WAV file: %ld bytes\n\n", numSamples);

	for (i = 0; i < signalSize; i++)
	{

		for (xChannels = 0; xChannels < header.channels; xChannels++)
		{
			fread((unsigned char *)&dataBuffer, sizeof(dataBuffer), 1, inFile);
			if (xChannels == 0)
			{
				signal[i] = dataBuffer / 32768.0;
			}
		}
	}
	fclose(inFile);
	return;
}

//calculate coeffs without gain
//Reference : dsp book chapter 20 (page 333 - 342)
int calculateCoeffs(float *aCoeffs, float *bCoeffs, int coeffSize, int gainRequest)
{
	float *tempACoeffs = (float *)malloc(sizeof(float) * MAXCOEFFSIZE);
	float *tempBCoeffs = (float *)malloc(sizeof(float) * MAXCOEFFSIZE);

	float cutOffFrequency = 0.0f; // cutoff frequency => Example _ 48Khz wav file --> 0.5 => 48KHz , 0.1 =>12KHz
	int passType = 0;			  // 0 low pass - 1 high pass
	int percentRipple = 0;		  // percent ripple value (if percent ripple higher than 0, waves dissipate)
	int numOfPole = 0;			  // like num of coeffs
	//Get input from user and control
	int status = FALSE;
	while(status != TRUE)
	{
		printf("\nEnter cutoff frequency 0 to .5\n");
		scanf("%f", &cutOffFrequency);
		if(cutOffFrequency >= 0.0 &&  cutOffFrequency <= 0.5)
			status = TRUE;
		else
			printf("Input Error!!");
	}
	status = FALSE;
	while(status != TRUE)
	{
		printf("\nEnter 0 For Low Pass - 1 For High Pass\n");
		scanf("%d", &passType);
		if(passType == 0 || passType == 1)
			status = TRUE;
		else
			printf("Input Error!!");
	}
	status = FALSE;
	while(status != TRUE)
	{
		printf("\nEnter percent ripple (0 to 29)\n");
		scanf("%d", &percentRipple);
		if(percentRipple >= 0 && percentRipple <= 29)
			status = TRUE;
		else
			printf("Input Error!!");
	}
	status = FALSE;
	while (status != TRUE)
	{
		printf("\nEnter number of poles ( 2, 4, ....20)\n");
		scanf("%d", &numOfPole);
		if (numOfPole % 2 == 0 && numOfPole > 0 && numOfPole < 21)
			status = TRUE;
		else
			printf("Input Error!!");
	}

	//initilize coeffs Array
	for (int i = 0; i < MAXCOEFFSIZE; i++)
	{
		aCoeffs[i] = 0.0f;
		bCoeffs[i] = 0.0f;
	}

	aCoeffs[2] = 1.0f;
	bCoeffs[2] = 1.0f;
	float *needElements;
	needElements = (float *)malloc(sizeof(float) * MAXCOEFFSIZE);

	//Substitution and Creation of Coeffs
	for (int p = 1; p < numOfPole / 2; p++)
	{
		needElements = needForCalculateCoefficient(cutOffFrequency, passType, percentRipple, numOfPole, p, needElements);

		for (int i = 0; i < MAXCOEFFSIZE; i++)
		{
			tempACoeffs[i] = aCoeffs[i];
			tempBCoeffs[i] = bCoeffs[i];
		}

		for (int i = 2; i < MAXCOEFFSIZE; i++)
		{
			aCoeffs[i] = needElements[0] * tempACoeffs[i] + needElements[1] * tempACoeffs[i - 1] + needElements[2] * tempACoeffs[i - 2];
			bCoeffs[i] = tempBCoeffs[i] - needElements[3] * tempBCoeffs[i - 1] - needElements[4] * tempBCoeffs[i - 2];
		}
	}

	bCoeffs[2] = 0;
	for (int i = 0; i < MAXCOEFFSIZE - 2; i++)
	{
		aCoeffs[i] = aCoeffs[i + 2];
		bCoeffs[i] = -(bCoeffs[i + 2]);
	}

	//Gain
	if (gainRequest)
	{
		float gainUpper = 0.0f;
		float gainLower = 0.0f;
		for (int i = 0; i < MAXCOEFFSIZE - 2; i++)
		{
			if (passType == 0)
			{
				gainUpper += aCoeffs[i];
				gainLower += bCoeffs[i];
			}
			else
			{
				gainUpper += aCoeffs[i] * pow(-1, i);
				gainLower += bCoeffs[i] * pow(-1, i);
			}
		}

		float gain = gainUpper / (1 - gainLower);

		for (int i = 0; i < MAXCOEFFSIZE - 2; i++)
		{
			aCoeffs[i] = aCoeffs[i] / gain;
		}
	}

	coeffSize = numOfPole + 1;

	free(tempACoeffs);
	free(tempBCoeffs);
	free(needElements);
	return coeffSize;
}

//Calculation of necessary elements to calculate coefficient
//Reference : dsp book chapter 20 (page 333 - 342)
float *needForCalculateCoefficient(float cutOffFrequency, int passType, int percentRipple, int numOfPole, int numOfStep, float *needElements)
{
	float rp = -(cos(PI / (numOfPole * 2) + (numOfStep - 1) * (PI / numOfPole)));
	float ip = sin(PI / (numOfPole * 2) + (numOfStep - 1) * (PI / numOfStep));

	float es = 0.0f;
	float vx = 0.0f;
	float kx = 0.0f;
	if (percentRipple != 0)
	{
		es = sqrt(pow((float)100 / (float)(100 - percentRipple), 2) - 1);
		vx = ((float)1 / (float)numOfPole) * log((float)(1 / es) + sqrt((1 / pow(es, 2) + 1)));
		kx = ((float)1 / (float)numOfPole) * log((float)(1 / es) + sqrt((1 / pow(es, 2) - 1)));
		kx = (exp(kx) + exp(-(kx))) * 0.5f;
		rp *= ((exp(vx) - exp(-(vx))) * 0.5f) / kx;
		ip *= ((exp(vx) + exp(-(vx))) * 0.5f) / kx;
	}
	// printf("%f,   %f,   %f,    %f\n\n\n\n", rp, ip, es, vx);

	float t = 2 * tan(0.5);
	float w = 2 * PI * cutOffFrequency;
	float m = pow(rp, 2) + pow(ip, 2);
	float d = 4 - 4 * rp * t + m * pow(t, 2);
	float X0 = pow(t, 2) / d;
	float X1 = 2 * pow(t, 2) / d;
	float X2 = pow(t, 2) / d;
	float Y1 = (8 - 2 * m * pow(t, 2)) / d;
	float Y2 = (-4 - 4 * rp * t - m * pow(t, 2)) / d;
	// printf("%f,   %f,   %f,    %f\n\n\n\n", t, w, m, d);
	// printf("%f,   %f,   %f,    %f\n\n\n\n", X0, X1, X2, Y1);

	float k = 0;
	if (passType == 1)
	{
		k = -cos(0.5 * w + 0.5) / cos(0.5 * w - 0.5);
	}
	else
	{
		k = sin(0.5 - 0.5 * w) / sin(0.5 + 0.5 * w);
	}

	d = 1 + Y1 * k - Y2 * pow(k, 2);
	float A0 = (X0 - X1 * k + X2 * pow(k, 2)) / d;
	float A1 = (-2 * X0 * k + X1 + X1 * pow(k, 2) - 2 * X2 * k) / d;
	float A2 = (X0 * pow(k, 2) - X1 * k + X2) / d;
	float B1 = (2 * k + Y1 + Y1 * pow(k, 2) - 2 * Y2 * k) / d;
	float B2 = (-(pow(k, 2)) - Y1 * k + Y2) / d;

	if (passType == 1)
	{
		A1 = -(A1);
		B1 = -(B1);
	}

	needElements[0] = A0;
	needElements[1] = A1;
	needElements[2] = A2;
	needElements[3] = B1;
	needElements[4] = B2;
	// printf("%f,   %f,   %f,    %f\n\n\n\n", A0, A1, A2, B1);

	return needElements;
}

//Applying filter
void filter(const float *signalIn, float *signalOut, const float *aCoeffs,
			const float *bCoeffs, const int signalSize, int coeffSize)
{
	int i;
	int j;

	for (i = 0; i < signalSize; i++)
	{
		signalOut[i] = 0;
		for (j = 0; j < coeffSize; j++)
		{
			if ((i - j) >= 0)
			{
				signalOut[i] += signalIn[i - j] * aCoeffs[j] + signalOut[i - j] * bCoeffs[j];
			}
		}
	}

	return;
}

//write wav file
void writeWavFile(const char *filename, float *signal, const int signalSize)
{
	unsigned int buffer4;
	unsigned short int buffer2;

	FILE *outfile;

	long numSamples;
	long sizeOfEachSample;
	long i = 0;
	short int dataBuffer;
	outfile = fopen(filename, "wb");

	if (outfile == NULL)
	{
		printf("Error opening output WAV file\n");
		system("pause");
		exit(1);
	}

	//create header
	printf("Signal Size to be output is %d Samples\n", signalSize);
	strcpy(header.riff, "RIFF");
	header.overallSize = tempHeader.overallSize;
	strcpy(header.wave, "WAVE");
	strcpy(header.fmtChunkMarker, "fmt ");
	header.lengthOfFmt = tempHeader.lengthOfFmt;
	header.formatType = tempHeader.formatType;
	header.channels = tempHeader.channels;
	header.sampleRate = tempHeader.sampleRate;
	header.byterate = tempHeader.byterate;
	header.blockAlign = tempHeader.blockAlign;
	header.bitsPerSample = tempHeader.bitsPerSample;
	strcpy(header.dataChunkHeader, "data");
	header.dataSize = tempHeader.dataSize;

	//writing wav file header parts
	fwrite((unsigned char *)&header.riff, sizeof(header.riff), 1, outfile);
	fwrite((unsigned char *)&header.overallSize, sizeof(buffer4), 1, outfile);
	fwrite((unsigned char *)&header.wave, sizeof(header.wave), 1, outfile);
	fwrite((unsigned char *)&header.fmtChunkMarker, sizeof(header.fmtChunkMarker), 1, outfile);
	fwrite((unsigned char *)&header.lengthOfFmt, sizeof(buffer4), 1, outfile);
	fwrite((unsigned char *)&header.formatType, sizeof(buffer2), 1, outfile);
	fwrite((unsigned char *)&header.channels, sizeof(buffer2), 1, outfile);
	fwrite((unsigned char *)&header.sampleRate, sizeof(buffer4), 1, outfile);
	fwrite((unsigned char *)&header.byterate, sizeof(buffer4), 1, outfile);
	fwrite((unsigned char *)&header.blockAlign, sizeof(buffer2), 1, outfile);
	fwrite((unsigned char *)&header.bitsPerSample, sizeof(buffer2), 1, outfile);
	fwrite((unsigned char *)&header.dataChunkHeader, sizeof(header.dataChunkHeader), 1, outfile);
	fwrite((unsigned char *)&header.dataSize, sizeof(buffer4), 1, outfile);

	for (i = 0; i < (signalSize); i++)
	{
		dataBuffer = (short int)(((*(signal + i))) * 32768);
		fwrite((unsigned char *)(&dataBuffer), 2, 1, outfile); //Data -> Left Channel
		fwrite((unsigned char *)(&dataBuffer), 2, 1, outfile); //Data -> Right Channel
	}
}