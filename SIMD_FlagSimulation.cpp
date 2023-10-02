#include "precomp.h"

VerletFlag::VerletFlag( int2 location, Surface* pattern )
{
	width = pattern->width;
	height = pattern->height;
	polePos = make_float2( location );
	posx = new float[width * height];
	posy = new float[width * height];
	prevPosx = new float[width * height];
	prevPosy = new float[width * height];
	color = new uint[width * height];
	backup = new uint[width * height * 4];
	memcpy( color, pattern->pixels, width * height * 4 );
	for (int x = 0; x < width; x++) for (int y = 0; y < height; y++)
	{
		posx[y + x * height] = location.x - x * 1.2f;
		posy[y + x * height] = y * 1.2f + location.y;
	}	
	memcpy( prevPosx, posx, width * height * 4 );
	memcpy(prevPosy, posy, width * height * 4);
}


void VerletFlag::Draw()
{
	for (int x = 0; x < width; x++) for (int y = 0; y < height; y++)
	{
		// plot each flag vertex bilinearly
		int index = y + x * height;
		int colorIndex = x + y * width;
		int2 intPos = make_int2( posx[index], posy[index] );
		backup[index * 4 + 0] = MyApp::map.bitmap->Read( intPos.x, intPos.y );
		backup[index * 4 + 1] = MyApp::map.bitmap->Read( intPos.x + 1, intPos.y );
		backup[index * 4 + 2] = MyApp::map.bitmap->Read( intPos.x, intPos.y + 1 );
		backup[index * 4 + 3] = MyApp::map.bitmap->Read( intPos.x + 1, intPos.y + 1 );
		hasBackup = true;
		MyApp::map.bitmap->PlotBilerp( posx[index], posy[index], color[colorIndex] );
	}
}

bool VerletFlag::Tick()
{
	uint64_t start = __rdtsc();
	int index;
	__m256 deltax8;
	__m256 deltay8;
	// move vertices
	for (int x = 0; x < width; x++) for (int y = 0; y < height/8; y++)
	{
		index = y + x * height / 8;
		deltax8 = _mm256_sub_ps(posx8[index], prevPosx8[index]);
		deltay8 = _mm256_sub_ps(posy8[index], prevPosy8[index]);
		prevPosx8[index] = posx8[index];
		prevPosy8[index] = posy8[index];
		posx8[index] = _mm256_add_ps(posx8[index], deltax8);
		posy8[index] = _mm256_add_ps(posy8[index], deltay8);
	}
	// apply forces
	float windForce = 0.1f + 0.05f * RandomFloat();
	float2 wind = windForce * normalize(make_float2(-1, (RandomFloat() * 0.5f) - 0.25f));
	__m256 windx8 = _mm256_set1_ps(wind.x);
	__m256 windy8 = _mm256_set1_ps(wind.y);
	for (int x = 1; x < width; x++) for (int y = 0; y < height/8; y++)
	{
		index = y + x * height/8;
		posx8[index] = _mm256_add_ps(posx8[index], windx8);
		posy8[index] = _mm256_add_ps(posy8[index], windy8);
	}
	float2 nudge;
	for (int x = 1; x < width; x++) for (int y = 0; y < height; y++)
	{
		index = y + x * height;
		if ((RandomUInt() & 31) == 31)
		{
			// small chance of a random nudge to add a bit of noise to the animation
			nudge = make_float2(RandomFloat() - 0.5f, RandomFloat() - 0.5f);
			posx[index] += nudge.x;
			posy[index] += nudge.y;
		}
	}
	__m256 aHalf8 = _mm256_set1_ps(0.5f);
	__m256 maxLength8 = _mm256_set1_ps(1.15f);
	__m256 maxLengthsquared8 = _mm256_mul_ps(maxLength8, maxLength8);
	union { int counter[8]; __m256 counter8; };
	__m256 rightx8;
	__m256 righty8;
	__m256 halfExcessx8;
	__m256 halfExcessy8;
	__m256 sqrLength8;
	__m256 mask;
	__m256 lengthInv8;
	// constraints: limit distance
	for (int i = 0; i < 25; i++)
	{
		counter8 = _mm256_setzero_ps();
		for (int x = 1; x < width; x++) for (int y = 0; y < height/8; y++)
		{
			index = y + x * height/8;

			rightx8 = _mm256_sub_ps(posx8[index - height / 8], posx8[index]);
			righty8 = _mm256_sub_ps(posy8[index - height / 8], posy8[index]);

			sqrLength8 = _mm256_fmadd_ps(rightx8, rightx8, _mm256_mul_ps(righty8, righty8));
			mask = _mm256_cmp_ps(sqrLength8, maxLengthsquared8, _CMP_GT_OQ);

			lengthInv8 = _mm256_invsqrt_ps(sqrLength8);

			halfExcessx8 = _mm256_and_ps(mask, _mm256_mul_ps(aHalf8, _mm256_fnmadd_ps(rightx8, _mm256_mul_ps(lengthInv8, maxLength8), rightx8)));
			halfExcessy8 = _mm256_and_ps(mask, _mm256_mul_ps(aHalf8, _mm256_fnmadd_ps(righty8, _mm256_mul_ps(lengthInv8, maxLength8), righty8)));

			posx8[index] = _mm256_add_ps(posx8[index], halfExcessx8);
			posy8[index] = _mm256_add_ps(posy8[index], halfExcessy8);
			posx8[index - height / 8] = _mm256_sub_ps(posx8[index - height / 8], halfExcessx8);
			posy8[index - height / 8] = _mm256_sub_ps(posy8[index - height / 8], halfExcessy8);

			counter8 = _mm256_or_ps(counter8, mask);
		}
		for (int y = 0; y < height; y++)
		{
			posx[y] = polePos.x;
			posy[y] = polePos.y + y * 1.2f;
		}
		for (int j = 1; j < 8; j++)
		{
			counter[0] = counter[0] | counter[j];
		}
		if (counter[0] == 0)
			break;
	}
	// all done
	uint64_t clocks_elapsed = __rdtsc() - start;
	//cout << clocks_elapsed << endl;
	static uint64_t runningAvg = 1000000;
	static int runs = 1;
	runs++;
	if (runs > 2000)
		return true;
	runningAvg = (995 * runningAvg + 5 * clocks_elapsed)/1000;
	cout << runningAvg << " : " << runs << endl;
	return true; // flags don't die
}

void VerletFlag::Remove()
{
	if (hasBackup) for (int x = width - 1; x >= 0; x--) for (int y = height - 1; y >= 0; y--)
	{
		int index = y + x * height;
		int2 intPos = make_int2( posx[index], posy[index] );
		MyApp::map.bitmap->Plot( intPos.x, intPos.y, backup[index * 4 + 0] );
		MyApp::map.bitmap->Plot( intPos.x + 1, intPos.y, backup[index * 4 + 1] );
		MyApp::map.bitmap->Plot( intPos.x, intPos.y + 1, backup[index * 4 + 2] );
		MyApp::map.bitmap->Plot( intPos.x + 1, intPos.y + 1, backup[index * 4 + 3] );
	}
}