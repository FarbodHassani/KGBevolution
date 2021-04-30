#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "chealpix.h"

using namespace std;

int main(int argc, char **argv)
{
	FILE * infile = NULL;
	float * map = NULL;
	//float * reversemap = NULL;
	int Nside = 0;
	char coordsys = 'G';
	
	if (argc < 4) return -1;
	
	Nside = atoi(argv[3]);
	
	cout << " Nside = " << Nside << endl;
	
	map = (float *) malloc(sizeof(float) * 12 * Nside * Nside);
	//reversemap = (float *) malloc(sizeof(float) * 12 * Nside * Nside);
	
	infile = fopen(argv[1], "rb");
	
	if (infile == NULL)
	{
		cout << " error opening file " << argv[1] << " for reading!" << endl;
		return -1;
	}
	
	if (fread((void *) map, sizeof(float), 12 * Nside * Nside, infile) != 12 * Nside * Nside)
	{
		cout << " error reading data from " << argv[1] << "!" << endl;
		return -1;
	}
	
	fclose(infile);
	
	/*cout << " reversing map... " << endl;
	
	for (int ring = 1; ring <= Nside; ring++)
	{
		for (long i = 0; i < 4*ring; i++)
			reversemap[2*ring*(ring-1)+i] = map[12*Nside*Nside-2*ring*(ring+1)+i];
	}
	
	for (int ring = Nside+1; ring < 3*Nside; ring++)
	{
		for (long i = 0; i < 4*Nside; i++)
			reversemap[2*Nside*(Nside+1)+(ring-(Nside+1))*4*Nside+i] = map[12*Nside*Nside-2*Nside*(Nside+1)-(ring-Nside)*4*Nside+i];
	}
	
	for (int ring = Nside; ring >= 1; ring--)
	{
		for (long i = 0; i < 4*ring; i++)
			reversemap[12*Nside*Nside-2*ring*(ring+1)+i] = map[2*ring*(ring-1)+i];
	}*/
	
	cout << " writing data to " << argv[2] << endl;
	
	write_healpix_map(map, Nside, argv[2], 0, &coordsys);
	
	free(map);
	//free(reversemap);
	
	return 0;
}

