#include <iostream>
#include <stdio.h>

using namespace std;

// readinputs: This quick program reads the inputs.dat file and spits out what's in there, 
// just so I can check it without re-running the setup in matlab.

// this 2D version is the same as the 3D version, except it does not have Bvec.

// ---------------------------------------------

int main()
{

  FILE * inputFile;
  inputFile = fopen("inputs.dat","rb");

  int dopml_top;
  int dopml_wall;
  int doioniz;
  int doelve;
  int dodetach;
  float maxalt;
  float stepalt;
  float dr1;
  float dr2;
  float range;
  float dt;
  int tsteps;
  float sig;
  float sigm;
  float Bmag;
  float Bvec [3];
  float camdist;
  float camalt;
  float L;
  int numfiles;

  fread(&dopml_top,sizeof(int),1,inputFile);
  fread(&dopml_wall,sizeof(int),1,inputFile);
  fread(&doioniz,sizeof(int),1,inputFile);
  fread(&doelve,sizeof(int),1,inputFile);
  fread(&dodetach,sizeof(int),1,inputFile);
  fread(&maxalt,sizeof(float),1,inputFile);
  fread(&stepalt,sizeof(float),1,inputFile);
  fread(&dr1,sizeof(float),1,inputFile);
  fread(&dr2,sizeof(float),1,inputFile);
  fread(&range,sizeof(float),1,inputFile);  
  fread(&dt,sizeof(float),1,inputFile);
  fread(&tsteps,sizeof(int),1,inputFile);
  fread(&sig,sizeof(float),1,inputFile);
  fread(&sigm,sizeof(float),1,inputFile);
  fread(&Bmag,sizeof(float),1,inputFile);
  fread(&Bvec,sizeof(float),3,inputFile);
  fread(&camdist,sizeof(float),1,inputFile);
  fread(&camalt,sizeof(float),1,inputFile);
  
  float Isource [tsteps];
  fread(&Isource,sizeof(float),tsteps,inputFile);
  fread(&L,sizeof(float),1,inputFile);
  fread(&numfiles,sizeof(int),1,inputFile);
  fclose(inputFile);

  // find the maximum of the current pulse, and find the time between where it crosses 10% on either side

  float maxcurrent = 0;
  for (int i = 0; i < tsteps; i++) {
    if (Isource[i] > maxcurrent) {
	maxcurrent = Isource[i];
      }
  }

  int riseflag = 0;
  int fallflag = 0;

  for (int i = 0; i < tsteps; i++) {
    if (Isource[i] > maxcurrent/10) {
      riseflag = i;
      break;
    }
  }
  for (int i = riseflag+1; i < tsteps; i++) {
    if (Isource[i] < maxcurrent/10) {
      fallflag = i;
      break;
    }
  }
  float pulsetime = (fallflag - riseflag)*dt;

  // okay now spit it all out to terminal:

  cout << "-----------------------------------------------\n";
  cout << "Inputs for 2D EMP Code:\n";
  cout << "-----------------------------------------------\n";
  cout << "Do PML = " << dopml_top << " (top), " << dopml_wall << " (wall)\n";
  cout << "Do Ionization = " << doioniz << "\n";
  cout << "Do Elve = " << doelve << "\n";
  cout << "Do detachment process = " << dodetach << "\n\n";
  cout << "Maximum Altitude = " << maxalt/1e3 << " km\n";
  cout << "Step Altitude = " << stepalt/1e3 << " km\n";
  cout << "dr1 = " << dr1/1e3 << " km\n";
  cout << "dr2 = " << dr2/1e3 << " km\n";
  cout << "range = " << range/1e3 << " km\n\n";
  cout << "dt = " << dt*1e6 << " us\n";
  cout << "time steps = " << tsteps << "\n\n";
  cout << "sigma = " << sig << "\n";
  cout << "sigma (magnetic) = " << sigm << "\n";
  cout << "B-field magnitude = " << Bmag*1e9 << " nT, direction: " << Bvec[0] << ", " << Bvec[1] << ", " << Bvec[2] << "\n\n";
  cout << "Peak current = " << maxcurrent/1e3 << " kA\n";
  cout << "Current pulse time = " << pulsetime*1e6 << " us\n";
  cout << "Altitude of current source = " << L/1e3 << " km\n";
  cout << "Camera located at " << camdist/1e3 << " km distance, " << camalt/1e3 << " km altitude\n";
  cout << "Will print " << numfiles << " output files\n";
  cout << "-----------------------------------------------\n";

  return 0;
}
