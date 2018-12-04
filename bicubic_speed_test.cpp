//compile with g++ bicubic_speed_test.cpp -o place_bets -lgsl -lm

#include "vegas_tables.h"

using namespace std;

int main()
{
  unsigned int Ntests;
  double xi,T,T_min,T_max,xi_min,xi_max;
  double rate;
  double eps = 1e-5;

  T_min  = 1e4; 
  T_max  = 1e7; 
  xi_min = 1e0; 
  xi_max = 1e3; 
  //const double T_bbox[2] = {10.*T_min,0.1*T_max};     // optional: won't use the full table
  //const double xi_bbox[2] = {10.*xi_min,0.9*xi_max};  // may save on memory
  
  /* A couple things to try
  - setting Ntests = 1000 is useful for testing speed
  - setting Ntests = 10 is useful for sampling random
    values to judge %-errors in interpolation schemes
  */

  Ntests = 1000;

  // Instantiate two different table instances 
  VEGAS_LUT hc_xstar_bicubic("Blondin_100by200.dat","bicubic");		
  VEGAS_LUT hc_xstar_bilinear("Blondin_100by200.dat","bilinear");
  hc_xstar_bicubic.initialize_table();
  hc_xstar_bilinear.initialize_table();
  //hc_xstar_bicubic.initialize_table(T_bbox,xi_bbox);  // optional
  

  double ttable[Ntests];
  double xtable[Ntests];
  
  for(int i=0; i < Ntests; i++)
  {
    ttable[i] = rand()%(int)(T_max-T_min)+T_min + eps;
    xtable[i] = rand()%(int)(xi_max-xi_min)+xi_min + eps;
  }
  

  clock_t start1 = clock();
  for(int i=0; i < Ntests; i++)
  {
    T = ttable[i];
    xi = xtable[i];
   
   cout << setprecision(9) 
    << "Bicubic calcs: T = " << T << "  xi = " << xi;

    rate = hc_xstar_bicubic.get_rate(T,xi);
    
    cout << setprecision(9) 
    << "  VEGAS_LUT rate = " << rate << endl;
  }
  cout << endl;

  clock_t start2 = clock();
  for(int i=0; i < Ntests; i++)
  {
    T = ttable[i];
    xi = xtable[i];

    cout << setprecision(9)
    << "Bilinear calcs: T = " << T << "  xi = " << xi;

	  rate = hc_xstar_bilinear.get_rate(T,xi);
	
    cout << setprecision(9)
    << "  VEGAD_LUT rate = " << rate << endl;
  }
  clock_t end = clock();
  
  // compare times
  double time1 = (start2 - start1) / (double)CLOCKS_PER_SEC;
  double time2 = (end - start2) / (double)CLOCKS_PER_SEC;
  cout << endl;
  cout<<"Bilinear took " << time1 << " seconds." << endl;
  cout<<"Bicubic took " << time2 << " seconds." << endl;
  cout<<"Bilinear is " << setprecision(3) << time1/time2
  << " times faster than bicubic." << endl;

  return 0;
}