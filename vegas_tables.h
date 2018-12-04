/*
VEGAS_TABLE Version 1.  5/27/16
report issues to Tim Waters,
waters@lanl.gov

**** Basic code description: 
Input: 	An ascii table with entries
		on a uniform grid of (x,y) values
		and a (x,y) value for which an
    interpolated value is sought
Output: The interpolated value 
     corresponding to that (x,y) pair

**** Detailed description:  

>> Step 1. Instantiation:
The class constructor calls a method
to parse an ascii text file which contains 
an (N+1)x(M+1) table:
The first row contains the M values of x.
The first column contains the N values of y.
The remainder of the table contains the NxM 
values of z.  

Here is a simple example of an input file:

4	1e0     1e1     1e2     1e3
1e4	1e-20	2e-20	3e-20	4e-20
1e5	5e-20	6e-20	7e-20	8e-20
1e6	1e-19	2e-19	3e-19	4e-19
1e7	5e-19	6e-19	7e-19	8e-19
1e8	1e-18	2e-18	3e-18	4e-18

In this example, N = 5 and M = 4.  
M must be specified as the first element.
The remaining elements in row 1 are values
of the photoionization parameter, xi.
The left column below M are the temperature 
values.  The remaining elements are net heating 
rates that are functions of T and xi.  

By default, the input file is limited to size 
1000 x 1000.  To increase this size limit, 
increase Nmax in the default constructor.

>> Step 2. Table initialization:
The file parsed above is read into a 
table, but this is not the actual table 
used.  The actual table is made upon calling
an initialization function, and optionally
array arguments can be passed specifying the 
min/max of both T and xi. This may be useful
when combining multiple overlapping tables
and will allow for flexibility if tables 
generated from non-uniform grids are needed,
although a table-search algorithm has not 
been implemented.  The input table is then 
reduced to a size that contains only the 
values between these limits.  

>> Step 3. Invocation:
The code returns the interpolated value
of whatever rate was in the input file. 
The code uses the bicubic and bilinear
interpolation algorithms from the GSL
library.  By default, bicubic is used.

Importantly, the code assumes that a 
nearly uniform grid was produced via 
x = j*(xmax-xmin)/(M-1) + xmin 
y = i*(ymax-ymin)/(N-1) + ymin
This allows avoiding an expensive
table search for the lookup operation.
For a uniform grid generated from the
above formulas, the table indices 
for a desired value of (x,y) are
i_x = (M-1)*(x - xmin)/(xmax - xmin)
i_y = (N-1)*(y - ymin)/(ymax - ymin)
Since xstar may slightly modify the
grid locations (x,y) originally 
specified, our lookup operation
consists of first guessing that the
lookup value indices are (i_x,i_y)
and then allowing for some wiggle 
room by searching 1 index left, right, 
above, and below upon finding that our 
guess was incorrect.  The code will 
terminate if this search fails.  

**** Example usages: 
1. Bicubic interpolation on a table that was
generated from a log-log spaced grid
VEGAS_LUT hc_xstar("input_table.tab","bicubic");			
hc_xstar.initialize_table();
hc_rate = hc_xstar.get_rate(T,xi);

1a. If lin-lin spacing was used instead, use
VEGAS_LUT hc_xstar("input_table.tab","bicubic",false,false);

1c. To apply a bounding box, use instead
const double T_bbox[2] = {1e5,1e7};
const double xi_bbox[2] = {1e1,5e2};
hc_xstar.initialize_table(T_bbox,xi_bbox);
*/

#ifndef H_VEGAS_TABLE
#define H_VEGAS_TABLE

#include <iostream>		// std::cout, std::endl
#include <iomanip>      // std::setprecision
#include <fstream>		// std::ifstream
#include <string>
#include <vector>
#include <time.h> 
#include <cmath>

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>


class VEGAS_LUT
{
public:
    //default constructor will parse the input file
    //log-log spacing for (T,xi) is assumed by default
	VEGAS_LUT(	const std::string a="VEGAS_table.tab", 
				const std::string b="bicubic", 
				bool c=true, bool d=true): 
	filename(a),gsl_routine(b),logT(c),logxi(d) 
	{
	 // allocate a (Nmax x Nmax) table for reading in the file
	  int Nmax = 1000;
	  input_table = new double*[Nmax];
	  for (int i = 0; i < Nmax; i++)
    	input_table[i] = new double[Nmax];
      xivals = new double[Nmax];
      Tvals = new double[Nmax];
    
      gsl_init();
      parse_file();
	}
	~VEGAS_LUT()  { free_memory(); } 
	
	//user must call the initialize_table method
	void initialize_table(const double*, const double*); 
	
	//user function that returns an interpolated value
	double get_rate(double, double);
	
	// actual table (T,xi) bounds: 
	// hydro code needs to know Temp bounds when doing implicit root finding
	double y_min,y_max,x_min,x_max; 	

private:
	bool logT,logxi;					// determines log gridding
	double* xlim;						// input table xi bounds
	double* ylim;						// input table T bounds
	const double* T_bounds;				// user specified T bounds
	const double* xi_bounds;			// user specified xi bounds
	double** 	input_table; 			// pointer to input table
	double** 	table;					// pointer to actual table
	double* 	xivals; 				// xi-vals from input table 
	double*		Tvals; 					// T-vals from input table 
	double* 	xvals; 					// xi-vals for actual table
	double*		yvals; 					// T-vals for actual table
	unsigned int 	NT,	  // # rows of input table (# of Ts)
					Nxi,  // # cols of input table (# of xis)
					N,	  // # rows of actual table (# of Ts)
					M,	  // # cols of actual table (# of xis)
					i_T,  // T-offset between the input and actual tables
					i_xi; // xi-offset between the input and actual tables
	typedef std::vector<double> array_t;	//dynamic array data type
	std::string filename, gsl_routine;		//input strings
	
	//member functions for handling the input file
	void parse_file();
	bool read_value(std::ifstream&, const int, array_t&);
	void free_input_table(); 
	void free_memory(); 
	
	//other member functions
	void reduce_input_table();
	void make_table(unsigned int, unsigned int);
	void determine_box_indices(double, double, unsigned int&, unsigned int&);
	
	//GSL variables and member functions
	bool use_gsl_bicubic;
	const gsl_interp2d_type *gsl_scheme;
	typedef gsl_spline2d* spline_t;
	spline_t **spline_table;
    gsl_interp_accel *xacc, *yacc;
    size_t nx,ny;
    void gsl_init();
    void make_spline_table();
};


void VEGAS_LUT::gsl_init()
{
    if (gsl_routine == "bilinear")
    {
      gsl_scheme = gsl_interp2d_bilinear;
      use_gsl_bicubic = false;
      nx = 2;
      ny = 2;
    }
    else if (gsl_routine == "bicubic")
    {
      gsl_scheme = gsl_interp2d_bicubic;
      use_gsl_bicubic = true;
      nx = 4;
      ny = 4;
    }
    else
    {
      std::cout << "ERROR: GSL interpolation scheme " << gsl_routine 
      << " not recognized. Choose from:\n"
      << "bilinear\n"
      << "bicubic\n"
      << std::endl;
      exit(0);
    }
}

bool VEGAS_LUT::read_value(std::ifstream& fin, const int j, array_t& input_data)
{
  double value;
  if (fin >> value) 
  { 
    input_data[j]=value; 
    return true; 
  }
  else return false;
}

void VEGAS_LUT::parse_file()
{
  std::ifstream fin(filename.c_str());
  if (!fin)
  {
    std::cout << "File " << filename << " does not exist!" << std::endl;
	exit(0);
  }	
  bool file_is_open = false;
  bool end_of_line = false;
  
  // read in first value: Nxi
  array_t val1(1);
  file_is_open = read_value(fin,0,val1);
  Nxi = (unsigned int) val1[0];
  
  // read in remainder of first row: all xi-values
  array_t this_row(Nxi+1);
  int i=1,j;
  while(!end_of_line)
  {
    file_is_open = read_value(fin,i,this_row);
    xivals[i-1] = this_row[i];
    if( i!=0 && i%Nxi == 0 )
        end_of_line = true;
    else i++; //next row element
  }
    
  // read in 1st element of 2nd row: the 1st T value
  file_is_open = read_value(fin,0,this_row);
  Tvals[0] = this_row[0];
  
  // read in the rest of the file
  end_of_line = false;
  i = 1; j=0;
  while(file_is_open)
  {
    file_is_open = read_value(fin,i,this_row);

    if( i!=0 && i%Nxi == 0 )
    {
        end_of_line = true;
        i=0; // j++; //reset column, increment row
    }
    else i++; //next row element

    if(end_of_line)
    {
      // 1st element of this_row is T
      Tvals[j] = this_row[0];

      // remaining elements are rates
      for (int ii=0; ii<Nxi; ii++)
        //std::cout << "this_row[" << ii << "] = " << this_row[ii] << std::endl;
        input_table[j][ii] = this_row[ii+1];
      
      j++;
    }
    end_of_line = false;
  }

  NT = j; //number of temperature values
  fin.close();
  
  if (Nxi < 2 || NT < 2)
    {
    std::cout << "ERROR: Table dimensions must be at least 2x2!" << std::endl;
    exit(0);
    }
}

void VEGAS_LUT::make_table(unsigned int i0, unsigned int j0)
{
    
    // allocate memory
	table = new double*[N];
	for (int j = 0; j < N; j++)
    	table[j] = new double[M];
    
    // assign values
	for (int j = 0; j < N; j++) 
	  for (int i = 0; i < M; i++)
	    table[j][i] = input_table[j0+j][i0+i];
}

void VEGAS_LUT::make_spline_table() 
{
    
    // allocate memory for spline_table
	spline_table = new spline_t*[N];
	for (int j = 0; j < N; j++)
    	spline_table[j] = new spline_t[M];
    
    /* Each location in our table is to be thought of as a 4-corner box.
     * The corners have locations (xb[], yb[]) and function values (zb[])
     * and this data is collected ahead of time so that GSL only needs to 
     * access it to carry out the interpolation.  For the case of bicubic 
     * interpolation, GSL will also take numerical derivatives and store 
     * all that data, so a lot of expense is saved by collecting this at
     * initialization. */
    double *xb = new double[nx];
	double *yb = new double[ny];
	double *zb = new double[nx*ny];
		
	/* Associations: 
	T  = y(j), with j from (0,N-1)
	xi = x(i), with i from (0,M-1)     
	*/
	for (int j=1; j < (N-2); j++)
      for (int i=1; i < (M-2); i++)
      {
      	// allocate storage for interpolation data 
        gsl_spline2d *spline = gsl_spline2d_alloc(gsl_scheme, nx, ny);
	    
        // evaluate position data (xb[],yb[]) and store function values (zb[]) to spline
	    for (int jb=0; jb < nx; jb++)
		  for (int ib=0; ib < ny; ib++)
		  {
		    if (use_gsl_bicubic)
		    {
		      xb[ib] = xvals[i+ib-1];
		      yb[jb] = yvals[j+jb-1];
			  gsl_spline2d_set(spline, zb, ib, jb, table[j+jb-1][i+ib-1]);
		    }
		    else 
		    {
		      xb[ib] = xvals[i+ib];
		      yb[jb] = yvals[j+jb];
		      gsl_spline2d_set(spline, zb, ib, jb, table[j+jb][i+ib]);
		    }
		  }
        
          /* add position data 
          * In the case of bicubic, this also internally calculates 
          * derivative data and stores that as well. */
          gsl_spline2d_init(spline, xb, yb, zb, nx, ny);
        
          // store this interpolation data in our spline table
          spline_table[j][i] = spline;
      } //end loop
    
    /* The above loop excludes all boxes on table edges, which always require
     * doing bilinear interpolation, so add that data now */
    int j_edges[2] = {0, N-2};
    int i_edges[2] = {0, M-2};
    int i,j;
    
    for (int jj=0; jj < 2; jj++)
      for (int i=0; i < (M-1); i++)
      {
        j = j_edges[jj];
      
		gsl_spline2d *spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, 2, 2);
		
		// evaluate position data (xb[],yb[]) and store function values (zb[]) to spline
		  for (int jb=0; jb < 2; jb++)
			for (int ib=0; ib < 2; ib++)
			{
				xb[ib] = xvals[i+ib];
				yb[jb] = yvals[j+jb];
				gsl_spline2d_set(spline, zb, ib, jb, table[j+jb][i+ib]);
			}
		
		/* add position data 
        * In the case of bicubic, this also internally calculates 
        * derivative data and stores that as well. */
        gsl_spline2d_init(spline, xb, yb, zb, 2, 2);
        
        // store this interpolation data in our spline table
        //std::cout<<  "(i,j) = (" << i << ", " << j << ")" << std::endl;
        spline_table[j][i] = spline;
      }	
      
    for (int j=0; j < (N-1); j++)
      for (int ii=0; ii < 2; ii++)
      {
        i = i_edges[ii];
      
		gsl_spline2d *spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, 2, 2);
		
		// evaluate position data (xb[],yb[]) and store function values (zb[]) to spline
		  for (int jb=0; jb < 2; jb++)
			for (int ib=0; ib < 2; ib++)
			{
				xb[ib] = xvals[i+ib];
				yb[jb] = yvals[j+jb];
				gsl_spline2d_set(spline, zb, ib, jb, table[j+jb][i+ib]);
			}
		
		/* add position data 
        * In the case of bicubic, this also internally calculates 
        * derivative data and stores that as well. */
        gsl_spline2d_init(spline, xb, yb, zb, 2, 2);
        
        // store this interpolation data in our spline table
        //std::cout<<  "(i,j) = (" << i << ", " << j << ")" << std::endl;
        spline_table[j][i] = spline;
      }	
    
    // These I think are used when a lookup table must be searched
    // but since we know the lookup locations there is no point
    xacc = NULL; //gsl_interp_accel_alloc();
    yacc = NULL; //gsl_interp_accel_alloc(); 

}


void VEGAS_LUT::free_memory()
{
  for (int i = 0; i < N; i++)
  {
    delete[] table[i];
    delete[] spline_table[i];
  }
  
  //gsl_interp_accel_free(xacc); 
  //gsl_interp_accel_free(yacc);
  
}


void VEGAS_LUT::free_input_table()
{
  for (int i = 0; i < N; i++)
    delete[] input_table[i];
}

void VEGAS_LUT::reduce_input_table()
{
    
    double Tmin,Tmax,ximin,ximax;
    Tmin = T_bounds[0];
    Tmax = T_bounds[1];
    ximin = xi_bounds[0];
    ximax = xi_bounds[1];
    
    if (Tmin < Tvals[0] || Tmax > Tvals[NT-1])
    {
    std::cout << "ERROR: Tmin or Tmax is not contained within file " 
    << filename << "!"
    << std::endl;
    exit(0);
    }
     
    if (ximin < xivals[0] || ximax > xivals[Nxi-1])
    {
    std::cout << "ERROR: x_min or x_max is not contained within file " 
    << filename << "!"
    << std::endl;
    exit(0);
    }
      
    // find the indices of the table corresponding to input min/max values
    int i = 0;
    while (Tvals[i] < Tmax && i < NT) 
      i++;
    y_max = Tvals[i];  //so y_max is slightly greater than Tmax
    
    int ii = 0;
    while (Tvals[ii] <= Tmin && ii < i) 
      ii++;
    y_min = Tvals[ii-1]; //so y_min is slightly less than Tmin
    
    // assign private variables
    i_T = ii-1; // start position in input table
    N = i - ii + 2; // # of T-vals in desired table
    
    i = 0;
    while (xivals[i] < ximax && i < Nxi) 
      i++;
    x_max = xivals[i]; //so x_max is slightly greater than ximax
    
    ii = 0;
    while (xivals[ii] <= ximin && ii < i) 
      ii++;
    x_min = xivals[ii-1]; //so x_min is slightly less than ximin
    
    // assign private variables
    i_xi = ii-1; // start position in input table
    M = i - ii + 2; // # of xi-vals in desired table
}

void VEGAS_LUT::initialize_table(const double * T_bbox = NULL, const double * xi_bbox = NULL)
{	
	int i0,j0; 
	
	if (T_bbox != NULL || xi_bbox != NULL)
	{
	  T_bounds = T_bbox;
	  xi_bounds = xi_bbox;
	  reduce_input_table(); // this sets i_xi, i_T, N, M, y_min, Tmax, x_min, x_max
	  i0 = i_xi; 
	  j0 = i_T;
	}
	else // the actual table will be the same size as the input table
	{
	  y_min = Tvals[0];
	  y_max = Tvals[NT-1];
	  x_min = xivals[0];
	  x_max = xivals[Nxi-1];
	  N = NT;
	  M = Nxi;
	  i0 = 0;
	  j0 = 0;
	}
	
	// store the values of xi and T that will be used
    // they way they will be used (i.e. log or lin)
    xvals = new double[M];
    yvals = new double[N];
    double x,y;
    for (int i=0; i<M; i++) 
    {
      if(logxi)	x = log10(xivals[i0+i]); 
      else 		x = xivals[i0+i];
      xvals[i] = x;
    }
    for (int j=0; j<N; j++) 
    {
      if (logT)	y = log10(Tvals[j0+j]); 
      else 		y = Tvals[j0+j];
      yvals[j] = y;
    }
	
    // generate the actual table of H/C rates
	make_table(i0,j0);
	
	// alloc storage for max/min vals of T,xi
	xlim = new double[2]; // xi limits
  ylim = new double[2]; // T limits
    
    // set any log gridding: use limits of input table
	 if (logT) {
	   ylim[0] = log10(Tvals[0]);  
	   ylim[1] = log10(Tvals[NT-1]); }
	  else {
	   ylim[0] = Tvals[0];  
	   ylim[1] = Tvals[NT-1]; }
     
     if (logxi) {
       xlim[0] = log10(xivals[0]);
       xlim[1] = log10(xivals[Nxi-1]); }
     else {
       xlim[0] = xivals[0];
       xlim[1] = xivals[Nxi-1]; }

    // generate the spline table for GSL
	make_spline_table();

	// delete input_table
	free_input_table();
	
}

void VEGAS_LUT::determine_box_indices(double x, double y, unsigned int &i_x, unsigned int &i_y)
{
  unsigned int i,j;
  bool correct_guess = false;
  
  /* First guess the indices assuming the input table was
   * generated from a uniform grid: the actual xi and T values
   * that xstar produces data for may be slightly different 
   * as xstar's stopping criteria algorithms may slightly adjust 
   * the input values.  On rare occasions this guess can be off 
   * by 1 index, so here we allow for that. */ 
  i = (Nxi-1) * (x - xlim[0]) / (xlim[1] - xlim[0]) - i_xi;
  j = (NT-1) * (y - ylim[0]) / (ylim[1] - ylim[0]) - i_T ; 

  /* now we check our guess in x and if needed, adjust one box left/right */
  if (x > xvals[i] && x < xvals[i+1])  //then we guessed right
  {
    i_x = i;
    correct_guess = true;
  }
  
  if (!correct_guess)
  {
    if (i==0 || i==(M-1))
    {
      if (i==0 && x < xvals[2])
        i_x = 1;
      else if (i==(M-1) && x > xvals[M-2])
        i_x = M-2;
    }
    else if (x < xvals[i]) //then move to the left 1
      i_x = i-1;
    else if (x > xvals[i+1]) //then move to the right 1
      i_x = i+1;
    else // we give up
    {
      std::cout << "ERROR: search for the table xi-index failed!" 
      			<< "\nDouble check inputs or try debugging "
      			<< "function determine_box_indices()"
      			<< std::endl;
      exit(0);
    }
  }
  
  /* now we check our guess in y and if needed, adjust one box up/down */
  correct_guess = false; // reset 
  
  if (y > yvals[j] && y < yvals[j+1])  //then we guessed right
  {
    i_y = j;
    correct_guess = true;
  }
  
  if (!correct_guess)
  {
    if (j==0 || j==(N-1))
    {
      if (j==0 && y < yvals[2])
        i_y = 1;
      else if (j==(N-1) && y > yvals[N-2])
        i_y = N-2;
    }
    else if (y < yvals[j]) //then move down 1
      i_y = j-1;
    else if (y > yvals[j+1]) //then move up 1
      i_y = j+1;
    else // we give up
    {
      std::cout << "ERROR: search for the table T-index failed!" 
      			<< "\nDouble check inputs or try debugging "
      			<< "function determine_box_indices()"
      			<< std::endl;
      exit(0);
    }
  }
  
}

double VEGAS_LUT::get_rate(double y_val, double x_val)
{   
	unsigned int i_x,i_y;
	double x,y;

	/* terminate program if values exceed bounding box 
    if (y_val < y_min || y_val > y_max || x_val < x_min || x_val > x_max)
    {
      std::cout 
      << "\nFATAL ERROR: (y_val,x_val) = (" << y_val << ", " << x_val << ") " 
      << "lies outside of table's bounding box!\n" 
      << ">> Bounding Box:\n" 
      << "(y_min,y_max) = (" << y_min << ", " << y_max << ")\n" 
      << "(x_min,x_max) = (" << x_min << ", " << x_max << ")\n"
      << std::endl;
	  exit(0);
    } 
    */
    
    /* use boundary rates for values outside bounding box */
    double eps = 1e-10;
    if (y_val < y_min) y_val = y_min*(1. + eps);
    if (y_val > y_max) y_val = y_max*(1. - eps);
    if (x_val < x_min) x_val = x_min*(1. + eps);
    if (x_val > x_max) x_val = x_max*(1. - eps);
	
	// apply any log scaling 
	if (logT)  y = log10(y_val);
	else       y = y_val;
  if (logxi) x = log10(x_val);
  else       x = x_val;
    
	 determine_box_indices(x,y,i_x,i_y);
  
	return gsl_spline2d_eval(spline_table[i_y][i_x], x, y, xacc, yacc);   
}

#endif //H_VEGAS_TABLE
