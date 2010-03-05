

#include "MCMC.h"

using namespace std;


void ArrayForm2D(std::string ContourSlice, uint32_t contour_length, vtkPoints ** pPoints)
{
	int height = 512;
	int width = 512;	
	///////////////////  Get the array of points, load into FloatArray
	double tuple[3];
	int tuple_index = 0;
	// value_char/value_float : current float in char or float format
	string coord_string ;
	double coord = 0;

	double pixel_spacing = 0.70313;
	// loop all the chars, load into tuple, then into contour	
	for (int i=0; i<contour_length; i++ ){
			if (ContourSlice[i]!= '\\') 
				coord_string.push_back(ContourSlice[i]);
			else{
				coord = atof(coord_string.c_str());
				tuple[tuple_index] = coord;
				tuple_index ++;
				if (tuple_index==3)
				{
					// cor2pixel
					tuple[0] = tuple[0]/pixel_spacing + height;
					tuple[1] = tuple[1]/pixel_spacing + width;

					// insert current tuple
					tuple_index = 0;					
					(*pPoints)->InsertNextPoint( tuple );
					
				}
				coord_string.clear();
			}
	}

}

