#ifndef __SIMAGE_H__
#define __SIMAGE_H__

#include <DTwoDimArray.h>
#include <string.h>


// A very simple image class.
//
class SDoublePlane : public _DTwoDimArray<double>
{
 public:
  SDoublePlane() { }
  SDoublePlane(int _rows, int _cols)  : _DTwoDimArray<double>(_rows, _cols)
    { 
      // be nice and initialize plane to all 0's
      memset(data_ptr(), 0, sizeof(double) * rows() * cols());
    }

  ///////////////////////////////////////////////////////
  // Addition operator

  SDoublePlane &operator+(const SDoublePlane &a)
    {
      // profiler->begin(4);
      if( _rows == a.rows() && _cols == a.cols())
	{
	  for(int i=0;i<_rows;i++)
	  {
	    for(int j=0;j<_cols;j++)
	    {
	        data[i][j] += a[i][j];
	    }
	  }
	}

      // profiler->end(4);
      return *this;
    }


};

#endif
