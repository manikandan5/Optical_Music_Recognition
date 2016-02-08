#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <DrawText.h>
#include <cmath>

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedSymbol class may be helpful!
//  Feel free to modify.
//
typedef enum {NOTEHEAD=0, QUARTERREST=1, EIGHTHREST=2} Type;
class DetectedSymbol {
public:
  int row, col, width, height;
  Type type;
  char pitch;
  double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<struct DetectedSymbol> &symbols)
{
  ofstream ofs(filename.c_str());

  for(int i=0; i<symbols.size(); i++)
    {
      const DetectedSymbol &s = symbols[i];
      ofs << s.row << " " << s.col << " " << s.width << " " << s.height << " ";
      if(s.type == NOTEHEAD)
	ofs << "filled_note " << s.pitch;
      else if(s.type == EIGHTHREST)
	ofs << "eighth_rest _";
      else 
	ofs << "quarter_rest _";
      ofs << " " << s.confidence << endl;
    }
}

// Function that outputs a visualization of detected symbols
void  write_detection_image(const string &filename, const vector<DetectedSymbol> &symbols, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];
  for(int i=0; i<3; i++)
    output_planes[i] = input;

  for(int i=0; i<symbols.size(); i++)
    {
      const DetectedSymbol &s = symbols[i];

      overlay_rectangle(output_planes[s.type], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 255, 2);
      overlay_rectangle(output_planes[(s.type+1) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);
      overlay_rectangle(output_planes[(s.type+2) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);

      if(s.type == NOTEHEAD)
	{
	  char str[] = {s.pitch, 0};
	  draw_text(output_planes[0], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[1], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[2], str, s.row, s.col+s.width+1, 0, 2);
	}
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  SDoublePlane output(input.rows(), input.cols());

  // Convolution code here
  
  return output;
}

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
  SDoublePlane output(input.rows(), input.cols());
    
    for (int i = 2; i < input.rows()-2; i++)
    {
        for (int j = 2; j < input.cols()-2; j++)
        {
            
            for (int ki = -2; ki < 3; ki++ )
            {
                for (int kj = -2; kj<3; kj++)
                {
                    output[i][j] = output[i][j] + filter[ki+2][kj+2] * input[i - ki][j - kj];
                    
                }
            }
        }
    }
    return output;
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement a sobel gradient estimation filter with 1-d filters
  

  return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
  
  return output;
}


SDoublePlane inverse(const SDoublePlane &input )
{
  SDoublePlane output(input.rows(),input.cols());
  //Logic for inverse has to be filled
  return output;
}

SDoublePlane flipper(const SDoublePlane &input )
{
  SDoublePlane output(input.rows(),input.cols());
  int i,j;
  for(i=0;i<input.rows();i++)
  {
    for(j=0;j<input.cols();j++)
    {
      output[i][j]= input[input.rows()-i-1][input.cols()-j-1];
    }
  }
  return output;
}

SDoublePlane createFilter()
{
    SDoublePlane gFilter(5,5);
    // initialization of standard deviation to 1.0
    double sigma = 1.0;
    double r, s = 2.0 * sigma * sigma;
    
    // sum is for normalization
    double sum = 0.0;
    
    for (int x = -2; x <= 2; x++)
    {
        for(int y = -2; y <= 2; y++)
        {
            r = sqrt(x*x + y*y);
            gFilter[x + 2][y + 2] = (exp(-(r*r)/s))/(M_PI * s);
            sum += gFilter[x + 2][y + 2];
        }
    }
    
    // normalize the Kernel
    for(int i = 0; i < 5; ++i)
        for(int j = 0; j < 5; ++j)
            gFilter[i][j] /= sum;
    
    return gFilter;
    
}

double maximum(SDoublePlane &input)
{
 double max = 0;
 for(int i=0;i<input.rows();i++)
 {
    for(int j=0;j<input.cols();j++)
    {
        if (max < input[i][j])
        {
            max = input[i][j];
        }
    }
 }
 return max;
}
SDoublePlane binaryImgGen(SDoublePlane input, int threshold)
{
    SDoublePlane output(input.rows(), input.cols());
    
    for (int i = 0; i < input.rows(); i++)
    {
        for (int j = 0; j < input.cols(); j++)
        {
            if(input[i][j]<threshold)
            {
                
                output[i][j] = 0;
            }
            else
            {
                output[i][j] = 1;
            }
            
        }
    }
    return output;
}

vector<DetectedSymbol> symDetectionByTemplate(SDoublePlane &input_image,SDoublePlane &template_1)
{
  vector<DetectedSymbol> symbols;
  SDoublePlane gFilter(5,5);
  gFilter = createFilter();
  double thresh = 100;

  SDoublePlane output_image = convolve_general(input_image, gFilter);

  SImageIO::write_png_file("Blurred.png", output_image, output_image, output_image);

  SDoublePlane output_image1 = binaryImgGen(output_image, 200);

  SImageIO::write_png_file("Threshold.png", output_image1, output_image1, output_image1);


  //SDoublePlane output_image = flipper(mean_filter); for testing the flipper function
  //compute the distance function
  SDoublePlane binary_template_1 = binaryImgGen(template_1, thresh);
  SDoublePlane binary_input_image = binaryImgGen(input_image, thresh);
  SDoublePlane inverse_template_1 = inverse(binary_template_1);
  SDoublePlane flipped_inverse_template_1 = flipper(inverse_template_1);
  SDoublePlane flipped_template_1 = flipper(binary_template_1);
  SDoublePlane inverse_input_image = inverse(binary_input_image);
  SDoublePlane F = convolve_general(input_image,flipped_template_1) + convolve_general(inverse_input_image,flipped_inverse_template_1) ;

  //find the maximum value in F
  int max ;
  max=maximum(F);
  //find indexes of maxima in F
  return symbols;
}

class Dimensions
{
 public:
  Dimensions()
  {
    row_coordinate = 0;
    spacing = 0;
  }
  int getCoordinate()
  {
    return row_coordinate;
  }

  int getSpacing()
  {
    return spacing;
  }
 protected:
  int row_coordinate;
  int spacing;
};

Dimensions findStaff(SDoublePlane &input)
{
  int inp_row = input.rows();
  int inp_col = input.cols();
  SDoublePlane output(inp_row, inp_col);
  int accum_r = inp_row;
  int accum_d = inp_row/10;
  _DTwoDimArray<double> accum(accum_r, accum_d);
  for(int i=0;i<(inp_row-(4*accum_d));i++)
  {
    for(int j=0;j<inp_col;j++)
    {
      for(int h=1;h<=accum_d;h++)
      {
        if (input[i][j] == 1 && input[i+h][j] == 1 && input[i+(2*h)][j] == 1 && input[i+(3*h)][j] == 1 && input[i+(4*h)][j] == 1 )
        {
            accum[i][h]++;
        }
      }
    }
  }



  return ;
}
//
// This main file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//
int main(int argc, char *argv[])
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
  
  // test step 2 by applying gaussian filters to the input image
  SDoublePlane mean_filter(3,3);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      mean_filter[i][j] = 1/9.0;

  SDoublePlane gFilter(5,5);
  gFilter = createFilter();


   
  // randomly generate some detected symbols -- you'll want to replace this
  //  with your symbol detection code obviously!
  vector<DetectedSymbol> symbols;
  for(int i=0; i<10; i++)
    {
      DetectedSymbol s;
      s.row = rand() % input_image.rows();
      s.col = rand() % input_image.cols();
      s.width = 20;
      s.height = 20;
      s.type = (Type) (rand() % 3);
      s.confidence = rand();
      s.pitch = (rand() % 7) + 'A';
      symbols.push_back(s);
    }

  write_detection_txt("detected.txt", symbols);
  write_detection_image("detected.png", symbols, input_image);
}
