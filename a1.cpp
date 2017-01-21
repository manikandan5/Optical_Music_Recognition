#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <DrawText.h>
#include <cmath>
#include <limits>
#include <ctime>

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
    for(int w=-width/2; w<=width/2; w++)
    {
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

void overlay_line(SDoublePlane &input, int _top, int _space, double graylevel, int width)
{
    for(int w=-width/2; w<=width/2; w++)
    {
        int top = _top+w, space=_space+w;
        
        // if any of the coordinates are out-of-bounds, truncate them
        top = min( max( top, 0 ), input.rows()-1);
        space = min( max( space, 0 ), input.rows()-1);
        
        // draw all the Staff lines
        for(int j=0; j<=input.cols(); j++)
            input[top][j] = input[top+space][j] = input[top+(2*space)][j] = input[top+(3*space)][j] = input[top+(4*space)][j] = graylevel;
    }
}
// DetectedSymbol class may be helpful!
//  Feel free to modify.
//
typedef enum {NOTEHEAD=0, QUARTERREST=1, EIGHTHREST=2} Type;
class DetectedSymbol
{
public:
    int row, col, width, height;
    Type type;
    char pitch;
    double confidence;
};

class Dimensions
{
public:
    int row_coordinate;
    int spacing;
    bool trebble;
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

void  write_staff_detection_image(const string &filename, const vector<Dimensions> &symbols, const SDoublePlane &input)
{
    SDoublePlane output_planes[3];
    
    for(int i=0; i<3; i++)
        output_planes[i] = input;
    
    for(int i=0; i<symbols.size(); i++)
    {
        const Dimensions &s = symbols[i];
        
        overlay_line(output_planes[0], s.row_coordinate, s.spacing ,0, 1);
        overlay_line(output_planes[1], s.row_coordinate, s.spacing ,0, 1);
        overlay_line(output_planes[2], s.row_coordinate, s.spacing ,255, 1);
    }
    
    SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}

SDoublePlane flipper(const SDoublePlane &input )
{
    SDoublePlane output(input.rows(),input.cols());
    
    for(int i=0;i<input.rows();i++)
    {
        for(int j=0;j<input.cols();j++)
        {
            output[i][j]= input[input.rows()-i-1][input.cols()-j-1];
        }
    }
    return output;
}

// The rest of these functions are incomplete. These are just suggestions to
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
    SDoublePlane output(input.rows(), input.cols());
    
    SDoublePlane flipped_filter = flipper(row_filter);
    
    int filterSize = flipped_filter.rows();
    
    int k = (filterSize-1)/2;
    
    for (int i = 0; i < input.rows(); i++)
    {
        for (int j = 0; j < input.cols(); j++)
        {
            if(i >= k && j >= k && ((i+k) < input.rows()) && (j+k) < input.cols())
            {
                for (int ki = -k; ki <= k; ki++ )
                {
                    for (int kj = -k; kj<= k; kj++)
                    {
                        output[i][j] = output[i][j] + flipped_filter[ki+k][kj+k] * input[i - ki][j - kj];
                    }
                }
            }
            else
            {
                output[i][j] = 255;
            }
            
        }
    }
    return output;
}



SDoublePlane convolve_edge(const SDoublePlane &input, const SDoublePlane &filter)
{
    SDoublePlane output(input.rows(), input.cols());
    
    SDoublePlane flipped_filter = flipper(filter);
    
    int filterSize = flipped_filter.rows();
    
    int k = (filterSize-1)/2;
    
    for (int i = 0; i < input.rows(); i++)
    {
        for (int j = 0; j < input.cols(); j++)
        {
            if(i >= k && j >= k && ((i+k) < input.rows()) && (j+k) < input.cols())
            {
                for (int ki = -k; ki <= k; ki++ )
                {
                    for (int kj = -k; kj<= k; kj++)
                    {
                        output[i][j] = output[i][j] + flipped_filter[ki+k][kj+k] * input[i - ki][j - kj];
                    }
                }
            }
            else
            {
                output[i][j] = 255;
            }
            
        }
    }
    return output;
}

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
    SDoublePlane output(input.rows(), input.cols());
    
    SDoublePlane flipped_filter = flipper(filter);
    
    SDoublePlane temp(input.rows()+2,input.cols()+2);
    
    for (int i = 0; i < temp.rows(); i++)
    {
        for (int j = 0; j < temp.cols(); j++)
        {
            
            if( i>= 1 && j>= 1 && (i < (temp.rows()-1) && (j<temp.cols()-1)))
            {
                temp[i][j] = input[i-1][j-1];
            }
            else
            {
                temp[i][j] = 255;
            }
        }
    }
    
    for (int i = 1; i < input.rows()+1; i++)
    {
        for (int j = 1; j < input.cols()+1; j++)
        {
            
            output[i-1][j-1] = temp[i-1][j-1] * flipped_filter[0][0] + temp[i][j-1] * flipped_filter[1][0] + temp[i+1][j-1] * flipped_filter[2][0] + temp[i-1][j] * flipped_filter[0][1] + temp[i][j] * flipped_filter[1][1] + temp[i+1][j] * flipped_filter[2][1] + temp[i-1][j+1] * flipped_filter[0][2] + temp[i][j+1] * flipped_filter[0][2] + temp[i+1][j+1] * flipped_filter[2][2];
            
        }
    }
    
    return output;
}

SDoublePlane convolve_general_4(const SDoublePlane &input, const SDoublePlane &filter)
{
    SDoublePlane output(input.rows(), input.cols());
    
    SDoublePlane flipped_filter = flipper(filter);
    
    int filterSize = flipped_filter.rows();
    
    int k = (filterSize-1)/2;
    
    for (int i = k; i < input.rows()-k; i++)
    {
        for (int j = k; j < input.cols()-k; j++)
        {
            
            for (int ki = -k; ki <= k; ki++ )
            {
                for (int kj = -k; kj<= k; kj++)
                {
                    output[i][j] = output[i][j] + flipped_filter[ki+k][kj+k] * input[i - ki][j - kj];
                }
            }
        }
    }
    return output;
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
                output[i][j] = 255;
            }
            
        }
    }
    return output;
}


SDoublePlane sobelSqRt(SDoublePlane inputX, SDoublePlane inputY)
{
    SDoublePlane output(inputX.rows(), inputX.cols());
    
    for (int i = 0; i < inputX.rows(); i++)
    {
        for (int j = 0; j < inputX.cols(); j++)
        {
            output[i][j] = sqrt((pow(inputX[i][j],2)+pow(inputY[i][j], 2)));
        }
    }
    return output;
    
}

// Apply a sobel operator to an image, returns the result
//
SDoublePlane sobel_gradient_filter(const SDoublePlane &input)
{
    SDoublePlane output(input.rows(), input.cols());
    
    SDoublePlane sobelFilterX(3,3);
    
    sobelFilterX[0][0]=-3;
    sobelFilterX[0][1]=0;
    sobelFilterX[0][2]=3;
    sobelFilterX[1][0]=-10;
    sobelFilterX[1][1]=0;
    sobelFilterX[1][2]=10;
    sobelFilterX[2][0]=-3;
    sobelFilterX[2][1]=0;
    sobelFilterX[2][2]=3;
    
    SDoublePlane sobelFilterY(3,3);
    
    sobelFilterY[0][0]=-3;
    sobelFilterY[0][1]=-10;
    sobelFilterY[0][2]=-3;
    sobelFilterY[1][0]=0;
    sobelFilterY[1][1]=0;
    sobelFilterY[1][2]=0;
    sobelFilterY[2][0]=3;
    sobelFilterY[2][1]=10;
    sobelFilterY[2][2]=3;
    
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            sobelFilterX[i][j] = sobelFilterX[i][j] /32.0;
            sobelFilterY[i][j] = sobelFilterY[i][j] /32.0;
            
        }
    }
    
    SDoublePlane sobelX = convolve_edge(input, sobelFilterX);
    
    SDoublePlane sobelY = convolve_edge(input, sobelFilterY);
    
    SDoublePlane sobelOuput = sobelSqRt(sobelX,sobelY);
    
    
    int maxValue = 0;
    
    for(int i = 0 ; i < sobelOuput.rows(); i++)
    {
        for(int j = 0 ; j < sobelOuput.cols(); j++)
        {
            if(sobelOuput[i][j]> maxValue)
            {
                maxValue = sobelOuput[i][j];
            }
            
        }
    }
    
    output = binaryImgGen(sobelOuput, maxValue/5);
    
    return output;
}

SDoublePlane convolve_template(const SDoublePlane &input, const SDoublePlane &filter)
{
    SDoublePlane output(input.rows(), input.cols());
    
    SDoublePlane flipped_filter = flipper(filter);
    
    SDoublePlane temp(input.rows()+2,input.cols()+2);
    
    
    for (int i = 0; i < temp.rows(); i++)
    {
        for (int j = 0; j < temp.cols(); j++)
        {
            if((i>=1 || j>=1) && (i < (temp.rows()-2) || (j<temp.cols()-2)))
            {
                temp[i][j] = input[i-1][j-1];
            }
            else
            {
                temp[i][j] = 255;
            }
        }
    }
    /*
     int filterSize = flipped_filter.rows();
     
     int k = (filterSize-1)/2;
     
     for (int i = 0; i < input.rows(); i++)
     {
     for (int j = 0; j < input.cols(); j++)
     {
     if(i >= k && j >= k && ((i+k) < input.rows()) && (j+k) < input.cols())
     {
     for (int ki = -k; ki <= k; ki++ )
     {
     for (int kj = -k; kj<= k; kj++)
     {
     output[i][j] = output[i][j] + flipped_filter[ki+k][kj+k] * input[i - ki][j - kj];
     }
     }
     }
     else
     {
     output[i][j] = 255;
     }
     
     }
     }*/
    return output;
}


SDoublePlane inverse(const SDoublePlane &input )
{
    SDoublePlane output(input.rows(),input.cols());
    //Logic for inverse has to be filled
    for(int i=0;i<input.rows();i++)
    {
        for(int j=0;j<input.cols();j++)
        {
            if (input[i][j] == 255)
                output[i][j] = 0;
            else
                output[i][j] = 255;
        }
    }
    return output;
}

// Apply an edge detector to an image, returns the binary edge map
//
SDoublePlane find_edges(const SDoublePlane &input, double thresh, const SDoublePlane &gFilter)
{
    SDoublePlane output(input.rows(), input.cols());
    
    // Implement an edge detector of your choice, e.g.
    // use your sobel gradient operator to compute the gradient magnitude and threshold
    
    SDoublePlane output_image = convolve_general(input, gFilter);
    
    SImageIO::write_png_file("Blurred.png", output_image, output_image, output_image);
    
    SDoublePlane output_image1 = binaryImgGen(output_image, thresh);
    
    SImageIO::write_png_file("Threshold.png", output_image1, output_image1, output_image1);
    
    SDoublePlane output_image2 = sobel_gradient_filter(output_image1);
    
    output = inverse(output_image2);
    
    return output;
}




SDoublePlane createFilter()
{
    SDoublePlane gFilter(3,3);
    // initialization of standard deviation to 1.0
    double sigma = 3.0;
    double r, s = 1.0 * sigma * sigma;
    
    // sum is for normalization
    double sum = 0.0;
    
    for (int x = -1; x <= 1; x++)
    {
        for(int y = -1; y <= 1; y++)
        {
            r = sqrt(x*x + y*y);
            gFilter[x + 1][y + 1] = (exp(-(r*r)/s))/(M_PI * s);
            sum += gFilter[x + 1][y + 1];
        }
    }
    
    // normalize the Kernel
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            gFilter[i][j] /= sum;
    
    return gFilter;
    
}

int maximum(SDoublePlane &input)
{
    int max = 0;
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

vector<DetectedSymbol> symDetectionByTemplate(const SDoublePlane &input_image, const SDoublePlane &template_gen, const string &str, const vector<Dimensions> &dim)
{
    vector<DetectedSymbol> symbols;
    SDoublePlane gFilter(5,5);
    gFilter = createFilter();
    double thresh = 100;
    int inp_row = input_image.rows();
    int inp_col = input_image.cols();
    int templ_row = template_gen.rows();
    int templ_col = template_gen.cols();
    
    //compute the Hamming distance function
    SDoublePlane binary_template_1 = binaryImgGen(template_gen, thresh);
    SDoublePlane binary_input_image = binaryImgGen(input_image, thresh);
    SDoublePlane inverse_template_1 = inverse(binary_template_1);
    SDoublePlane flipped_inverse_template_1 = flipper(inverse_template_1);
    SDoublePlane flipped_template_1 = flipper(binary_template_1);
    SDoublePlane inverse_input_image = inverse(binary_input_image);
    SDoublePlane F = convolve_general_4(input_image,flipped_template_1) + convolve_general_4(inverse_input_image,flipped_inverse_template_1);
    
    DetectedSymbol s;
    //find the maximum value in F
    int max;
    max = maximum(F);
    int prev_row = 0;
    int prev_col = 0;
    //find indexes of maxima in F
    // These indexes will be the location where the template is most likely to be present.
    
    for(int i= 0; i< inp_row; i++)
    {
        for(int j=0; j< inp_col; j++)
        {
            if (F[i][j] >= (0.95*max))
            {
                //cout << " Row: " << i << " Col: " << j << " " << F[i][j] << endl;
                s.row = i - int(ceil(templ_row/2)) - 2;
                s.col = j - int(ceil(templ_col/2)) - 4;
                s.width = templ_col + 3;
                s.height = templ_row + 3;
                if (str == "NOTEHEAD")   //There should be a factor added to this when we scale up or scale down the template images
                {
                    s.type = (Type) 0;
                }
                else if (str == "QUARTERREST")  //There should be a factor added to this when we scale up or scale down the template images
                {
                    s.type = (Type) 1;
                }
                else
                {
                    s.type = (Type) 2;
                }
                
                s.confidence = 1 - ((max - F[i][j])/(0.1*max)) ;
                
                if ((abs(prev_row-s.row) + abs(prev_col - s.col))> 2)
                {
                    int starting_trebble_staff = 0,starting_bass_staff = 0;
                    int h = 0,min = 10000000,dist = 0;
                    
                    for(int k = 0;k<dim.size();k++)
                    {
                        const Dimensions &d = dim[k];
                        dist = abs(i - d.row_coordinate);
                        if (dist < min)
                        {
                            min = dist;
                            if (d.trebble)
                            {
                                starting_trebble_staff = d.row_coordinate;
                                h = d.spacing;
                                starting_bass_staff = 0;
                            }
                            else
                            {
                                starting_bass_staff = d.row_coordinate;
                                h = d.spacing;
                                starting_trebble_staff = 0;
                            }
                        }
                    }
                    if (starting_trebble_staff)
                    {
                        if (((i >(starting_trebble_staff - (4*h))) && (i < starting_trebble_staff - (3*h))) \
                            || (starting_trebble_staff == i) \
                            || ((i > (starting_trebble_staff + (3*h))) && (i < starting_trebble_staff + (4*h))))
                        {
                            s.pitch= 'F';
                        }
                        
                        if ((i == (starting_trebble_staff - (3*h))) \
                            || ((i > starting_trebble_staff) && (i < (starting_trebble_staff+h))) \
                            || (i == (starting_trebble_staff + (4*h))))
                        {
                            s.pitch= 'E';
                        }
                        
                        if (((i >(starting_trebble_staff - (5*h))) && (i < starting_trebble_staff - (4*h))) \
                            || (i == starting_trebble_staff + h ) \
                            || ((i > (starting_trebble_staff + (4*h))) && (i < starting_trebble_staff + (5*h))))
                        {
                            s.pitch= 'D';
                        }
                        
                        if ((i == (starting_trebble_staff - (2*h))) \
                            || ((i > starting_trebble_staff + h) && (i < (starting_trebble_staff+(2*h)))) \
                            || (i == (starting_trebble_staff + (5*h))))
                        {
                            s.pitch= 'C';
                        }
                        
                        if (((i >(starting_trebble_staff - (6*h))) && (i < starting_trebble_staff - (5*h))) \
                            || (i == starting_trebble_staff + (2*h)) \
                            || ((i > (starting_trebble_staff + (5*h))) && (i < starting_trebble_staff + (6*h))))
                        {
                            s.pitch= 'B';
                        }
                        
                        if ((i == (starting_trebble_staff - h)) \
                            || ((i > (starting_trebble_staff + (2*h))) && (i < (starting_trebble_staff + (3*h)))) )
                        {
                            s.pitch = 'A';
                        }
                        
                        if (((i >(starting_trebble_staff - (7*h))) && (i < starting_trebble_staff - (6*h))) \
                            || (i == starting_trebble_staff + (3*h)) \
                            || ((i > (starting_trebble_staff + (6*h))) && (i < starting_trebble_staff + (7*h))))
                        {
                            s.pitch= 'G';
                        }
                    }
                    else
                    {
                        if (((i > (starting_bass_staff+(4*h))) && (i < (starting_bass_staff+(5*h)))) \
                            || ( i == (starting_bass_staff + h)) \
                            || ((i > (starting_bass_staff- (2*h))) && (i < (starting_bass_staff- (3*h)))))
                        {
                            s.pitch= 'F';
                        }
                        
                        if (((i > (starting_bass_staff + h)) && (i < (starting_bass_staff+(2*h)))) \
                            || ( i == (starting_bass_staff + (6*h))) \
                            || ( i == (starting_bass_staff - (2*h))))
                        {
                            s.pitch= 'E';
                        }
                        
                        if (((i > (starting_bass_staff+(5*h))) && (i < (starting_bass_staff+(6*h)))) \
                            || ( i == (starting_bass_staff + (2*h))) \
                            || ((i > (starting_bass_staff - h)) && (i < (starting_bass_staff- (2*h)))))
                        {
                            s.pitch= 'D';
                        }
                        
                        if (( i == (starting_bass_staff + h)) \
                            || (( i > (starting_bass_staff - (2*h))) && (i < (starting_bass_staff-(3*h)))) \
                            || (i == (starting_bass_staff- (6*h))))
                        {
                            s.pitch= 'C';
                        }
                        
                        if (((i > (starting_bass_staff - h)) && (i < (starting_bass_staff))) \
                            || ( i == (starting_bass_staff + (4*h))) \
                            || ((i > (starting_bass_staff + (7*h))) && (i < (starting_bass_staff + (8*h)))))
                        {
                            s.pitch= 'B';
                        }
                        
                        if ((i == starting_bass_staff) \
                            || ((i > (starting_bass_staff - (3*h))) && (i < (starting_bass_staff - (4*h)))) \
                            || ((i > (starting_bass_staff + (3*h))) && (i < (starting_bass_staff + (4*h)))))
                        {
                            s.pitch = 'A';
                        }
                        
                        if (((i > (starting_bass_staff)) && (i < (starting_bass_staff+h))) \
                            || ( i == (starting_bass_staff + (5*h))) \
                            || ((i == (starting_bass_staff - (3*h))) ))
                        {
                            s.pitch= 'G';
                        }
                    }
                    
                    symbols.push_back(s);
                    
                }
                
                prev_row = s.row;
                prev_col = s.col;
            }
        }
    }
    
    return symbols;
}


vector<DetectedSymbol> symDetectionAfterEdges(const SDoublePlane &input_image, const SDoublePlane &template_gen, const string &str, const vector<Dimensions> &dim)
{
    vector<DetectedSymbol> symbols;
    
    SDoublePlane gammaFunctionOutput(input_image.rows(),input_image.cols());
    
    SDoublePlane dOutput(input_image.rows(),input_image.cols());
    
    SDoublePlane F(input_image.rows(),input_image.cols());
    
    double thresh = 100;
    int inp_row = input_image.rows();
    int inp_col = input_image.cols();
    int templ_row = template_gen.rows();
    int templ_col = template_gen.cols();
    

    
    for(int i=0; i<input_image.rows();i++)
    {
        for(int j=0; j<input_image.cols();j++)
        {
            if(input_image[i][j]>0)
            {
                gammaFunctionOutput[i][j] = 0;
            }
            else
            {
                gammaFunctionOutput[i][j] = std::numeric_limits<double>::infinity();
            }
        }
    }
    
    
    double tempValue = std::numeric_limits<double>::infinity();
    
    double Value = 0;
    
    
    for(int i=0; i<input_image.rows();i++)
    {
        for(int j=0; j<input_image.cols();j++)
        {
            tempValue = std::numeric_limits<double>::infinity();
            
            if(gammaFunctionOutput[i][j] != 0)
            {
                for(int a=0; a<input_image.rows(); a++)
                {
                    for(int b=0; b<input_image.cols(); b++)
                    {
                        Value = gammaFunctionOutput[a][b] + sqrt(pow((i-a),2) + pow((j-b), 2) );
                        
                        if(Value < tempValue)
                        {
                            tempValue = Value;
                        }
                    }
                }
                dOutput[i][j] = tempValue;
            }
        }
    }
    
    for(int i=0; i<input_image.rows()-template_gen.rows();i++)
    {
        for(int j=0; j<input_image.cols()-template_gen.cols();j++)
        {
            
            for(int k = 0; k<template_gen.rows()-1; k++)
            {
                for(int l = 0; l<template_gen.cols()-1; l++)
                {
                    F[i][j] = F[i][j] + template_gen[k][l]*dOutput[i+k][j+l];
                }
            }
        }
    }
    
    
    DetectedSymbol s;
    //find the maximum value in F
    int max;
    max = maximum(F);
    int prev_row = 0;
    int prev_col = 0;
    //find indexes of maxima in F
    // These indexes will be the location where the template is most likely to be present.
    
    for(int i= 0; i< inp_row; i++)
    {
        for(int j=0; j< inp_col; j++)
        {
            if (F[i][j] >= (0.95*max))
            {
                //cout << " Row: " << i << " Col: " << j << " " << F[i][j] << endl;
                s.row = i - int(ceil(templ_row/2)) - 2;
                s.col = j - int(ceil(templ_col/2)) - 4;
                s.width = templ_col + 3;
                s.height = templ_row + 3;
                if (str == "NOTEHEAD")   //There should be a factor added to this when we scale up or scale down the template images
                {
                    s.type = (Type) 0;
                }
                else if (str == "QUARTERREST")  //There should be a factor added to this when we scale up or scale down the template images
                {
                    s.type = (Type) 1;
                }
                else
                {
                    s.type = (Type) 2;
                }
                
                s.confidence = 1 - ((max - F[i][j])/(0.1*max)) ;
                
                if ((abs(prev_row-s.row) + abs(prev_col - s.col))> 2)
                {
                    int starting_trebble_staff = 0,starting_bass_staff = 0;
                    int h = 0,min = 10000000,dist = 0;
                    
                    for(int k = 0;k<dim.size();k++)
                    {
                        const Dimensions &d = dim[k];
                        dist = abs(i - d.row_coordinate);
                        if (dist < min)
                        {
                            min = dist;
                            if (d.trebble)
                            {
                                starting_trebble_staff = d.row_coordinate;
                                h = d.spacing;
                                starting_bass_staff = 0;
                            }
                            else
                            {
                                starting_bass_staff = d.row_coordinate;
                                h = d.spacing;
                                starting_trebble_staff = 0;
                            }
                        }
                    }
                    if (starting_trebble_staff)
                    {
                        if (((i >(starting_trebble_staff - (4*h))) && (i < starting_trebble_staff - (3*h))) \
                            || (starting_trebble_staff == i) \
                            || ((i > (starting_trebble_staff + (3*h))) && (i < starting_trebble_staff + (4*h))))
                        {
                            s.pitch= 'F';
                        }
                        
                        if ((i == (starting_trebble_staff - (3*h))) \
                            || ((i > starting_trebble_staff) && (i < (starting_trebble_staff+h))) \
                            || (i == (starting_trebble_staff + (4*h))))
                        {
                            s.pitch= 'E';
                        }
                        
                        if (((i >(starting_trebble_staff - (5*h))) && (i < starting_trebble_staff - (4*h))) \
                            || (i == starting_trebble_staff + h ) \
                            || ((i > (starting_trebble_staff + (4*h))) && (i < starting_trebble_staff + (5*h))))
                        {
                            s.pitch= 'D';
                        }
                        
                        if ((i == (starting_trebble_staff - (2*h))) \
                            || ((i > starting_trebble_staff + h) && (i < (starting_trebble_staff+(2*h)))) \
                            || (i == (starting_trebble_staff + (5*h))))
                        {
                            s.pitch= 'C';
                        }
                        
                        if (((i >(starting_trebble_staff - (6*h))) && (i < starting_trebble_staff - (5*h))) \
                            || (i == starting_trebble_staff + (2*h)) \
                            || ((i > (starting_trebble_staff + (5*h))) && (i < starting_trebble_staff + (6*h))))
                        {
                            s.pitch= 'B';
                        }
                        
                        if ((i == (starting_trebble_staff - h)) \
                            || ((i > (starting_trebble_staff + (2*h))) && (i < (starting_trebble_staff + (3*h)))) )
                        {
                            s.pitch = 'A';
                        }
                        
                        if (((i >(starting_trebble_staff - (7*h))) && (i < starting_trebble_staff - (6*h))) \
                            || (i == starting_trebble_staff + (3*h)) \
                            || ((i > (starting_trebble_staff + (6*h))) && (i < starting_trebble_staff + (7*h))))
                        {
                            s.pitch= 'G';
                        }
                    }
                    else
                    {
                        if (((i > (starting_bass_staff+(4*h))) && (i < (starting_bass_staff+(5*h)))) \
                            || ( i == (starting_bass_staff + h)) \
                            || ((i > (starting_bass_staff- (2*h))) && (i < (starting_bass_staff- (3*h)))))
                        {
                            s.pitch= 'F';
                        }
                        
                        if (((i > (starting_bass_staff + h)) && (i < (starting_bass_staff+(2*h)))) \
                            || ( i == (starting_bass_staff + (6*h))) \
                            || ( i == (starting_bass_staff - (2*h))))
                        {
                            s.pitch= 'E';
                        }
                        
                        if (((i > (starting_bass_staff+(5*h))) && (i < (starting_bass_staff+(6*h)))) \
                            || ( i == (starting_bass_staff + (2*h))) \
                            || ((i > (starting_bass_staff - h)) && (i < (starting_bass_staff- (2*h)))))
                        {
                            s.pitch= 'D';
                        }
                        
                        if (( i == (starting_bass_staff + h)) \
                            || (( i > (starting_bass_staff - (2*h))) && (i < (starting_bass_staff-(3*h)))) \
                            || (i == (starting_bass_staff- (6*h))))
                        {
                            s.pitch= 'C';
                        }
                        
                        if (((i > (starting_bass_staff - h)) && (i < (starting_bass_staff))) \
                            || ( i == (starting_bass_staff + (4*h))) \
                            || ((i > (starting_bass_staff + (7*h))) && (i < (starting_bass_staff + (8*h)))))
                        {
                            s.pitch= 'B';
                        }
                        
                        if ((i == starting_bass_staff) \
                            || ((i > (starting_bass_staff - (3*h))) && (i < (starting_bass_staff - (4*h)))) \
                            || ((i > (starting_bass_staff + (3*h))) && (i < (starting_bass_staff + (4*h)))))
                        {
                            s.pitch = 'A';
                        }
                        
                        if (((i > (starting_bass_staff)) && (i < (starting_bass_staff+h))) \
                            || ( i == (starting_bass_staff + (5*h))) \
                            || ((i == (starting_bass_staff - (3*h))) ))
                        {
                            s.pitch= 'G';
                        }
                    }
                    
                    symbols.push_back(s);
                    
                }
                
                prev_row = s.row;
                prev_col = s.col;
            }
        }
    }
    return symbols;
}

vector<Dimensions> findStaff(const SDoublePlane &input)
{
    int inp_row = input.rows();
    int inp_col = input.cols();
    int accum_d = inp_row/10;
    int accum_r = inp_row-10;
    vector<Dimensions> staves;
    _DTwoDimArray<int> accum(accum_r, accum_d);
    
    for(int i=0;i<accum_r;i++)
    {
        for(int j=0;j<accum_d;j++)
        {
            accum[i][j] = 0;
        }
    }
    //Voting the Accumulator Space based on the pixel values.
    for(int i=1;i<accum_r;i++)
    {
        for(int j=0;j<inp_col;j++)
        {
            for(int h=1;h<=accum_d;h++)
            {
                if ((i +(4*h)+1) < inp_row-1)
                {
                    if ((input[i][j] == 0 || input[i+1][j] ==0 || input[i-1][j] ==0) && \
                        (input[i+h][j] == 0 || input[i+h+1][j] == 0 || input[i+h-1][j] == 0) && \
                        (input[i+(2*h)][j] == 0 || input[i+(2*h)+1][j] == 0 || input[i+(2*h)-1][j] == 0) && \
                        (input[i+(3*h)][j] == 0 || input[i+(3*h)+1][j] == 0 || input[i+(3*h)-1][j] == 0) && \
                        (input[i+(4*h)][j] == 0 || input[i+(4*h)+1][j] == 0 || input[i+(4*h)-1][j] == 0))
                    {
                        accum[i][h]++;
                    }
                }
            }
        }
    }
    //Finding the Parameters surpassing a threshold
    Dimensions dim,prev_dim;
    bool trebble = 1;
    for(int i=0;i<accum_r;i++)
    {
        for(int j=0;j<accum_d;j++)
        {
            if (accum[i][j] >= 0.9*inp_col)
            {
                dim.row_coordinate =i;
                dim.spacing = j;
                dim.trebble = trebble;
                
                if ((abs(prev_dim.row_coordinate - dim.row_coordinate) > 5) && (dim.spacing > 3))
                {
                    staves.push_back(dim);
                    trebble = !trebble;
                }
                //cout << dim.row_coordinate << " " << dim.spacing << " " << abs(prev_dim.row_coordinate - dim.row_coordinate) << endl;
                prev_dim.row_coordinate = dim.row_coordinate;
                prev_dim.spacing = dim.spacing;
                break;
            }
        }
    }
    return staves;
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
    
    // + Template File names
    string template_1_str = "template1.png";
    string template_2_str = "template2.png";
    string template_3_str = "template3.png";
    // - Template File names
    
    // + Read data from the image file as well as templates
    SDoublePlane input_image = SImageIO::read_png_file(input_filename.c_str());
    SDoublePlane template_1 = SImageIO::read_png_file(template_1_str.c_str());
    SDoublePlane template_2 = SImageIO::read_png_file(template_2_str.c_str());
    SDoublePlane template_3 = SImageIO::read_png_file(template_3_str.c_str());
    // - Read data from the image file as well as templates
    
    // + Mean filter given part of skeleton code
    SDoublePlane mean_filter(3,3);
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            mean_filter[i][j] = 1/9.0;
    // - Mean filter given part of skeleton code
    
    // + Gaussian filter generated for smoothening image and reducing noise
    SDoublePlane gFilter(3,3);
    gFilter = createFilter();
    // - Gaussian filter generated for smoothening image and reducing noise
    
    // + Finding the edge map for all the images
    SDoublePlane output_image1 = find_edges(input_image, 200, gFilter);
    
    SImageIO::write_png_file("edges.png", output_image1, output_image1, output_image1);
    
    SDoublePlane template1 = find_edges(template_1, 200, gFilter);
    
    SImageIO::write_png_file("Edge-template1.png", template1, template1, template1);
    
    SDoublePlane template2 = find_edges(template_2, 200, gFilter);
    
    SImageIO::write_png_file("Edge-template2.png", template2, template2, template2);
    
    SDoublePlane template3 = find_edges(template_3, 200, gFilter);
    
    SImageIO::write_png_file("Edge-template3.png", template3, template3, template3);
    // - Finding the edge map for all the images
    
    // + Finding Staves
    SDoublePlane output_image3 = binaryImgGen(input_image, 150);
    vector<Dimensions> staff_lines = findStaff(output_image3);
    write_staff_detection_image("staves.png", staff_lines, input_image);
    // - Finding Staves
    
    // + Step 4 - Template matching method
    vector<DetectedSymbol> symbols_4 = symDetectionByTemplate(input_image,template_1,"NOTEHEAD", staff_lines);
    write_detection_image("detected4.png", symbols_4, input_image);
    // - Step 4 - Template matching method
    
    // + Step 5 - Template matching method after edge detection
    vector<DetectedSymbol> symbols1 = symDetectionAfterEdges(input_image,template_1,"NOTEHEAD", staff_lines);
    write_detection_image("detected5.png", symbols1, input_image);
    // - Step 5 - Template matching method after edge detection
    
    
    /*for(int i=0; i<10; i++)
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
     */
    
    // + Writing final output
    write_detection_txt("detected7.txt", symbols_4);
    write_detection_image("detected7.png", symbols_4, input_image);
    // - Writing final output
}
