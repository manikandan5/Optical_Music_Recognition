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

void overlay_line(SDoublePlane &input, int _top, int _space, double graylevel, int width)
{
    for(int w=-width/2; w<=width/2; w++) {
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
class DetectedSymbol {
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

SDoublePlane convolve_edge(const SDoublePlane &input, const SDoublePlane &filter)
{
    SDoublePlane output(input.rows(), input.cols());
    
    
    return output;
}
// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
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
    
    SDoublePlane sobelX = convolve_general(input, sobelFilterX);
    
    SDoublePlane sobelY = convolve_general(input, sobelFilterY);
    
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

vector<DetectedSymbol> symDetectionByTemplate(const SDoublePlane &input_image,const SDoublePlane &template_1,const vector<Dimensions> &dim)
{
    vector<DetectedSymbol> symbols;
    SDoublePlane gFilter(5,5);
    gFilter = createFilter();
    double thresh = 100;
    int inp_row = input_image.rows();
    int inp_col = input_image.cols();

    //compute the Hamming distance function
    SDoublePlane binary_template_1 = binaryImgGen(template_1, thresh);
    SDoublePlane binary_input_image = binaryImgGen(input_image, thresh);
    SDoublePlane inverse_template_1 = inverse(binary_template_1);
    SDoublePlane flipped_inverse_template_1 = flipper(inverse_template_1);
    SDoublePlane flipped_template_1 = flipper(binary_template_1);
    SDoublePlane inverse_input_image = inverse(binary_input_image);
    SDoublePlane F = convolve_general(input_image,flipped_template_1) + convolve_general(inverse_input_image,flipped_inverse_template_1);
    
    //find the maximum value in F
    int max;
    max = maximum(F);
    
    //find indexes of maxima in F
    // These indexes will be the location where the template is most likely to be present.

    for(int i= 0; i< inp_row; i++)
    {
        for(int j=0; j< inp_col; j++)
        {
            if (F[i][j] >= (0.95*max))
            {
              cout << " Row: " << i << " Col: " << j << " " << F[i][j] << endl;
            }
        }
    }
    //find indexes of maxima in F and affect them to symbols
    /*int u,v,w;
    for (int u=0; u<input_image.rows(); u++)
    {
        for (int v=0; v < input_image.cols(); v++)
            if (F[u][v]=max )
            {
                w=w+1;
                symbols[w].row = v;
                symbols[w].col = u ;
                //this part of the code is to identify the pitch and type
                //if ((starting__trebble_staff=u || starting_trebble_staff+3*h<u<starting_trebble_staff+4*h) && (type condition))
                //{symbols[w].pitch= "F";}
                //symbols[w].type = NOTEHEAD;
                //if ((starting_trebble_staff<u<starting_trebble_staff+h || u=starting_trebble_staff+4*h) && (type condition))
                //{symbols[w].pitch= "E" ;
                //symbols[w].type = NOTEHEAD ;}
                //if ((u=starting_trebble_staff+h || starting_trebble_staff+4*h <u< starting_trebble_staff+5*h) && (type condition))
                //{symbols[w].pitch="D";
                //symbols[w].type=NOTEHEAD ;}
                //if ((starting_trebble_staff+h<u<starting_trebble_staff+2*h || u=starting_trebble_staff+5*h) && (type condition))
                //{symbols[w].pitch= "C" ;
                //symbols[w].type = NOTEHEAD ;}
                //if ((u=starting_trebble_staff+2*h || starting_trebble_staff+5*h <u< starting_trebble_staff+6*h) && (type condition))
                //{symbols[w].pitch="B";
                //symbols[w].type=NOTEHEAD ;}
                //if ((starting_trebble_staff+2*h<u<starting_trebble_staff+3*h || u=starting_trebble_staff+6*h) && (type condition))
                //{symbols[w].pitch= "A" ;
                //symbols[w].type = NOTEHEAD ;}
                //if ((u=starting_trebble_staff+3*h || starting_trebble_staff+6*h <u< starting_trebble_staff+7*h) && (type condition))
                //{symbols[w].pitch="G";
                //symbols[w].type=NOTEHEAD ;}
            }
    }*/
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
    string template_1_str = "template1.png";
    SDoublePlane input_image = SImageIO::read_png_file(input_filename.c_str());
    SDoublePlane template_1 = SImageIO::read_png_file(template_1_str.c_str());
    // test step 2 by applying gaussian filters to the input image
    SDoublePlane mean_filter(3,3);
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            mean_filter[i][j] = 1/9.0;
    
    SDoublePlane gFilter(5,5);
    gFilter = createFilter();
    
    
    SDoublePlane output_image = convolve_general(input_image, gFilter);
    
    SImageIO::write_png_file("Blurred.png", output_image, output_image, output_image);
    
    SDoublePlane output_image1 = binaryImgGen(output_image, 200);
    
    SImageIO::write_png_file("Threshold.png", output_image1, output_image1, output_image1);
    
    SDoublePlane output_image2 = sobel_gradient_filter(output_image1);
    output_image2 = inverse(output_image2);
    SImageIO::write_png_file("Edge.png", output_image2, output_image2, output_image2);
    
    SDoublePlane output_image3 = binaryImgGen(input_image, 150);
    vector<Dimensions> staff_lines = findStaff(output_image3);
    write_staff_detection_image("staves_withoutsobel.png", staff_lines, input_image);

    //staff_lines = findStaff(output_image2);
    //write_staff_detection_image("staves_withsobel.png", staff_lines, input_image);

    //- Temporary Edge Detedtion
    
    // randomly generate some detected symbols -- you'll want to replace this
    //  with your symbol detection code obviously!
    vector<DetectedSymbol> symbols = symDetectionByTemplate(input_image,template_1,staff_lines);
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
