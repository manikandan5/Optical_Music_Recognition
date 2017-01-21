# Optical Music Recognition

This code was written for the course Computer Vision (CSCI-B 657) at Indiana University handled by Professor David Crandall. Skeleton code was provided by the Professor to get us started with the assignment.


**What does the program do?** <br/>
* The program tries to find the musical notes from the music sheet.

**How does it find it?** <br/>
* Initially we try to remove the noise in the image by convoluting the image with a Gaussian filter.
* Then we find the edges in the image by convoluting with a Sobel filter.
* Then the staff lines are found using Hough Transformation.
* Once the staff lines are found, we find on which staff line the nodes are present using template matching.
* With this data, the program tries to estimate the correct note.


Detailed explanation about how the code works and the reason why we chose this implementation could be found [here](https://github.com/manikandan5/Optical_Music_Recognition/blob/master/Report.pdf).

**How to run the program?** 

This command compiles the program:
* make 

To run the program:
* ./omr image.png
