#include "ppm.h"

using namespace std;


int create_ppm(int N, int M, double const * const data, std::string filename) {
    int img_N = 256;
    int img_M = 256;
    size_t SIZE = 256 * 256 * 3;
    int * image_data = new int [SIZE];

    ofstream myImage;
    myImage.open(filename);

    myImage << "P3" << endl;
    myImage << img_N << " " << img_M << endl;
    myImage << "255" << endl;

    int pixel = 0;
    for (int row = 0; row < img_N; row++) {
        for (int col = 0; col < img_M; col++) {
            int red{0};
            int green{0};
            int blue{0};

            int row_id = (int) (((float) row/ (float) img_N) * N);
            int col_id = (int) (((float) col/ (float) img_M) * M);
            double color = data[col_id+row_id*M];
            if (data[col_id+row_id*M] < 0.0001)
            {
                green = (-4./log10(color+1e-8) * 255);
                blue = (-log10(color+1e-8)/8 * 255);
            }
            else
            {
                red = ((log10(color)/4.+1.) * 255);
                green = (-log10(color)/4. * 255);
            }
            image_data[pixel*3] = red > 255 ? 255 : red;
            image_data[pixel*3 + 1] = green > 255 ? 255 : green;
            image_data[pixel*3 + 2] = blue > 255 ? 255 : blue; 

            pixel++;
        }
    }
    for (int x = 0; x < SIZE; x++) {
        int value = image_data[x];
        myImage << value << " " << endl;
    }

    myImage.close();

    delete [] image_data;

    return 0;
}
