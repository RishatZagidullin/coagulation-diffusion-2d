#include "ppm.h"

using namespace std;


int create_ppm(int N, int M, double const * const data, std::string filename, double maxi) {
    int img_N = 256;
    int img_M = 256;
    size_t SIZE = img_N * img_M;
    int * image_data = new int [SIZE];

    ofstream myImage;
    myImage.open(filename);

    myImage << "P2" << endl;
    myImage << img_N << " " << img_M << endl;
    myImage << "255" << endl;

    int pixel = 0;
    for (int row = 0; row < img_N; row++) {
        for (int col = 0; col < img_M; col++) {

            int row_id = (int) (((float) row/ (float) img_N) * N);
            int col_id = (int) (((float) col/ (float) img_M) * M);

            double val = data[col_id+row_id*M]/maxi;

            int res = (int) ( (log10(val+1e-12)+12)/12 * 255);

            //if (res > 255)
            //    std::cout << "pixel more than 255, shouldn't happen\n";
            image_data[pixel] = res > 255 ? 255 : res;

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
