#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include<algorithm>
#include<random>
#include<iostream>

static std::default_random_engine engine(10);
static std::normal_distribution<double> normal(0,1);


void color_matching(double* image_double, double* image_target_double, double* image_result_double, int W, int H){

	memcpy(&image_result_double, &image_double, W*H*3*sizeof(double));

	for (int nslices=0; nslices < 1; nslices++){

		double x = normal(engine);
		double y = normal(engine);
		double z = normal(engine);
		std::vector<std::pair<double,int>> proj1(W*H);
		std::vector<std::pair<double,int>> proj2(W*H);

		double n = sqrt(x * x + y * y + z * z);
		x/=n;
		y/=n;
		z/=n;
		for (int i = 0; i < W * H; i ++){
			double p1 = image_result_double[i*3+0]*x + image_result_double[i*3+1]*y + image_result_double[i*3+2]*z;
			double p2 = image_target_double[i*3+0]*x + image_target_double[i*3+1]*y + image_target_double[i*3+2]*z;
			proj1[i] = std::pair<double,int>(p1,1);
			proj2[i] = std::pair<double,int>(p2,1);
		}
		std::sort(proj1.begin(), proj1.end());
		std::sort(proj2.begin(), proj2.end());



			for (int i=0; i < W*H; i++){
				double distance_along_line = (proj2[i].first- proj1[i].first);

				image_result_double[proj1[i].second * 3 + 0] += x * distance_along_line;
				image_result_double[proj1[i].second * 3 + 1] += y * distance_along_line;
				image_result_double[proj1[i].second * 3 + 2] += z * distance_along_line;





			}




	}
}




int main() {

	int W, H, C;
	
	//stbi_set_flip_vertically_on_load(true);
	unsigned char *image = stbi_load("8733654151_b9422bb2ec_k.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);


	unsigned char *image_target = stbi_load("target_resized.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);




	std::vector<double> image_double(W*H*3);
	std::vector<double> image_target_double(W*H*3);

	for (int i=0; i<W*H*3; i++){
		image_double[i] = image[i];
		image_target_double[i] = image_target[i];
	}

	
	std::vector<double> image_result_double(W*H * 3);
	std::cout << "test";

	color_matching((double*) &image_double, (double*) &image_target,(double*) &image_result_double, W, H);
	std::cout << "test";


	std::vector<unsigned char> image_result(W*H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

	image_result[(i*W + j) * 3 + 0] = std::min(255., std::max(0., image_result_double[(i*W+j)*3+0]));
    image_result[(i*W + j) * 3 + 1] = std::min(255., std::max(0., image_result_double[(i*W+j)*3+1]));
    image_result[(i*W + j) * 3 + 2] = std::min(255., std::max(0., image_result_double[(i*W+j)*3+2]));
		}
	}

	stbi_write_png("test_image.png", W, H, 3, &image_result[0], 0);

	return 0;
}