//============================================================================
// Name        : Barcode.cpp
// Author      : Gerard Tobin
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================

#include <stdio.h>

#include <future>
#include <string>
#include <iostream>
#include <time.h>
#include <exception>

#include <opencv2/opencv.hpp>

#include <stdio.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <list>
#include <thread>

using namespace std;
typedef unsigned char byte;
struct point
{
    int x;
    int y;
    int m;
};

void Line(byte **input, float x1, float y1, float x2, float y2, const int color)
{
    // Bresenham's line algorithm
    const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
    if (steep)
    {
        std::swap(x1, y1);
        std::swap(x2, y2);
    }

    if (x1 > x2)
    {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }

    const float dx = x2 - x1;
    const float dy = fabs(y2 - y1);

    float error = dx / 2.0f;
    const int ystep = (y1 < y2) ? 1 : -1;
    int y = (int)y1;

    const int maxX = (int)x2;

    for (int x = (int)x1; x <= maxX; x++)
    {
        if (steep)
        {
            input[y][x] = color;
        }
        else
        {
            input[x][y] = color;
        }

        error -= dy;
        if (error < 0)
        {
            y += ystep;
            error += dx;
        }
    }
}

void floodFill(byte **input, int x, int y, int k, int height, int width)
{

    vector<vector<int>> M;

    for (int i = 0; i < height; i++)
    {
        vector<int> row;
        for (int j = 0; j < width; j++)
        {
            row.push_back(input[i][j]);
        }
        M.push_back(row);
    }
    queue<pair<int, int>> nodeQ;
    nodeQ.push({x, y});
    int oldCol = M[x][y];
    while (!nodeQ.empty())
    {
        pair<int, int> currNode = nodeQ.front();
        nodeQ.pop();
        if (M[currNode.first][currNode.second] == oldCol)
        {
            M[currNode.first][currNode.second] = k;
            if (currNode.first > 0)
                nodeQ.push({currNode.first - 1, currNode.second});
            if (currNode.first < (M.size() - 1))
                nodeQ.push({currNode.first + 1, currNode.second});
            if (currNode.second > 0)
                nodeQ.push({currNode.first, currNode.second - 1});
            if (currNode.second < (M[0].size() - 1))
                nodeQ.push({currNode.first, currNode.second + 1});
        }
    }

    for (int i = 0; i < height; i++)
    {
        vector<int> row = M.at(i);
        for (int j = 0; j < width; j++)
        {
            input[i][j] = row.at(j);
        }
    }
}

void cropborder(byte **dilation3, int border, int height, int width)
{

    for (int i = height - border; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            dilation3[i][j] = 255;
        }
    }
    for (int i = 0; i <= border; i++)
    {
        for (int j = 0; j < width; j++)
        {
            dilation3[i][j] = 255;
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = width - border; j < width; j++)
        {
            dilation3[i][j] = 255;
        }
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < border; j++)
        {
            dilation3[i][j] = 255;
        }
    }
}

void printArray(byte **matrix, int height, int width, string new_dir_name, string new_file_name)
{
    cv::Mat output_temp(height, width, CV_8UC1);
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {

            output_temp.at<byte>(i, j) = matrix[i][j];
        }
    }

    imwrite(new_dir_name + new_file_name, output_temp);
}

/**
 * @brief 
 * 
 * @param orig 
 * @param img 
 * @param height 
 * @param width 
 * @param threshold_limit 
 * @param he 
 * @param wi 
 * @param new_dir_name 
 * @param print 
 */
cv::Rect image_process(byte **img, int height, int width, double threshold_limit, int &he, int &wi, string new_dir_name, byte **(&arr2), bool print)
{

    byte **rot = new byte *[height];
    byte **dilation = new byte *[height];
    byte **erosion = new byte *[height];
    byte **sobel = new byte *[height];
    byte **matrix = new byte *[height];
    int **label = new int *[height];

    for (int i = 0; i < height; i++)
    {
        rot[i] = new byte[width];
    }
    for (int i = 0; i < height; i++)
    {
        dilation[i] = new byte[width];
        erosion[i] = new byte[width];
        sobel[i] = new byte[width];
        label[i] = new int[width];
        matrix[i] = new byte[width];
    }

    for (int i = 0; i < height; i++)
    {

        for (int j = 0; j < width; j++)
        {

            matrix[i][j] = img[i][j];
        }
    }
    if (print)
    {
        printArray(matrix, height, width, new_dir_name, "greyscale.jpg");
    }

    int vertical = 0;
    int horizontal = 0;
    int grad = 0;

    int sobel_h[] = {-1, 0, 1,
                     -2, 0, 2,
                     -1, 0, 1};
    int sobel_v[] = {1, 2, 1,
                     0, 0, 0,
                     -1, -2, -1};

    for (int i = 1; i < height - 1; i++)
    {
        for (int j = 1; j < width - 1; j++)
        {
            int ul = matrix[i - 1][j - 1];
            int um = matrix[i - 1][j];
            int ur = matrix[i - 1][j + 1];
            int ml = matrix[i][j - 1];

            int mm = matrix[i][j];
            int mr = matrix[i][j + 1];
            int ll = matrix[i + 1][j - 1];
            int lm = matrix[i + 1][j];
            int lr = matrix[i + 1][j + 1];

            horizontal = ul * -1 + um * 0 + ur * 1 + ml * -2 + mm * 0 + mr * 2 + ll * -1 + lm * 0 + lr * 1;
            vertical = ul * 1 + 2 * um + ur * 1 + ml * 0 + mm * 0 + mr * 0 - 1 * ll + -2 * lm + lr * -1;
            grad = sqrt(horizontal * horizontal + vertical * vertical);
            sobel[i][j] = grad;
        }
    }
    if (print)
    {
        printArray(sobel, height, width, new_dir_name, "sobel.jpg");
    }

    for (int i = 1; i < height - 1; i++)
    {
        for (int j = 1; j < width - 1; j++)
        {
            if (sobel[i][j] > threshold_limit)
            {
                sobel[i][j] = 0;
            }
            else
            {
                sobel[i][j] = 255;
            }
        }
    }
    if (print)
    {
        printArray(sobel, height, width, new_dir_name, "binary_threshold.jpg");
    }
    // erosion

    int vertical_size = height / 15;

    for (int i = 0; i < height; i++)
    {
        for (int j = vertical_size / 2; j < width - vertical_size / 2; j++)
        {
            int counter = 0;
            for (int k = -vertical_size / 2; k < vertical_size / 2; k++)
            {

                if (
                    sobel[i][j + k] == 0)
                {
                    counter++;
                }
            }
            if (counter != 0)
            {
                erosion[i][j] = 0;
            }
            else
            {
                erosion[i][j] = 255;
            }
        }
    }
    if (print)
    {
        printArray(erosion, height, width, new_dir_name, "first_erosion.jpg");
    }
    // dilation 1
    for (int i = 0; i < height; i++)
    {
        for (int j = vertical_size / 2; j < width - vertical_size / 2; j++)
        {
            int counter = 0;
            for (int k = -vertical_size / 2; k < vertical_size / 2; k++)
            {

                if (
                    erosion[i][j + k] == 0)
                {
                    counter++;
                }
            }
            if (counter == vertical_size)
            {
                dilation[i][j] = 0;
            }
            else
            {
                dilation[i][j] = 255;
            }
        }
    }
    if (print)
    {
        printArray(dilation, height, width, new_dir_name, "first_dilation.jpg");
    }
    // dilation 2
    for (int i = 0; i < height; i++)
    {
        for (int j = vertical_size / 2; j < width - vertical_size / 2; j++)
        {
            int counter = 0;
            for (int k = -vertical_size / 2; k < vertical_size / 2; k++)
            {

                if (
                    dilation[i][j + k] == 0)
                {
                    counter++;
                }
            }
            if (counter == vertical_size)
            {
                erosion[i][j] = 0;
            }
            else
            {
                erosion[i][j] = 255;
            }
        }
    }
    if (print)
    {
        printArray(erosion, height, width, new_dir_name, "second_dilation.jpg");
    }
    // dilation 3
    int size = 7;
    for (int i = size; i < height - size; i++)
    {
        for (int j = size; j < width - size; j++)
        {

            if (
                erosion[i - 3][j] == 0 &&
                erosion[i - 2][j] == 0 &&
                erosion[i - 1][j] == 0 &&
                erosion[i][j - 3] == 0 &&
                erosion[i][j - 2] == 0 &&
                erosion[i][j - 1] == 0 &&
                erosion[i][j + 1] == 0 &&
                erosion[i][j + 2] == 0 &&
                erosion[i][j + 3] == 0 &&
                erosion[i + 1][j] == 0 &&
                erosion[i + 2][j] == 0 &&
                erosion[i + 3][j] == 0)
            {
                dilation[i][j] = 0;
            }
            else
            {
                dilation[i][j] = 255;
            }
        }
    }
    if (print)
    {
        printArray(dilation, height, width, new_dir_name, "morph_cross.jpg");
    }

    cropborder(dilation, vertical_size, height, width);
    if (print)
    {
        printArray(dilation, height, width, new_dir_name, "crop.jpg");
    }
    int finalx = 0;
    int finaly = 0;
    for (int m = 5; m < 50; m = m + 1)
    {

        int box_size = m;

        for (int i = 0; i < height - box_size; i = i + box_size)
        {

            for (int j = 0; j < width - box_size; j = j + box_size)
            {
                int counter = 0;
                for (int k = 0; k < box_size; k++)
                {
                    for (int l = 0; l < box_size; l++)
                    {
                        if (dilation[i + k][j + l] == 0)
                        {
                            counter++;
                        }
                    }
                }
                if (counter == box_size * box_size)
                {

                    finalx = i;
                    finaly = j;
                }
            }
        }
    }

    floodFill(dilation, finalx, finaly, 125, height, width);

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (dilation[i][j] == 125)
            {
                dilation[i][j] = 0;
            }
            else
            {
                dilation[i][j] = 255;
            }
        }
    }
    if (print)
    {
        printArray(dilation, height, width, new_dir_name, "flood_fill.jpg");
    }

    int ulx = 0;
    int uly = 0;
    int llx = 0;
    int lly = 0;
    int urx = 0;
    int ury = 0;
    int lrx = 0;
    int lry = 0;
    int cx = 0;
    int cy = 0;
    // counter=0;
    for (int i = height - 1; i >= 0; i--)
    {
        for (int j = width - 1; j >= 0; j--)
        {
            if (dilation[i][j] == 0)
            {
                ury = j;
                urx = i;
                break;
            }
        }
    }

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (dilation[i][j] == 0)
            {
                lly = j;
                llx = i;
                break;
            }
        }
    }

    for (int i = height - 1; i >= 0; i--)
    {
        for (int j = 0; j < width; j++)
        {
            if (dilation[i][j] == 0 && abs(i - urx) > 25)
            {
                uly = j;
                ulx = i;

                break;
            }
        }
    }

    for (int i = 0; i < height; i++)
    {
        for (int j = width - 1; j >= 0; j--)
        {

            if (dilation[i][j] == 0 && abs(j - uly) > 25)
            {
                lry = j;
                lrx = i;
                break;
            }
        }
    }

    cx = ((urx + lrx) / 2 - (ulx + llx) / 2) / 2 + (ulx + llx) / 2;
    cy = ((lly + lry) / 2 - (uly + ury) / 2) / 2 + (uly + ury) / 2;

    int cx_left = cx - 25, cx_right = cx + 25, cy_upper = cy + 25, cy_lower = cy - 25;

    for (int i = 0; i < 25; i++)
    {

        if (dilation[cx + i][cy] == 255)
        {
            cx_right = cx + i - 1;
            break;
        }
    }
    for (int i = 0; i < 25; i++)
    {

        if (dilation[cx - i][cy] == 255)
        {
            cx_left = cx - i + 1;
            break;
        }
    }
    for (int j = 0; j < 25; j++)
    {

        if (dilation[cx][cy - j] == 255)
        {
            cy_upper = cy - j + 1;
            break;
        }
    }

    for (int j = 0; j < 25; j++)
    {

        if (dilation[cx][cy + j] == 255)
        {
            cy_lower = cy + j - 1;
            break;
        }
    }
    double sum_intensity = 0;
    int counter = 0;
    double average_intensity = 0.0;
    double sum_background_intensity = 0;
    double average_background_intensity = 0.0;
    double sum_barcodelines_intensity = 0;

    for (int i = cx_left; i <= cx_right; i++)
    {
        for (int j = cy_lower; j <= cy_upper; j++)
        {
            sum_intensity = sum_intensity + matrix[i][j];
            counter++;
        }
    }
    average_intensity = sum_intensity / counter;
    if (print)
    {
        printf("\n average intensity %f", average_intensity);
    }

    counter = 0;
    for (int i = cx_left; i <= cx_right; i++)
    {
        for (int j = cy_lower; j <= cy_upper; j++)
        {
            if (matrix[i][j] >= average_intensity)
            {
                sum_background_intensity = sum_background_intensity + matrix[i][j];
                counter++;
            }
        }
    }
    average_background_intensity = sum_background_intensity / counter;
    if (print)
    {
        printf("\n average background intensity %f", average_background_intensity);
    }
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (matrix[i][j] > average_intensity)
            {
                //  printf("\n matrix[i][j] %d",matrix[i][j]);
                erosion[i][j] = 255;
            }
            else
            {
                erosion[i][j] = 0;
                ;
            }
        }
    }

    int cx_copy = cx;
    int cy_copy = cy;
    for (int i = cy_copy - 25; i < cy_copy + 25; i++)
    {
        for (int j = cx_copy - 25; j < cx_copy + 25; j++)
        {
            if (erosion[i][j] == 0)
            {
                cx = j;
                cy = i;
                break;
            }
        }
    }

    if (print)
    {

        printf("\n height %d cy %d width %d cx %d erosion[cy][cx] %d", height, cy, width, cx, erosion[cy][cx]);
    }

    vector<point> arr;

    // point dummy;
    // dummy.x=cx;
    // dummy.y=cy;
    // arr.push_back(dummy);
    // arr.push_back(dummy);
    // switch initial conditions cx for cy etc
    int i = cx;
    // for( int j= cy ;j< width ;j++)

    int j = cy;

    while (j < width)
    {

        int minimum_y = cx;
        int minimum_x = cy;

        j++;
        i = cx;

        int blank_space_counter = 0;
        for (int k = 0; k < 15; k++)
        {

            if (erosion[i][j + k] == 255)
            {

                blank_space_counter++;
            }
        }
        if (blank_space_counter == 15)
        {

            break;
        }

        while (i < height)
        {

            if ((erosion[i][j] == 0 || erosion[i][j + 1] == 0))
            {
                if (i >= minimum_y)
                {
                    point temp;
                    minimum_y = i;
                    minimum_x = j;

                    temp.x = minimum_x;
                    temp.y = minimum_y;

                    arr.push_back(temp);
                }
                i++;
            }
            else
            {
                break;
            }
        }
    }

    j = cy;

    while (j < width)
    {
        int minimum_y = cx;
        int minimum_x = cy;
        j++;
        i = cx;

        int blank_space_counter = 0;
        for (int k = 0; k < 15; k++)
        {

            if (erosion[i][j + k] == 255)
            {

                blank_space_counter++;
            }
        }
        if (blank_space_counter == 15)
        {

            break;
        }

        while (i > 0)
        {

            if (

                (erosion[i][j] == 0 || erosion[i][j + 1] == 0))
            {
                if (i <= minimum_y)
                {
                    point temp;
                    minimum_y = i;
                    minimum_x = j;

                    temp.x = minimum_x;
                    temp.y = minimum_y;

                    arr.push_back(temp);
                }
                i--;
            }
            else
            {

                break;
            }
        }
    }

    j = cy;

    while (j > 0)
    {
        int minimum_y = cx;
        int minimum_x = cy;
        j--;
        i = cx;

        int blank_space_counter = 0;
        for (int k = 0; k < 15; k++)
        {

            if (erosion[i][j - k] == 255)
            {

                blank_space_counter++;
            }
        }
        if (blank_space_counter == 15)
        {

            break;
        }

        while (i < height)
        {

            if (

                (erosion[i][j] == 0 || erosion[i][j + 1] == 0))
            {
                if (i >= minimum_y)
                {
                    point temp;
                    minimum_y = i;
                    minimum_x = j;

                    temp.x = minimum_x;
                    temp.y = minimum_y;

                    arr.push_back(temp);
                }
                i++;
            }
            else
            {

                break;
            }
        }
    }

    j = cy;

    while (j > 0)
    {

        int minimum_y = cx;
        int minimum_x = cy;

        j--;
        i = cx;

        int blank_space_counter = 0;
        for (int k = 0; k < 15; k++)
        {

            if (erosion[i][j - k] == 255)
            {

                blank_space_counter++;
            }
        }
        if (blank_space_counter == 15)
        {

            break;
        }

        while (i > 0)
        {

            if ((erosion[i][j] == 0 || erosion[i][j + 1] == 0))
            {
                if (i <= minimum_y)
                {
                    point temp;
                    minimum_y = i;
                    minimum_x = j;

                    temp.x = minimum_x;
                    temp.y = minimum_y;

                    arr.push_back(temp);
                }
                i--;
            }
            else
            {

                break;
            }
        }
    }

    for (int i = 0; i < arr.size() - 1; i++)
    {
        Line(erosion, arr[i].y, arr[i].x, arr[i].y, arr[i].x, 125);
    }
    if (print)
    {
        printArray(erosion, height, width, new_dir_name, "average_intensity.jpg");
    }

    int max_y = 0;
    int max_x = 0;
    int min_y = 10000;
    int min_x = 10000;

    for (int i = 0; i < arr.size() - 1; i++)
    {
        if (arr[i].y > max_y)
        {
            max_y = arr[i].y;
        }
        if (arr[i].x > max_x)
        {
            max_x = arr[i].x;
        }
        if (arr[i].y < min_y)
        {
            min_y = arr[i].y;
        }
        if (arr[i].x < min_x)
        {
            min_x = arr[i].x;
        }
    }

    int correction_x = 15;
    int colour1 = 0;
    int colour2 = 0;

    for (int i = 15; i <= 50; i++)
    {
        colour1 = erosion[((max_y - min_y) / 2 + min_y)][min_x - i];
        colour2 = erosion[((max_y - min_y) / 2 + min_y)][min_x - i - 1];
        if (colour2 != colour1)
        {
            correction_x = i - 1;
            break;
        }
    }
    if (print)
    {

        printf("\ncorrection_x %d colour1 %d colour2 %d", correction_x, colour1, colour2);
    }

    min_x = min_x - correction_x;
    max_x = max_x + 10;

    int h = max_y - min_y;
    int w = max_x - min_x;

    printf("\n croped size %d, %d, %d, %d", min_x, max_x, min_y, max_y);
    byte **crop = new byte *[(max_y - min_y)];
    for (int i = 0; i < (max_y - min_y); i++)
    {
        crop[i] = new byte[(max_x - min_x)];
    }

    for (int i = min_y; i < max_y; i++)
    {
        for (int j = min_x; j < max_x; j++)
        {
            crop[i - min_y][j - min_x] = img[i][j];
        }
    }

    byte **ret = new byte *[h];
    for (int i = 0; i < h; i++)
    {
        ret[i] = new byte[w];
    }
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
            ret[i][j] = crop[i][j];
        }
    }

    he = h;
    wi = w;
    arr2 = ret;
    if (print)
    {
        printArray(crop, h, w, new_dir_name, "cropped.jpg");
    }

    for (int i = 0; i < height; ++i)
    {

        delete[] erosion[i];
        delete[] dilation[i];
        delete[] sobel[i];
        delete[] label[i];
        delete[] matrix[i];
    }

    return cv::Rect(min_x, min_y, w, h);
    //			  		return new_dir_name+"cropped.jpg";
}

/**
 * @brief Get the Scanline object
 * 
 * @param matrix 假设矩阵里包含了条码信息
 * @param height 
 * @param width 
 * @param print 调试输出
 * @param k 扫描位置因子
 * @return vector<long double> 
 */
vector<long double> getScanline(byte **matrix, int height, int width, bool print, double k, int &left)
{

    vector<long double> scanline;

    for (int j = 0; j < width; j++)
    {

        int index = height * k;
        long double intensity = matrix[index][j];

        scanline.push_back(intensity);
        if (print)
        {
            // 调试扫描行信息
            printf("\n( %d , %f )", j, (double)intensity);
        }
    }

    double sum = 0;
    int counter = 0;
    ;
    long double average = 0.0;
    for (unsigned int i = 0; i < scanline.size(); i++)
    {

        // if(scanline[i]>0.5)
        //{
        sum = sum + scanline[i];
        counter++;
        //}
    }
    average = sum / counter;
    // 为什么做二值化处理
    for (unsigned int i = 0; i < scanline.size(); i++)
    {
        if (scanline[i] > 255)
        {
            scanline[i] = 255;
        }

        if (scanline[i] < 0)
        {
            scanline[i] = 0;
        }

        if (scanline[i] > average)
        {
            scanline[i] = 1.0;
        }
        else if (scanline[i] <= average)
        {
            scanline[i] = 0.0;
        }
    }
    int index_of_first_black_pixel = 0;
    int index_of_last_black_pixel = 639;  // 因为图像是 640 x 480 的

    // 最右边的条码线
    for (int j = scanline.size() - 1; j >= 0; j--)
    {
        if (scanline[j] < 0.5)
        {

            index_of_last_black_pixel = j;
            break;
        }
    }
    scanline.erase(scanline.begin() + index_of_last_black_pixel + 1, scanline.end());

    // 最左边的条码线
    for (int j = 0; j < scanline.size(); j++)
    {
        if (scanline[j] < 0.5)
        {
            index_of_first_black_pixel = j;

            break;
        }
    }
    scanline.erase(scanline.begin(), scanline.begin() + index_of_first_black_pixel - 1);

    if (print)
    {
        printf("\nindex_of_first_black_pixel %d", index_of_first_black_pixel);

        printf("\nindex_of_last_black_pixel %d", index_of_last_black_pixel);
        printf("\nscanline length %llu", scanline.size());
    }
    left = index_of_first_black_pixel;
    return scanline;
}

vector<vector<int>> bars(vector<long double> scanline, bool print)
{
    double sum_intensity = 0;
    double average_intensity = 0;
    for (int i = 0; i < scanline.size(); i++)
    {
        sum_intensity = sum_intensity + scanline[i];
    }
    average_intensity = sum_intensity / scanline.size();
    if (print)
    {
        printf("\naverage_intensity %f", average_intensity);
    }

    // 每个像素点标记
    vector<int> binarizedPoints2;
    for (unsigned int i = 0; i < scanline.size(); i++)
    {
        if (scanline[i] >= average_intensity)
        {
            binarizedPoints2.push_back(0);
        }
        else
        {
            binarizedPoints2.push_back(1);
        }
    }

    if (print)
    {
        for (unsigned int i = 0; i < binarizedPoints2.size(); i++)
        {
            printf("\n(%d,%d)", i, binarizedPoints2[i]);
        }
        printf("\n binarizedPoints2 length %llu", binarizedPoints2.size());
    }
    vector<vector<int>> bars;
    int l = binarizedPoints2.size();

    int i = 0;
    vector<int> templist0;
    vector<int> templist1;
    vector<int> templist;

    // 按照标记对线条分组，最宽不能超过 60 像素
    while (i <= l)
    {
    label0:
        if (binarizedPoints2[i] == 0)
        {
            templist0.push_back(0);

            i++;
            goto label0;
        }
        else
        {
            bars.push_back(templist0);
            templist0.clear();
        }

        if (bars.size() >= 60)
        {
            break;
        }
    label1:

        if (binarizedPoints2[i] == 1)
        {
            templist1.push_back(1);

            i++;
            goto label1;
        }
        else
        {
            bars.push_back(templist1);
            templist1.clear();
        }

        if (bars.size() >= 60)
        {
            break;
        }
    }
    if (print)
    {
        printf("\nnumber of  bars   %llu", bars.size());
    }
    for (unsigned int i = l; i > 0; i--)
    {
        if (binarizedPoints2[l] == binarizedPoints2[i])
        {
            templist.push_back(binarizedPoints2[l]);
        }
        else
        {
            break;
        }
    }
    bars.push_back(templist);

    return bars;
}

vector<int> leftPadding(vector<vector<int>> bars)
{

    return bars.at(0);
}

vector<int> rightPadding(vector<vector<int>> bars)
{

    return bars.at(60);
}

vector<vector<int>> leftGuard(vector<vector<int>> bars)
{

    vector<vector<int>> leftGuard;

    for (unsigned int i = 0; i < 3; i++)
    {
        leftGuard.push_back(bars.at(i + 1));
    }

    return leftGuard;
}

vector<vector<int>> rightGuard(vector<vector<int>> bars)
{
    vector<vector<int>> rightGuard;

    for (unsigned int i = 0; i < 3; i++)
    {
        rightGuard.push_back(bars.at(i + 57));
    }

    return rightGuard;
}

vector<vector<int>> midGuard(vector<vector<int>> bars)
{
    vector<vector<int>> midGuard;

    for (unsigned int i = 0; i < 5; i++)
    {
        midGuard.push_back(bars.at(i + 28));
    }

    return midGuard;
}
vector<vector<vector<int>>> lhsAndRhsBars(vector<vector<int>> bars)
{
    vector<vector<vector<int>>> lhsAndRhsBars;
    vector<vector<int>> temp;

    for (unsigned int j = 0; j < 6; j++)
    {
        temp.clear();
        for (unsigned int i = 0; i < 4; i++)
        {
            temp.push_back(bars.at(i + j * 4 + 4));
        }
        lhsAndRhsBars.push_back(temp);
    }

    for (unsigned int j = 0; j < 6; j++)
    {
        temp.clear();
        for (unsigned int i = 0; i < 4; i++)
        {
            temp.push_back(bars.at(i + j * 4 + 33));
        }
        lhsAndRhsBars.push_back(temp);
    }
    return lhsAndRhsBars;
}

vector<vector<vector<int>>> lhsBars(vector<vector<int>> bars)
{
    vector<vector<vector<int>>> lhsBars;
    vector<vector<int>> temp;

    for (unsigned int j = 0; j < 6; j++)
    {
        temp.clear();
        for (unsigned int i = 0; i < 4; i++)
        {
            temp.push_back(bars.at(i + j * 4 + 4));
        }
        lhsBars.push_back(temp);
    }

    return lhsBars;
}

vector<vector<vector<int>>> rhsBars(vector<vector<int>> bars)
{
    vector<vector<vector<int>>> rhsBars;
    vector<vector<int>> temp;

    for (unsigned int j = 0; j < 6; j++)
    {
        temp.clear();
        for (unsigned int i = 0; i < 4; i++)
        {
            temp.push_back(bars.at(i + j * 4 + 33));
        }
        rhsBars.push_back(temp);
    }
    return rhsBars;
}

double averageBarLength1(vector<vector<int>> bars)
{
    vector<vector<int>> lGuard = leftGuard(bars);
    vector<vector<int>> rGuard = rightGuard(bars);
    vector<vector<int>> mGuard = midGuard(bars);
    double avgBarLength1 = 0;

    avgBarLength1 = (lGuard[0].size() + lGuard[2].size() + rGuard[0].size() + rGuard[2].size() +
                     mGuard[1].size() + mGuard[3].size()) /
                    6;

    return avgBarLength1;
}

double averageBarLength0(vector<vector<int>> bars)
{
    vector<vector<int>> lGuard = leftGuard(bars);
    vector<vector<int>> rGuard = rightGuard(bars);
    vector<vector<int>> mGuard = midGuard(bars);
    double avgBarLength0 = 0;

    avgBarLength0 = (lGuard[1].size() + rGuard[1].size() + rGuard[2].size() +
                     mGuard[0].size() + mGuard[2].size() + mGuard[4].size()) /
                    5;

    return avgBarLength0;
}
long double maxDist_lhs(int j, long double avbl1, long double avbl0, vector<vector<vector<int>>> bars_lhs)
{
    long double avgBarLength0 = avbl0;
    long double avgBarLength1 = avbl1;
    vector<vector<vector<int>>> s = bars_lhs;

    vector<long double> r0L = {avgBarLength0 * 3, avgBarLength1 * 2, avgBarLength0,
                               avgBarLength1};
    vector<long double> r0G = {avgBarLength0, avgBarLength1, avgBarLength0 * 2,
                               avgBarLength1 * 3};

    vector<long double> r1L = {avgBarLength0 * 2, avgBarLength1 * 2, avgBarLength0 * 2,
                               avgBarLength1};
    vector<long double> r1G = {avgBarLength0, avgBarLength1 * 2, avgBarLength0 * 2,
                               avgBarLength1 * 2};

    vector<long double> r2L = {avgBarLength0 * 2, avgBarLength1, avgBarLength0 * 2,
                               avgBarLength1 * 2};
    vector<long double> r2G = {avgBarLength0 * 2, avgBarLength1 * 2, avgBarLength0,
                               avgBarLength1 * 2};

    vector<long double> r3L = {avgBarLength0, avgBarLength1 * 3, avgBarLength0, avgBarLength1};
    vector<long double> r3G = {avgBarLength0, avgBarLength1, avgBarLength0 * 4, avgBarLength1};

    vector<long double> r4L = {avgBarLength0, avgBarLength1, avgBarLength0 * 3,
                               avgBarLength1 * 2};
    vector<long double> r4G = {avgBarLength0 * 2, avgBarLength1 * 3, avgBarLength0,
                               avgBarLength1};

    vector<long double> r5L = {avgBarLength0, avgBarLength1 * 2, avgBarLength0 * 3,
                               avgBarLength1};
    vector<long double> r5G = {avgBarLength0, avgBarLength1 * 3, avgBarLength0 * 2,
                               avgBarLength1};

    vector<long double> r6L = {avgBarLength0, avgBarLength1, avgBarLength0, avgBarLength1 * 4};
    vector<long double> r6G = {avgBarLength0 * 4, avgBarLength1, avgBarLength0, avgBarLength1};

    vector<long double> r7L = {avgBarLength0, avgBarLength1 * 3, avgBarLength0,
                               avgBarLength1 * 2};
    vector<long double> r7G = {avgBarLength0 * 2, avgBarLength1, avgBarLength0 * 3,
                               avgBarLength1};

    vector<long double> r8L = {avgBarLength0, avgBarLength1 * 2, avgBarLength0,
                               avgBarLength1 * 3};
    vector<long double> r8G = {avgBarLength0 * 3, avgBarLength1, avgBarLength0 * 2,
                               avgBarLength1};

    vector<long double> r9L = {avgBarLength0 * 3, avgBarLength1, avgBarLength0,
                               avgBarLength1 * 2};
    vector<long double> r9G = {avgBarLength0 * 2, avgBarLength1, avgBarLength0,
                               avgBarLength1 * 3};

    vector<vector<long double>> r = {r0L, r0G, r1L, r1G, r2L, r2G, r3L, r3G,
                                     r4L, r4G, r5L, r5G, r6L, r6G, r7L, r7G, r8L,
                                     r8G, r9L, r9G};

    long double total = 0;

    for (unsigned int i = 0; i < 20; i++)
    {
        long double max = 0;
        max = pow(r.at(i).at(0) - s.at(j).at(0).size(), 2) +
              pow(r.at(i).at(1) - s.at(j).at(1).size(), 2) +
              pow(r.at(i).at(2) - s.at(j).at(2).size(), 2) +
              pow(r.at(i).at(3) - s.at(j).at(3).size(), 2);

        if (max >= total)
        {
            total = max;
        }
    }
    return total;
}

long double maxDist_rhs(int j, long double avbl1, long double avbl0, vector<vector<vector<int>>> bars_rhs)
{
    long double avgBarLength0 = avbl0;
    long double avgBarLength1 = avbl1;
    vector<vector<vector<int>>> s = bars_rhs;

    vector<long double> r0R = {avgBarLength1 * 3, avgBarLength0 * 2, avgBarLength1,
                               avgBarLength0};

    vector<long double> r1R = {avgBarLength1 * 2, avgBarLength0 * 2, avgBarLength1 * 2,
                               avgBarLength0};

    vector<long double> r2R = {avgBarLength1 * 2, avgBarLength0, avgBarLength1 * 2,
                               avgBarLength0 * 2};

    vector<long double> r3R = {avgBarLength1, avgBarLength0 * 3, avgBarLength1, avgBarLength0};

    vector<long double> r4R = {avgBarLength1, avgBarLength0, avgBarLength1 * 3,
                               avgBarLength0 * 2};

    vector<long double> r5R = {avgBarLength1, avgBarLength0 * 2, avgBarLength1 * 3,
                               avgBarLength0};

    vector<long double> r6R = {avgBarLength1, avgBarLength0, avgBarLength1, avgBarLength0 * 4};

    vector<long double> r7R = {avgBarLength1, avgBarLength0 * 3, avgBarLength1,
                               avgBarLength0 * 2};

    vector<long double> r8R = {avgBarLength1, avgBarLength0 * 2, avgBarLength1,
                               avgBarLength0 * 3};

    vector<long double> r9R = {avgBarLength1 * 3, avgBarLength0, avgBarLength1,
                               avgBarLength0 * 2};

    vector<vector<long double>> r = {r0R, r1R, r2R, r3R,
                                     r4R, r5R, r6R, r7R, r8R, r9R};

    long double total = 0;

    for (unsigned int i = 0; i < 10; i++)
    {
        long double max = 0;
        max = pow(r.at(i).at(0) - s.at(j).at(0).size(), 2) +
              pow(r.at(i).at(1) - s.at(j).at(1).size(), 2) +
              pow(r.at(i).at(2) - s.at(j).at(2).size(), 2) +
              pow(r.at(i).at(3) - s.at(j).at(3).size(), 2);

        if (max >= total)
        {
            total = max;
        }
    }
    return total;
}

long double pdash_lhs(int i, int j, long double avbl1, long double avbl0, vector<vector<vector<int>>> bars_lhs)
{
    long double avgBarLength0 = avbl0;
    long double avgBarLength1 = avbl1;
    vector<vector<vector<int>>> s = bars_lhs;

    vector<long double> r0L = {avgBarLength0 * 3, avgBarLength1 * 2, avgBarLength0,
                               avgBarLength1};
    vector<long double> r0G = {avgBarLength0, avgBarLength1, avgBarLength0 * 2,
                               avgBarLength1 * 3};

    vector<long double> r1L = {avgBarLength0 * 2, avgBarLength1 * 2, avgBarLength0 * 2,
                               avgBarLength1};
    vector<long double> r1G = {avgBarLength0, avgBarLength1 * 2, avgBarLength0 * 2,
                               avgBarLength1 * 2};

    vector<long double> r2L = {avgBarLength0 * 2, avgBarLength1, avgBarLength0 * 2,
                               avgBarLength1 * 2};
    vector<long double> r2G = {avgBarLength0 * 2, avgBarLength1 * 2, avgBarLength0,
                               avgBarLength1 * 2};

    vector<long double> r3L = {avgBarLength0, avgBarLength1 * 3, avgBarLength0, avgBarLength1};
    vector<long double> r3G = {avgBarLength0, avgBarLength1, avgBarLength0 * 4, avgBarLength1};

    vector<long double> r4L = {avgBarLength0, avgBarLength1, avgBarLength0 * 3,
                               avgBarLength1 * 2};
    vector<long double> r4G = {avgBarLength0 * 2, avgBarLength1 * 3, avgBarLength0,
                               avgBarLength1};

    vector<long double> r5L = {avgBarLength0, avgBarLength1 * 2, avgBarLength0 * 3,
                               avgBarLength1};
    vector<long double> r5G = {avgBarLength0, avgBarLength1 * 3, avgBarLength0 * 2,
                               avgBarLength1};

    vector<long double> r6L = {avgBarLength0, avgBarLength1, avgBarLength0, avgBarLength1 * 4};
    vector<long double> r6G = {avgBarLength0 * 4, avgBarLength1, avgBarLength0, avgBarLength1};

    vector<long double> r7L = {avgBarLength0, avgBarLength1 * 3, avgBarLength0,
                               avgBarLength1 * 2};
    vector<long double> r7G = {avgBarLength0 * 2, avgBarLength1, avgBarLength0 * 3,
                               avgBarLength1};

    vector<long double> r8L = {avgBarLength0, avgBarLength1 * 2, avgBarLength0,
                               avgBarLength1 * 3};
    vector<long double> r8G = {avgBarLength0 * 3, avgBarLength1, avgBarLength0 * 2,
                               avgBarLength1};

    vector<long double> r9L = {avgBarLength0 * 3, avgBarLength1, avgBarLength0,
                               avgBarLength1 * 2};
    vector<long double> r9G = {avgBarLength0 * 2, avgBarLength1, avgBarLength0,
                               avgBarLength1 * 3};

    vector<vector<long double>> r = {r0L, r0G, r1L, r1G, r2L, r2G, r3L, r3G,
                                     r4L, r4G, r5L, r5G, r6L, r6G, r7L, r7G, r8L,
                                     r8G, r9L, r9G};

    long double sum = 0;
    sum = pow(r.at(i).at(0) - s.at(j).at(0).size(), 2) +
          pow(r.at(i).at(1) - s.at(j).at(1).size(), 2) +
          pow(r.at(i).at(2) - s.at(j).at(2).size(), 2) +
          pow(r.at(i).at(3) - s.at(j).at(3).size(), 2);

    long double total = 0.0;
    total = 1 - (sum / maxDist_lhs(j, avgBarLength1, avgBarLength0, s));
    return total;
}

long double pdash_rhs(int i, int j, long double avbl1, long double avbl0, vector<vector<vector<int>>> bars_rhs)
{
    long double avgBarLength0 = avbl0;
    long double avgBarLength1 = avbl1;
    vector<vector<vector<int>>> s = bars_rhs;

    vector<long double> r0R = {avgBarLength1 * 3, avgBarLength0 * 2, avgBarLength1,
                               avgBarLength0};

    vector<long double> r1R = {avgBarLength1 * 2, avgBarLength0 * 2, avgBarLength1 * 2,
                               avgBarLength0};

    vector<long double> r2R = {avgBarLength1 * 2, avgBarLength0, avgBarLength1 * 2,
                               avgBarLength0 * 2};

    vector<long double> r3R = {avgBarLength1, avgBarLength0 * 3, avgBarLength1, avgBarLength0};

    vector<long double> r4R = {avgBarLength1, avgBarLength0, avgBarLength1 * 3,
                               avgBarLength0 * 2};

    vector<long double> r5R = {avgBarLength1, avgBarLength0 * 2, avgBarLength1 * 3,
                               avgBarLength0};

    vector<long double> r6R = {avgBarLength1, avgBarLength0, avgBarLength1, avgBarLength0 * 4};

    vector<long double> r7R = {avgBarLength1, avgBarLength0 * 3, avgBarLength1,
                               avgBarLength0 * 2};

    vector<long double> r8R = {avgBarLength1, avgBarLength0 * 2, avgBarLength1,
                               avgBarLength0 * 3};

    vector<long double> r9R = {avgBarLength1 * 3, avgBarLength0, avgBarLength1,
                               avgBarLength0 * 2};

    vector<vector<long double>> r = {r0R, r1R, r2R, r3R,
                                     r4R, r5R, r6R, r7R, r8R, r9R};

    long double sum = 0;
    sum = pow(r.at(i).at(0) - s.at(j).at(0).size(), 2) +
          pow(r.at(i).at(1) - s.at(j).at(1).size(), 2) +
          pow(r.at(i).at(2) - s.at(j).at(2).size(), 2) +
          pow(r.at(i).at(3) - s.at(j).at(3).size(), 2);

    long double total = 0.0;
    total = 1 - (sum / maxDist_rhs(j, avgBarLength1, avgBarLength0, s));
    return total;
}
long double p_lhs(int i, int j, double avbl1, double avbl0, vector<vector<vector<int>>> bars_lhs)
{
    long double avgBarLength0 = avbl0;
    long double avgBarLength1 = avbl1;
    vector<vector<vector<int>>> s = bars_lhs;

    vector<long double> r0L = {avgBarLength0 * 3, avgBarLength1 * 2, avgBarLength0,
                               avgBarLength1};
    vector<long double> r0G = {avgBarLength0, avgBarLength1, avgBarLength0 * 2,
                               avgBarLength1 * 3};

    vector<long double> r1L = {avgBarLength0 * 2, avgBarLength1 * 2, avgBarLength0 * 2,
                               avgBarLength1};
    vector<long double> r1G = {avgBarLength0, avgBarLength1 * 2, avgBarLength0 * 2,
                               avgBarLength1 * 2};

    vector<long double> r2L = {avgBarLength0 * 2, avgBarLength1, avgBarLength0 * 2,
                               avgBarLength1 * 2};
    vector<long double> r2G = {avgBarLength0 * 2, avgBarLength1 * 2, avgBarLength0,
                               avgBarLength1 * 2};

    vector<long double> r3L = {avgBarLength0, avgBarLength1 * 3, avgBarLength0, avgBarLength1};
    vector<long double> r3G = {avgBarLength0, avgBarLength1, avgBarLength0 * 4, avgBarLength1};

    vector<long double> r4L = {avgBarLength0, avgBarLength1, avgBarLength0 * 3,
                               avgBarLength1 * 2};
    vector<long double> r4G = {avgBarLength0 * 2, avgBarLength1 * 3, avgBarLength0,
                               avgBarLength1};

    vector<long double> r5L = {avgBarLength0, avgBarLength1 * 2, avgBarLength0 * 3,
                               avgBarLength1};
    vector<long double> r5G = {avgBarLength0, avgBarLength1 * 3, avgBarLength0 * 2,
                               avgBarLength1};

    vector<long double> r6L = {avgBarLength0, avgBarLength1, avgBarLength0, avgBarLength1 * 4};
    vector<long double> r6G = {avgBarLength0 * 4, avgBarLength1, avgBarLength0, avgBarLength1};

    vector<long double> r7L = {avgBarLength0, avgBarLength1 * 3, avgBarLength0,
                               avgBarLength1 * 2};
    vector<long double> r7G = {avgBarLength0 * 2, avgBarLength1, avgBarLength0 * 3,
                               avgBarLength1};

    vector<long double> r8L = {avgBarLength0, avgBarLength1 * 2, avgBarLength0,
                               avgBarLength1 * 3};
    vector<long double> r8G = {avgBarLength0 * 3, avgBarLength1, avgBarLength0 * 2,
                               avgBarLength1};

    vector<long double> r9L = {avgBarLength0 * 3, avgBarLength1, avgBarLength0,
                               avgBarLength1 * 2};
    vector<long double> r9G = {avgBarLength0 * 2, avgBarLength1, avgBarLength0,
                               avgBarLength1 * 3};

    vector<vector<long double>> r = {r0L, r0G, r1L, r1G, r2L, r2G, r3L, r3G,
                                     r4L, r4G, r5L, r5G, r6L, r6G, r7L, r7G, r8L,
                                     r8G, r9L, r9G};

    long double total = 0;
    for (unsigned int k = 0; k < 20; k++)
    {
        total = total + pdash_lhs(k, j, avgBarLength1, avgBarLength0, s);
    }

    long double ret = pdash_lhs(i, j, avgBarLength1, avgBarLength0, s) / total;
    return ret;
}

long double p_rhs(int i, int j, double avbl1, double avbl0, vector<vector<vector<int>>> bars_rhs)
{
    long double avgBarLength0 = avbl0;
    long double avgBarLength1 = avbl1;
    vector<vector<vector<int>>> s = bars_rhs;

    vector<long double> r0R = {avgBarLength1 * 3, avgBarLength0 * 2, avgBarLength1,
                               avgBarLength0};

    vector<long double> r1R = {avgBarLength1 * 2, avgBarLength0 * 2, avgBarLength1 * 2,
                               avgBarLength0};

    vector<long double> r2R = {avgBarLength1 * 2, avgBarLength0, avgBarLength1 * 2,
                               avgBarLength0 * 2};

    vector<long double> r3R = {avgBarLength1, avgBarLength0 * 3, avgBarLength1, avgBarLength0};

    vector<long double> r4R = {avgBarLength1, avgBarLength0, avgBarLength1 * 3,
                               avgBarLength0 * 2};

    vector<long double> r5R = {avgBarLength1, avgBarLength0 * 2, avgBarLength1 * 3,
                               avgBarLength0};

    vector<long double> r6R = {avgBarLength1, avgBarLength0, avgBarLength1, avgBarLength0 * 4};

    vector<long double> r7R = {avgBarLength1, avgBarLength0 * 3, avgBarLength1,
                               avgBarLength0 * 2};

    vector<long double> r8R = {avgBarLength1, avgBarLength0 * 2, avgBarLength1,
                               avgBarLength0 * 3};

    vector<long double> r9R = {avgBarLength1 * 3, avgBarLength0, avgBarLength1,
                               avgBarLength0 * 2};

    vector<vector<long double>> r = {r0R, r1R, r2R, r3R,
                                     r4R, r5R, r6R, r7R, r8R, r9R};

    long double total = 0;
    for (unsigned int k = 0; k < 10; k++)
    {
        total = total + pdash_rhs(k, j, avgBarLength1, avgBarLength0, s);
    }

    long double ret = pdash_rhs(i, j, avgBarLength1, avgBarLength0, s) / total;
    return ret;
}
bool checkSum(vector<int> bc)
{
    vector<int> finalBarcode = bc;
    int evenSum = finalBarcode[10] + finalBarcode[8] + finalBarcode[6] +
                  finalBarcode[4] + finalBarcode[2] + finalBarcode[0];
    int oddSum = finalBarcode[11] + finalBarcode[9] + finalBarcode[7] +
                 finalBarcode[5] + finalBarcode[3] + finalBarcode[1];

    int checkSumDigit = (10 - (3 * oddSum + evenSum) % 10) % 10;
    if (finalBarcode[12] == checkSumDigit)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void barcode(long double st, vector<vector<int>> bars, long double avbl1, long double avbl0, long double d)
{

    long double delta = d;
    long double step = st;

    long double avgBarLength0 = avbl0;
    long double avgBarLength1 = avbl1;
    vector<vector<int>> s = bars;
    vector<vector<vector<int>>> s_lhs;
    vector<vector<vector<int>>> s_rhs;

    s_lhs = lhsBars(s);

    s_rhs = rhsBars(s);

    vector<long double> maxProbabilities_lhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesDigits_lhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesIndices_lhs = {0, 0, 0, 0, 0, 0};
    int maxProbIndex_lhs = 0;
    long double maxProb_lhs = 0;
    vector<int> maxProbabilitiesDigits2_lhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesIndices2_lhs = {0, 0, 0, 0, 0, 0};
    vector<long double> maxProbabilities2_lhs = {0, 0, 0, 0, 0, 0};
    vector<long double> finalProbList_lhs = {};
    vector<double> parity_lhs = {};

    vector<long double> maxProbabilities_rhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesDigits_rhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesIndices_rhs = {0, 0, 0, 0, 0, 0};
    int maxProbIndex_rhs = 0;
    long double maxProb_rhs = 0;
    vector<int> maxProbabilitiesDigits2_rhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesIndices2_rhs = {0, 0, 0, 0, 0, 0};
    vector<long double> maxProbabilities2_rhs = {
        0,
        0,
        0,
        0,
        0,
        0,
    };
    vector<long double> finalProbList_rhs = {};
    vector<double> parity_rhs = {};
    vector<int> finalbarCode = {};

    long double sigma = 1;
    int digitNumber = -1;

    for (long double c = (avgBarLength1 - delta); c <= (avgBarLength1 + delta); c = c + step)
    {
        for (long double e = (avgBarLength0 - delta); e <= (avgBarLength0 + delta); e = e + step)
        {

            vector<long double> r0L = {c * 3, e * 2, c, e};
            ;
            vector<long double> r0G = {c, e, c * 2, e * 3};

            vector<long double> r1L = {c * 2, e * 2, c * 2, e};
            vector<long double> r1G = {c, e * 2, c * 2, e * 2};

            vector<long double> r2L = {c * 2, e, c * 2, e * 2};
            vector<long double> r2G = {c * 2, e * 2, c, e * 2};

            vector<long double> r3L = {c, e * 3, c, e};
            vector<long double> r3G = {c, e, c * 4, e};

            vector<long double> r4L = {c, e, c * 3, e * 2};
            vector<long double> r4G = {c * 2, e * 3, c, e};

            vector<long double> r5L = {c, e * 2, c * 3, e};
            vector<long double> r5G = {c, e * 3, c * 2, e};

            vector<long double> r6L = {c, e, c, e * 4};
            vector<long double> r6G = {c * 4, e, c, e};

            vector<long double> r7L = {c, e * 3, c, e * 2};
            vector<long double> r7G = {c * 2, e, c * 3, e};

            vector<long double> r8L = {c, e * 2, c, e * 3};
            vector<long double> r8G = {c * 3, e, c * 2, e};

            vector<long double> r9L = {c * 3, e, c, e * 2};
            vector<long double> r9G = {c * 2, e, c, e * 3};

            vector<vector<long double>> r = {r0L, r0G, r1L, r1G, r2L, r2G, r3L, r3G,
                                             r4L, r4G, r5L, r5G, r6L, r6G, r7L, r7G, r8L,
                                             r8G, r9L, r9G};

            for (int a = 0; a < 6; a++)
            {

                maxProb_lhs = 0;
                maxProbIndex_lhs = 0;

                for (int b = 0; b < 20; b++)
                {

                    if (p_lhs(b, a, c, e, s_lhs) >= maxProb_lhs)
                    {
                        maxProb_lhs = p_lhs(b, a, c, e, s_lhs);

                        maxProbIndex_lhs = b;
                    }
                }

                sigma = sigma * maxProb_lhs;

                digitNumber = ((maxProbIndex_lhs - maxProbIndex_lhs % 2) / 2);

                if (maxProb_lhs >= maxProbabilities_lhs[a])
                {
                    maxProbabilities_lhs[a] = maxProb_lhs;
                    maxProbabilitiesDigits_lhs[a] = digitNumber;
                    maxProbabilitiesIndices_lhs[a] = maxProbIndex_lhs;
                }
            }

            for (unsigned int i = 0; i < 6; i++)
            {
                if (maxProbabilities_lhs[i] >= maxProbabilities2_lhs[i])
                {
                    maxProbabilities2_lhs[i] = maxProbabilities_lhs[i];
                    maxProbabilitiesDigits2_lhs[i] = maxProbabilitiesDigits_lhs[i];
                    maxProbabilitiesIndices2_lhs[i] = maxProbabilitiesIndices_lhs[i];
                }
            }
        }
    }

    sigma = 1;
    for (long double c = (avgBarLength1 - delta); c <= (avgBarLength1 + delta); c = c + step)
    {
        for (long double e = (avgBarLength0 - delta); e <= (avgBarLength0 + delta); e = e + step)
        {

            vector<long double> r0R = {e * 3, c * 2, e, c};

            vector<long double> r1R = {e * 2, c * 2, e * 2, c};

            vector<long double> r2R = {e * 2, c, e * 2, c * 2};

            vector<long double> r3R = {e, c * 3, e, c};

            vector<long double> r4R = {e, c, e * 3, c * 2};

            vector<long double> r5R = {e, c * 2, e * 3, c};

            vector<long double> r6R = {e, c, e, c * 4};

            vector<long double> r7R = {e, c * 3, e, c * 2};

            vector<long double> r8R = {e, c * 2, e, c * 3};

            vector<long double> r9R = {e * 3, c, e, c * 2};

            vector<vector<long double>> r = {r0R, r1R, r2R, r3R,
                                             r4R, r5R, r6R, r7R, r8R, r9R};

            for (int a = 0; a < 6; a++)
            {

                maxProb_rhs = 0;
                maxProbIndex_rhs = 0;

                for (int b = 0; b < 10; b++)
                {

                    if (p_rhs(b, a, c, e, s_rhs) >= maxProb_rhs)
                    {
                        maxProb_rhs = p_rhs(b, a, c, e, s_rhs);

                        maxProbIndex_rhs = b;
                    }
                }

                sigma = sigma * maxProb_rhs;

                digitNumber = maxProbIndex_rhs;

                if (maxProb_rhs >= maxProbabilities_rhs[a])
                {
                    maxProbabilities_rhs[a] = maxProb_rhs;
                    maxProbabilitiesDigits_rhs[a] = digitNumber;
                    maxProbabilitiesIndices_rhs[a] = maxProbIndex_rhs;
                }
            }

            for (unsigned int i = 0; i < 6; i++)
            {
                if (maxProbabilities_rhs[i] >= maxProbabilities2_rhs[i])
                {
                    maxProbabilities2_rhs[i] = maxProbabilities_rhs[i];
                    maxProbabilitiesDigits2_rhs[i] = maxProbabilitiesDigits_rhs[i];
                    maxProbabilitiesIndices2_rhs[i] = maxProbabilitiesIndices_rhs[i];
                }
            }

            //		 }}

            vector<int> finalBarCode = {};
            for (unsigned int i = 0; i < 6; i++)
            {
                finalBarCode.push_back(maxProbabilitiesDigits_lhs[i]);
            }
            for (unsigned int i = 0; i < 6; i++)
            {
                finalBarCode.push_back(maxProbabilitiesDigits_rhs[i]);
            }

            for (unsigned int i = 0; i < 6; i++)
            {
                if ((maxProbabilitiesIndices_lhs[i]) % 2 == 0)
                {
                    parity_lhs.push_back(1);
                }
                if ((maxProbabilitiesIndices_lhs[i]) % 2 == 1)
                {
                    parity_lhs.push_back(0);
                }
            }

            vector<vector<int>> m = {{1, 1, 1, 1, 1, 1}, {1, 1, 0, 1, 0, 0}, {1, 1, 0, 0, 1, 0}, {1, 1, 0, 0, 0, 1}, {1, 0, 1, 1, 0, 0}, {1, 0, 0, 1, 1, 0}, {1, 0, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0}, {1, 0, 1, 0, 0, 1}, {1, 0, 0, 1, 0, 1}};

            int mNumber = -1;

            for (int i = 0; i < 10; i++)
            {
                if (m.at(i).at(0) == parity_lhs.at(0) &&
                    m.at(i).at(1) == parity_lhs.at(1) &&
                    m.at(i).at(2) == parity_lhs.at(2) &&
                    m.at(i).at(3) == parity_lhs.at(3) &&
                    m.at(i).at(4) == parity_lhs.at(4) &&
                    m.at(i).at(5) == parity_lhs.at(5))
                {
                    mNumber = i;
                    break;
                }
            }
            finalBarCode.insert(finalBarCode.begin(), mNumber);

            printf("\nfinal Bar Code method 1  ");
            for (unsigned int i = 0; i < 13; i++)
            {
                printf(" %d ", finalBarCode[i]);
            }

            printf(checkSum(finalBarCode) ? "true" : "false");
        }
    }
}

string barcode_rhs(long double st, vector<vector<int>> bars, long double avbl1, long double avbl0, long double d)
{

    long double delta = d;
    long double step = st;

    long double avgBarLength0 = avbl0;
    long double avgBarLength1 = avbl1;
    vector<vector<int>> s = bars;
    vector<vector<vector<int>>> s_lhs;
    vector<vector<vector<int>>> s_rhs;

    //	 s_lhs=lhsBars(s);

    s_rhs = rhsBars(s);

    /**
          vector<long double> maxProbabilities_lhs = {0, 0, 0, 0, 0, 0};
          vector<int> maxProbabilitiesDigits_lhs = {0, 0, 0, 0, 0, 0};
          vector<int> maxProbabilitiesIndices_lhs = {0, 0, 0, 0, 0, 0};
          int maxProbIndex_lhs = 0;
          long double maxProb_lhs=0;
          vector<int> maxProbabilitiesDigits2_lhs = {0, 0, 0, 0, 0, 0};
          vector<int> maxProbabilitiesIndices2_lhs = {0, 0, 0, 0, 0, 0};
          vector<long double> maxProbabilities2_lhs = {0, 0, 0, 0, 0, 0};
          vector<long double> finalProbList_lhs = {};
          vector<double> parity_lhs = {};
    */

    vector<long double> maxProbabilities_rhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesDigits_rhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesIndices_rhs = {0, 0, 0, 0, 0, 0};
    int maxProbIndex_rhs = 0;
    long double maxProb_rhs = 0;
    vector<int> maxProbabilitiesDigits2_rhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesIndices2_rhs = {0, 0, 0, 0, 0, 0};
    vector<long double> maxProbabilities2_rhs = {
        0,
        0,
        0,
        0,
        0,
        0,
    };
    vector<long double> finalProbList_rhs = {};
    vector<double> parity_rhs = {};
    vector<int> finalbarCode = {};

    long double sigma = 1;
    int digitNumber = -1;

    /**
         for(long double c = (avgBarLength1 - delta); c <= (avgBarLength1 + delta); c = c + step)
          {
             for(long double e = (avgBarLength0 - delta); e <= (avgBarLength0 + delta); e = e + step)
            {


              vector<long double> r0L ={c*3, e*2, c, e}; ;
              vector<long double> r0G = {c, e, c*2, e*3};

              vector<long double> r1L = {c*2, e*2, c*2, e};
              vector<long double> r1G = {c, e*2, c*2, e*2};

              vector<long double> r2L = {c*2, e, c*2, e*2};
              vector<long double> r2G = {c*2, e*2, c, e*2};

              vector<long double> r3L = {c, e*3, c, e};
              vector<long double> r3G = {c, e, c*4, e};

              vector<long double> r4L = {c, e, c*3, e*2};
              vector<long double> r4G = {c*2, e*3, c, e};

              vector<long double> r5L = {c, e*2, c*3, e};
              vector<long double> r5G = {c, e*3, c*2, e};

              vector<long double> r6L = {c, e, c, e*4};
              vector<long double> r6G = {c*4, e, c, e};

              vector<long double> r7L = {c, e*3, c, e*2};
              vector<long double> r7G = {c*2, e, c*3, e};

              vector<long double> r8L = {c, e*2, c, e*3};
              vector<long double> r8G = {c*3, e, c*2, e};

              vector<long double> r9L = {c*3, e, c, e*2};
              vector<long double> r9G = {c*2, e, c, e*3};


              vector<vector<long double>> r = {r0L, r0G, r1L, r1G,  r2L, r2G,  r3L, r3G,
                r4L, r4G,  r5L, r5G,  r6L, r6G,  r7L, r7G, r8L,
                 r8G,  r9L, r9G};





              for(int a = 0; a < 6; a++)
              {

               maxProb_lhs = 0;
               maxProbIndex_lhs = 0;

               for(int b = 0; b < 20; b++)
               {

                if(p_lhs(b, a, c, e, s_lhs) >= maxProb_lhs)
                {
                    maxProb_lhs = p_lhs(b, a, c, e, s_lhs);

                    maxProbIndex_lhs = b;
                }

               }


               sigma = sigma*maxProb_lhs;





               digitNumber = ((maxProbIndex_lhs - maxProbIndex_lhs%2)/2);




               if(maxProb_lhs >= maxProbabilities_lhs[a])
               {
                maxProbabilities_lhs[a] = maxProb_lhs;
                maxProbabilitiesDigits_lhs[a] = digitNumber;
                maxProbabilitiesIndices_lhs[a] = maxProbIndex_lhs;
               }



               }



              for(unsigned  int i = 0; i < 6; i++)
              {
               if(maxProbabilities_lhs[i] >= maxProbabilities2_lhs[i])
               {
                maxProbabilities2_lhs[i] = maxProbabilities_lhs[i];
                maxProbabilitiesDigits2_lhs[i] = maxProbabilitiesDigits_lhs[i];
                maxProbabilitiesIndices2_lhs[i] = maxProbabilitiesIndices_lhs[i];
               }

              }




             }}

    */

    sigma = 1;
    for (long double c = (avgBarLength1 - delta); c <= (avgBarLength1 + delta); c = c + step)
    {
        for (long double e = (avgBarLength0 - delta); e <= (avgBarLength0 + delta); e = e + step)
        {

            vector<long double> r0R = {e * 3, c * 2, e, c};

            vector<long double> r1R = {e * 2, c * 2, e * 2, c};

            vector<long double> r2R = {e * 2, c, e * 2, c * 2};

            vector<long double> r3R = {e, c * 3, e, c};

            vector<long double> r4R = {e, c, e * 3, c * 2};

            vector<long double> r5R = {e, c * 2, e * 3, c};

            vector<long double> r6R = {e, c, e, c * 4};

            vector<long double> r7R = {e, c * 3, e, c * 2};

            vector<long double> r8R = {e, c * 2, e, c * 3};

            vector<long double> r9R = {e * 3, c, e, c * 2};

            vector<vector<long double>> r = {r0R, r1R, r2R, r3R,
                                             r4R, r5R, r6R, r7R, r8R, r9R};

            for (int a = 0; a < 6; a++)
            {

                maxProb_rhs = 0;
                maxProbIndex_rhs = 0;

                for (int b = 0; b < 10; b++)
                {

                    if (p_rhs(b, a, c, e, s_rhs) >= maxProb_rhs)
                    {
                        maxProb_rhs = p_rhs(b, a, c, e, s_rhs);

                        maxProbIndex_rhs = b;
                    }
                }

                sigma = sigma * maxProb_rhs;

                digitNumber = maxProbIndex_rhs;

                if (maxProb_rhs >= maxProbabilities_rhs[a])
                {
                    maxProbabilities_rhs[a] = maxProb_rhs;
                    maxProbabilitiesDigits_rhs[a] = digitNumber;
                    maxProbabilitiesIndices_rhs[a] = maxProbIndex_rhs;
                }
            }

            for (unsigned int i = 0; i < 6; i++)
            {
                if (maxProbabilities_rhs[i] >= maxProbabilities2_rhs[i])
                {
                    maxProbabilities2_rhs[i] = maxProbabilities_rhs[i];
                    maxProbabilitiesDigits2_rhs[i] = maxProbabilitiesDigits_rhs[i];
                    maxProbabilitiesIndices2_rhs[i] = maxProbabilitiesIndices_rhs[i];
                }
            }

            //		 }}

            vector<int> finalBarCode = {};
            /**     for(unsigned  int i = 0; i < 6; i++)
                 {
                    finalBarCode.push_back(maxProbabilitiesDigits_lhs[i]);
                 }
                 */
            for (unsigned int i = 0; i < 6; i++)
            {
                finalBarCode.push_back(maxProbabilitiesDigits_rhs[i]);
            }

            /**
                     for(unsigned int i = 0; i < 6; i++)
                     {
                        if((maxProbabilitiesIndices_lhs[i])%2 == 0)
                        {
                          parity_lhs.push_back(1);

                        }
                        if((maxProbabilitiesIndices_lhs[i])%2 == 1)
                        {
                            parity_lhs.push_back(0);

                        }


                     }


                       vector<vector<int>> m = {{1, 1, 1, 1, 1, 1}, {1, 1, 0, 1, 0, 0}, {1, 1, 0, 0, 1, 0}, {1,
                           1, 0, 0, 0, 1}, {1, 0, 1, 1, 0, 0}, {1, 0, 0, 1, 1, 0}, {1, 0,
                          0, 0, 1, 1}, {1, 0, 1, 0, 1, 0}, {1, 0, 1, 0, 0, 1}, {1, 0, 0, 1,
                           0, 1}};

                       int mNumber = -1;

                       for(int i = 0; i < 10; i++)
                       {
                         if(m.at(i).at(0) == parity_lhs.at(0) &&
                            m.at(i).at(1) == parity_lhs.at(1) &&
                             m.at(i).at(2) == parity_lhs.at(2) &&
                             m.at(i).at(3) == parity_lhs.at(3) &&
                             m.at(i).at(4) == parity_lhs.at(4) &&
                             m.at(i).at(5) == parity_lhs.at(5) )
                         {
                             mNumber = i ; break;
                         }
                       }
                       finalBarCode.insert(finalBarCode.begin() , mNumber);
            */
            string rhs_barcode = to_string(finalBarCode[0]) + to_string(finalBarCode[1]) + to_string(finalBarCode[2]) + to_string(finalBarCode[3]) + to_string(finalBarCode[4]) + to_string(finalBarCode[5]);

            // printf("\nfinal Bar Code method ");
            for (unsigned int i = 0; i < 6; i++)
            {
                //  printf(" %d ",finalBarCode[i]);
            }

            // printf(checkSum(finalBarCode) ? "true" : "false");

            return rhs_barcode;
        }
    }
}

string barcode_lhs(long double st, vector<vector<int>> bars, long double avbl1, long double avbl0, long double d)
{

    long double delta = d;
    long double step = st;

    long double avgBarLength0 = avbl0;
    long double avgBarLength1 = avbl1;
    vector<vector<int>> s = bars;
    vector<vector<vector<int>>> s_lhs;
    //	 vector<vector<vector<int>>>  s_rhs;

    s_lhs = lhsBars(s);

    // s_rhs=rhsBars(s);

    vector<long double> maxProbabilities_lhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesDigits_lhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesIndices_lhs = {0, 0, 0, 0, 0, 0};
    int maxProbIndex_lhs = 0;
    long double maxProb_lhs = 0;
    vector<int> maxProbabilitiesDigits2_lhs = {0, 0, 0, 0, 0, 0};
    vector<int> maxProbabilitiesIndices2_lhs = {0, 0, 0, 0, 0, 0};
    vector<long double> maxProbabilities2_lhs = {0, 0, 0, 0, 0, 0};
    vector<long double> finalProbList_lhs = {};
    vector<double> parity_lhs = {};

    /**
          vector<long double> maxProbabilities_rhs = {0, 0, 0, 0, 0, 0};
          vector<int> maxProbabilitiesDigits_rhs = {0, 0, 0, 0, 0, 0};
          vector<int> maxProbabilitiesIndices_rhs = {0, 0, 0, 0, 0, 0};
          int maxProbIndex_rhs = 0;
          long double maxProb_rhs=0;
          vector<int> maxProbabilitiesDigits2_rhs = {0, 0, 0, 0, 0, 0};
          vector<int> maxProbabilitiesIndices2_rhs = {0, 0, 0, 0, 0, 0};
          vector<long double> maxProbabilities2_rhs = {0, 0, 0, 0, 0, 0,};
          vector<long double> finalProbList_rhs = {};
          vector<double> parity_rhs = {};
          */
    vector<int> finalbarCode = {};

    long double sigma = 1;
    int digitNumber = -1;

    for (long double c = (avgBarLength1 - delta); c <= (avgBarLength1 + delta); c = c + step)
    {
        for (long double e = (avgBarLength0 - delta); e <= (avgBarLength0 + delta); e = e + step)
        {

            vector<long double> r0L = {c * 3, e * 2, c, e};
            ;
            vector<long double> r0G = {c, e, c * 2, e * 3};

            vector<long double> r1L = {c * 2, e * 2, c * 2, e};
            vector<long double> r1G = {c, e * 2, c * 2, e * 2};

            vector<long double> r2L = {c * 2, e, c * 2, e * 2};
            vector<long double> r2G = {c * 2, e * 2, c, e * 2};

            vector<long double> r3L = {c, e * 3, c, e};
            vector<long double> r3G = {c, e, c * 4, e};

            vector<long double> r4L = {c, e, c * 3, e * 2};
            vector<long double> r4G = {c * 2, e * 3, c, e};

            vector<long double> r5L = {c, e * 2, c * 3, e};
            vector<long double> r5G = {c, e * 3, c * 2, e};

            vector<long double> r6L = {c, e, c, e * 4};
            vector<long double> r6G = {c * 4, e, c, e};

            vector<long double> r7L = {c, e * 3, c, e * 2};
            vector<long double> r7G = {c * 2, e, c * 3, e};

            vector<long double> r8L = {c, e * 2, c, e * 3};
            vector<long double> r8G = {c * 3, e, c * 2, e};

            vector<long double> r9L = {c * 3, e, c, e * 2};
            vector<long double> r9G = {c * 2, e, c, e * 3};

            vector<vector<long double>> r = {r0L, r0G, r1L, r1G, r2L, r2G, r3L, r3G,
                                             r4L, r4G, r5L, r5G, r6L, r6G, r7L, r7G, r8L,
                                             r8G, r9L, r9G};

            for (int a = 0; a < 6; a++)
            {

                maxProb_lhs = 0;
                maxProbIndex_lhs = 0;

                for (int b = 0; b < 20; b++)
                {

                    if (p_lhs(b, a, c, e, s_lhs) >= maxProb_lhs)
                    {
                        maxProb_lhs = p_lhs(b, a, c, e, s_lhs);

                        maxProbIndex_lhs = b;
                    }
                }

                sigma = sigma * maxProb_lhs;

                digitNumber = ((maxProbIndex_lhs - maxProbIndex_lhs % 2) / 2);

                if (maxProb_lhs >= maxProbabilities_lhs[a])
                {
                    maxProbabilities_lhs[a] = maxProb_lhs;
                    maxProbabilitiesDigits_lhs[a] = digitNumber;
                    maxProbabilitiesIndices_lhs[a] = maxProbIndex_lhs;
                }
            }

            for (unsigned int i = 0; i < 6; i++)
            {
                if (maxProbabilities_lhs[i] >= maxProbabilities2_lhs[i])
                {
                    maxProbabilities2_lhs[i] = maxProbabilities_lhs[i];
                    maxProbabilitiesDigits2_lhs[i] = maxProbabilitiesDigits_lhs[i];
                    maxProbabilitiesIndices2_lhs[i] = maxProbabilitiesIndices_lhs[i];
                }
            }
        }
    }

    sigma = 1;
    /**
for(long double c = (avgBarLength1 - delta); c <= (avgBarLength1 + delta); c = c + step)
{
  for(long double e = (avgBarLength0 - delta); e <= (avgBarLength0 + delta); e = e + step)
 {



   vector<long double> r0R = {e*3, c*2, e, c};

   vector<long double> r1R = {e*2, c*2, e*2, c};

   vector<long double> r2R = {e*2, c, e*2, c*2};

   vector<long double> r3R = {e, c*3, e, c};

   vector<long double> r4R = {e, c, e*3, c*2};

   vector<long double> r5R = {e, c*2, e*3, c};

   vector<long double> r6R = {e, c, e, c*4};

   vector<long double> r7R = {e, c*3, e, c*2};

   vector<long double> r8R = {e, c*2, e, c*3};

   vector<long double> r9R = {e*3, c, e, c*2};

   vector<vector<long double>> r = { r0R, r1R,  r2R,  r3R,
      r4R,  r5R,  r6R,  r7R,  r8R,  r9R};





   for(int a = 0; a < 6; a++)
   {

    maxProb_rhs = 0;
    maxProbIndex_rhs = 0;

    for(int b = 0; b < 10; b++)
    {

     if(p_rhs(b, a, c, e, s_rhs) >= maxProb_rhs)
     {
         maxProb_rhs = p_rhs(b, a, c, e, s_rhs);

         maxProbIndex_rhs = b;
     }

    }


    sigma = sigma*maxProb_rhs;





    digitNumber = maxProbIndex_rhs;



    if(maxProb_rhs >= maxProbabilities_rhs[a])
    {
     maxProbabilities_rhs[a] = maxProb_rhs;
     maxProbabilitiesDigits_rhs[a] = digitNumber;
     maxProbabilitiesIndices_rhs[a] = maxProbIndex_rhs;
    }



    }



   for(unsigned  int i = 0; i < 6; i++)
   {
    if(maxProbabilities_rhs[i] >= maxProbabilities2_rhs[i])
    {
     maxProbabilities2_rhs[i] = maxProbabilities_rhs[i];
     maxProbabilitiesDigits2_rhs[i] = maxProbabilitiesDigits_rhs[i];
     maxProbabilitiesIndices2_rhs[i] = maxProbabilitiesIndices_rhs[i];
    }

   }





*/

    vector<int> finalBarCode = {};
    for (unsigned int i = 0; i < 6; i++)
    {
        finalBarCode.push_back(maxProbabilitiesDigits_lhs[i]);
    }
    /**	     for(unsigned  int i = 0; i < 6; i++)
             {
                  finalBarCode.push_back(maxProbabilitiesDigits_rhs[i]);
             }
    */

    for (unsigned int i = 0; i < 6; i++)
    {
        if ((maxProbabilitiesIndices_lhs[i]) % 2 == 0)
        {
            parity_lhs.push_back(1);
        }
        if ((maxProbabilitiesIndices_lhs[i]) % 2 == 1)
        {
            parity_lhs.push_back(0);
        }
    }

    vector<vector<int>> m = {{1, 1, 1, 1, 1, 1}, {1, 1, 0, 1, 0, 0}, {1, 1, 0, 0, 1, 0}, {1, 1, 0, 0, 0, 1}, {1, 0, 1, 1, 0, 0}, {1, 0, 0, 1, 1, 0}, {1, 0, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0}, {1, 0, 1, 0, 0, 1}, {1, 0, 0, 1, 0, 1}};

    int mNumber = -1;

    for (int i = 0; i < 10; i++)
    {
        if (m.at(i).at(0) == parity_lhs.at(0) &&
            m.at(i).at(1) == parity_lhs.at(1) &&
            m.at(i).at(2) == parity_lhs.at(2) &&
            m.at(i).at(3) == parity_lhs.at(3) &&
            m.at(i).at(4) == parity_lhs.at(4) &&
            m.at(i).at(5) == parity_lhs.at(5))
        {
            mNumber = i;
            break;
        }
    }
    finalBarCode.insert(finalBarCode.begin(), mNumber);

    string lhs_string = "";

    //  printf("\nfinal Bar Code LHS %s\n",lhs_string);
    for (unsigned int i = 0; i < 7; i++)
    {
        lhs_string = lhs_string + to_string(finalBarCode[i]);
        // printf(" %d ",finalBarCode[i]);
    }

    return lhs_string;
    //	    	printf(checkSum(finalBarCode) ? "true" : "false");
}
string int_to_string(int i)
{
    switch (i)
    {
    case 0:
        return "0";
    case 1:
        return "1";
    case 2:
        return "2";
    case 3:
        return "3";
    case 4:
        return "4";
    case 5:
        return "5";
    case 6:
        return "6";
    case 7:
        return "7";
    case 8:
        return "8";
    case 9:
        return "9";
    default:
        return "";
    }
}

// void Line(byte** input, float x1, float y1, float x2, float y2, const int color );

// 1. 输入图片被缩放为 640 x 480 分辨率
// 2. 图像处理后找到条码的所在框（OpenCV处理）
// 3. 输出条码的第一条线坐标和最后一条线坐标
// 4. 按照条码的标准解析条码线
// 5. 解析错误，换个位置重做第4步
// 6. 解析成功则返回条码信息
int main(int argc, char *argv[])
{
    bool debug = false;

    string file_name = argv[1];

    if (strcmp(argv[2], "true") == 0)
    {
        debug = true;
    }

    //    string file_name="C:\\codeblocks\\barcodeLocalisationRead_linux\\barcodes_mbdb\\barcode_2.jpg";
    string cropped_file_name;

    string new_dir_name;
    size_t last_token_index = file_name.find_last_of("/\\");
    new_dir_name = file_name.substr(0, last_token_index) + "\\";

    //   printf("\n new_dir_name %s",new_dir_name.c_str());
    //   printf("\n file_name %s",file_name.c_str());
    cv::Mat img;
    img = imread(file_name, cv::IMREAD_COLOR);
    cv::resize(img, img, cv::Size(640, 480));

    int height = img.rows;
    int width = img.cols;
    byte **matrix = new byte *[height];

    for (int i = 0; i < height; i++)
    {
        matrix[i] = new byte[width];
    }

    for (int i = 0; i < height; i++)
    {

        for (int j = 0; j < width; j++)
        {

            // BGR not RGB.
            long double blue = img.at<cv::Vec3b>(i, j)[0];
            long double green = img.at<cv::Vec3b>(i, j)[1];
            long double red = img.at<cv::Vec3b>(i, j)[2];
            ;
            // Y (x) = 0.299R(x) + 0.587G(x) + 0.114B(x)
            long double intensity = (red * 0.299 + green * 0.587 + blue * 0.114);

            matrix[i][j] = intensity;
        }
    }

    int he = 0;
    int wi = 0;
    byte **ret = new byte *[he];
    for (int i = 0; i < he; i++)
    {
        ret[i] = new byte[wi];
    }
    clock_t start = clock();
    // The following loop is necessary as not all barcodes are rectangular ona flat surface.
    // Some barcodes are wider at the top and narrower at the bottom as the camera and barcode are at an angle.
    // Some barcodes are also on curved surfaces whcih distorts the image.

    string lhs = "0000000";
    string rhs = "000000";
    // 图像处理一次
    // crop存在偏差可能性，即可能不是恰好的条码框
    cv::Rect crop = image_process(matrix, height, width, 150, he, wi, new_dir_name, ret, debug);
    // 条码只需要一次行扫描
    // 如果扫描失败，换个位置尝试。这儿循环相当于取了三行尝试
    for (double k = 0.25; k < 1.0; k = k + 0.25)
    {
        int left = 0;
        vector<long double> scanline = getScanline(ret, he, wi, debug, k, left);
        crop.x = crop.x + left;
        // cv::Point pt1(min_x, min_y), pt2(max_x, max_y);
        cv::rectangle(img, crop.tl(), crop.br(), cv::Scalar(0,255,0), 3);
        cv::imwrite(new_dir_name + "orig.jpg", img);

        vector<vector<int>> Bars = bars(scanline, debug);

        vector<int> lPadding = leftPadding(Bars);
        vector<int> rPadding = rightPadding(Bars);
        vector<vector<int>> lGuard = leftGuard(Bars);
        vector<vector<int>> rGuard = rightGuard(Bars);
        vector<vector<int>> mGuard = midGuard(Bars);
        vector<vector<vector<int>>> lAndRBars = lhsAndRhsBars(Bars);
        vector<vector<vector<int>>> lBars = lhsBars(Bars);
        vector<vector<vector<int>>> rBars = rhsBars(Bars);
        double avbl1 = averageBarLength1(Bars);
        double avbl0 = averageBarLength0(Bars);

        if (debug)
        {

            printf("\naverage bar length 0 %f\naverage bar length 1 %f", avbl0, avbl1);

            printf("\nleft padding %llu right padding %llu", lPadding.size(), rPadding.size());

            printf("\nleft guard");
            for (unsigned int i = 0; i < lGuard.size(); i++)
            {
                printf("\n");
                for (unsigned int j = 0; j < lGuard.at(i).size(); j++)
                {
                    printf("%d", lGuard.at(i).at(j));
                }
            }
            printf("\nright guard");
            for (unsigned int i = 0; i < rGuard.size(); i++)
            {
                printf("\n");
                for (unsigned int j = 0; j < rGuard.at(i).size(); j++)
                {
                    printf("%d", rGuard.at(i).at(j));
                }
            }
            printf("\nmid guard");
            for (unsigned int i = 0; i < mGuard.size(); i++)
            {
                printf("\n");
                for (unsigned int j = 0; j < mGuard.at(i).size(); j++)
                {
                    printf("%d", mGuard.at(i).at(j));
                }
            }
        }

        std::future<string> lhsDigits = std::async(&barcode_lhs, 0.1, Bars, avbl1, avbl0, 0);
        lhs = lhsDigits.get();

        std::future<string> rhsDigits = std::async(&barcode_rhs, 0.1, Bars, avbl1, avbl0, 0);
        rhs = rhsDigits.get();

        vector<int> barcode;
        for (int i = 0; i < 7; i++)
        {

            string temp;
            if (temp.compare("-1"))
                temp = "0";

            temp.push_back(lhs.at(i));
            barcode.push_back(stoi(temp));
        }
        for (int i = 0; i < 6; i++)
        {
            string temp;
            temp.push_back(rhs.at(i));
            barcode.push_back(stoi(temp));
        }

        string csum = "blank";
        if (checkSum(barcode))
        {
            csum = "true";
            printf("\nchecksum %s", csum.c_str());
            break;
        }
        else
        {
            // csum="false";
            // printf("\nchecksum %s",csum.c_str());
        }
    }
    cout << "\nlhs " << lhs << "\nrhs " << rhs;
    clock_t end = clock();
    double time = (double)(end - start);
    printf("\ntime %f\n", time / 1000);

    return 0;
}
