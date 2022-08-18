/*
RUN
./build/pose_estimation_3d2d_rs435_solvePnP 1.png 2.png 1_depth.png

Test move camera with aruco as pose ground truth  (scene fixe)
  /home/lc/datasets/pxEnvVO/datasetRealsenseD435iTool76/datamarkerDist50cm_Mov30cmRot30deg
  /6.png            origin
  /109.png  aprox traslation 30cm rot 30deg
RUN
  ./build/pose_estimation_3d2d_rs435_solvePnP /home/lc/datasets/pxEnvVO/datasetRealsenseD435iTool76/datamarkerDist50cm_Mov30cmRot30deg/MonoImg/6.png /home/lc/datasets/pxEnvVO/datasetRealsenseD435iTool76/datamarkerDist50cm_Mov30cmRot30deg/MonoImg/109.png /home/lc/datasets/pxEnvVO/datasetRealsenseD435iTool76/datamarkerDist50cm_Mov30cmRot30deg/DepthImgAligned/6.pgm
*/
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <chrono>

using namespace std;
using namespace cv;



class VisualOdometry
{
private:
  /* data */
public:
  VisualOdometry(const cv::Mat& _K):K(_K){};
  ~VisualOdometry(){};



void PoseEstimation3d2d(const cv::Mat& img_1,const cv::Mat& img_2, const cv::Mat& d1, 
  cv::Mat& R, cv::Mat& t )
{
  vector<KeyPoint> keypoints_1, keypoints_2;
  vector<DMatch> matches;
  find_feature_matches(img_1, img_2, keypoints_1, keypoints_2, matches);
  cout << "一共找到了" << matches.size() << "组匹配点" << endl;
   
  vector<Point3f> pts_3d;
  vector<Point2f> pts_2d;
  for (DMatch m:matches) {
    ushort d = d1.ptr<unsigned short>(int(keypoints_1[m.queryIdx].pt.y))[int(keypoints_1[m.queryIdx].pt.x)];
    if (d == 0)   // bad depth
      continue;
    //float dd = d / 5000.0; //dataset original
    float dd = d / 1000.0; //dataset rs d435i  in mm
    Point2d p1 = pixel2cam(keypoints_1[m.queryIdx].pt, K);
    pts_3d.push_back(Point3f(p1.x * dd, p1.y * dd, dd));
    pts_2d.push_back(keypoints_2[m.trainIdx].pt);
  }

  cout << "3d-2d pairs: " << pts_3d.size() << endl;

  chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
  Mat r;
  solvePnP(pts_3d, pts_2d, K, Mat(), r, t, false); // 调用OpenCV 的 PnP 求解，可选择EPNP，DLS等方法
  
  cv::Rodrigues(r, R); // r为旋转向量形式，用Rodrigues公式转换为矩阵
  chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
  chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
  cout << "solve pnp in opencv cost time: " << time_used.count() << " seconds." << endl;

}

void find_feature_matches(const Mat &img_1, const Mat &img_2,
                          std::vector<KeyPoint> &keypoints_1,
                          std::vector<KeyPoint> &keypoints_2,
                          std::vector<DMatch> &matches) {
  //-- 初始化
  Mat descriptors_1, descriptors_2;
  // used in OpenCV3
  Ptr<FeatureDetector> detector = ORB::create();
  Ptr<DescriptorExtractor> descriptor = ORB::create();
  // use this if you are in OpenCV2
  // Ptr<FeatureDetector> detector = FeatureDetector::create ( "ORB" );
  // Ptr<DescriptorExtractor> descriptor = DescriptorExtractor::create ( "ORB" );
  Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("BruteForce-Hamming");
  //-- 第一步:检测 Oriented FAST 角点位置
  detector->detect(img_1, keypoints_1);
  detector->detect(img_2, keypoints_2);

  //-- 第二步:根据角点位置计算 BRIEF 描述子
  descriptor->compute(img_1, keypoints_1, descriptors_1);
  descriptor->compute(img_2, keypoints_2, descriptors_2);

  //-- 第三步:对两幅图像中的BRIEF描述子进行匹配，使用 Hamming 距离
  vector<DMatch> match;
  // BFMatcher matcher ( NORM_HAMMING );
  matcher->match(descriptors_1, descriptors_2, match);

  //-- 第四步:匹配点对筛选
  double min_dist = 10000, max_dist = 0;

  //找出所有匹配之间的最小距离和最大距离, 即是最相似的和最不相似的两组点之间的距离
  for (int i = 0; i < descriptors_1.rows; i++) {
    double dist = match[i].distance;
    if (dist < min_dist) min_dist = dist;
    if (dist > max_dist) max_dist = dist;
  }

  printf("-- Max dist : %f \n", max_dist);
  printf("-- Min dist : %f \n", min_dist);

  //当描述子之间的距离大于两倍的最小距离时,即认为匹配有误.但有时候最小距离会非常小,设置一个经验值30作为下限.
  for (int i = 0; i < descriptors_1.rows; i++) {
    if (match[i].distance <= max(2 * min_dist, 30.0)) {
      matches.push_back(match[i]);
    }
  }
}

Point2d pixel2cam(const Point2d &p, const Mat &K) {
  return Point2d
    (
      (p.x - K.at<double>(0, 2)) / K.at<double>(0, 0),
      (p.y - K.at<double>(1, 2)) / K.at<double>(1, 1)
    );
}

  cv::Mat K;

};






int main(int argc, char **argv) {
  if (argc != 4) {
    cout << "usage: pose_estimation_3d2d img1 img2 depth1" << endl;
    return 1;
  }

//Mat K = (Mat_<double>(3, 3) << 1397.61133, 0, 976.10999, 0, 1395.06567, 532.28210, 0, 0, 1);  // calib Realsense D435i resolution 1920x1080
  Mat calibK = (Mat_<double>(3, 3) << 520.9, 0, 325.1, 0, 521.0, 249.7, 0, 0, 1); //original example
  VisualOdometry vo = VisualOdometry(calibK);

  Mat img_1 = imread(argv[1], cv::IMREAD_COLOR);
  Mat img_2 = imread(argv[2], cv::IMREAD_COLOR);
  assert(img_1.data && img_2.data && "Can not load images!");

  Mat imgd_1 = imread(argv[3], cv::IMREAD_UNCHANGED);       // 深度图为16位无符号数，单通道图像
  
  Mat R, t;
  vo.PoseEstimation3d2d(img_1,img_2,imgd_1,R,t);

  cout << "R=" << endl << R << endl;
  cout << "t=" << endl << t << endl;

  return 0;
}

