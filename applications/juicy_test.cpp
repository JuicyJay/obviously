/*
 * juicy_test.cpp
 *
 *      Author: jasmin
 */

#include "math.h"
#include "obcore/base/System.h"
#include "obcore/math/linalg/eigen/Matrix.h"
#include "obcore/math/mathbase.h"
#include "obgraphic/Obvious3D.h"
#include "obvision/reconstruct/Sensor.h"
#include <iostream>

double azimuthResolution = 3.0;
double inclResolution    = 2.0;
int**  _indexMap;

/**
 * returns azimuth angle and inclination angle
 * @param[in] xCoord x coordinate of a 3D point
 * @param[in] yCoord y coordinate of a 3D point
 * @param[in] zCoord z coordinate of a 3D point
 * @param[out] inclAngle inclination angle in x-z-plane, 16 layers of rays from -15° to +15°, 2° resolution (VLP16 PUCK)
 * @param[out] azimAngle azimuth angle in x-y-plane, 0° to 360°
 */
void returnAngles(double xCoord, double yCoord, double zCoord, double* inclAngle, double* azimAngle)
{
  // INCLINATION
  double theta = 0.0; // careful - this is the angle between z-axis and x-y-plane as defined in polar coords --> not inclination angle

  // das macht May mit acos
  //  double r = sqrt(coords3D(0,i) * coords3D(0,i) + coords3D(1,i) * coords3D(1,i) + coords3D(2,i) * coords3D(2,i));
  //  double theta = acos(coords3D(1,i) / r);
  //  if(coords3D(2,i)>0)
  //    theta = -theta;catkin

  double length     = sqrt(xCoord * xCoord + yCoord * yCoord + zCoord * zCoord);
  double lengthInv  = 1.0 / length;
  double valForAcos = zCoord * lengthInv;
  theta             = obvious::rad2deg(acos(valForAcos));
  std::cout << __PRETTY_FUNCTION__ << " theta = " << theta << std::endl;

  if(theta > 90.0)
  {
    *inclAngle = -(theta - 90.0);
  }
  else
  {
    *inclAngle = 90.0 - theta;
  }
  std::cout << __PRETTY_FUNCTION__ << " inclAngle = " << *inclAngle << std::endl;

  // AZIMUTH
  *azimAngle = obvious::rad2deg(atan2(yCoord, xCoord));
  if(*azimAngle < 0)
  {
    *azimAngle += 360.0; // express angles positively in 3rd and 4th quadrant
  }

  std::cout << __PRETTY_FUNCTION__ << " azimAngle = " << *azimAngle << std::endl;
}

/**
 * returns ray index
 * @param[in] azimAngle azimuth angle calculated in returnAngles()
 * @param[in] inclAngle inclination angle calculated in returnAngles()
 * @param[out] azimIndex index of azimuth ray number which is closest to 3D point to assign laser data to 3D point
 * @param[out] inclIndex index of inclination ray number which is closest to 3D point to assign laser data to 3D point
 * @todo make this function abstract so each velodyne sensor has to implement it since it may differ from sensor to sensor
 */
void returnRayIndex(double azimAngle, double inclAngle, unsigned int* azimIndex, unsigned int* inclIndex)
{
  // AZIMUTH
  *azimIndex = round(azimAngle / azimuthResolution);
  std::cout << __PRETTY_FUNCTION__ << " azimIndex = " << *azimIndex << std::endl;

  // INCLINATION --> may differ for different sensors
  // assignment will always be the same for VLP16 - inclination resolution is fixed 2° and nbr of rays is 16
  double mapInclination = inclAngle + 15.0;                       // map inclination angles (-15° -> +15°) up to positive range 0° - 30°
  *inclIndex            = round(mapInclination / inclResolution); // maybe exchange round with floor ? or ceiling? idk

  std::cout << __PRETTY_FUNCTION__ << " inclIndex = " << *inclIndex << std::endl;
}

// nicht richtig?
void setIndexMap(unsigned int azimuthRays, unsigned int verticalRays)
{
  unsigned int column = 0;
  // evtl mit _height u _width ersetzen
  obvious::System<int>::allocate(azimuthRays, verticalRays, _indexMap);
  for(unsigned int row = 0; row < azimuthRays; row++) // r row
  {
    for(column = 0; column < verticalRays; column++) // c column
    {

      _indexMap[row][column] = row * (verticalRays) + column;

      //      std::cout << "_indexMap = " << _indexMap[row][column] << std::endl;
    }
    column = 0; // iterate over 16 vertical rays for each azimuth ray
  }
}

void setDistanceMap(vector<float> phiAzim, vector<float> thetaIncl, vector<float> dist) {}

// indexSensormodel: indices aus dem Sensormodell 0 = -15, 1 = -13 ...
// indexVelodyneROS: indices aus Abfeuerungsreihenfolge der Laserstrahlen 0 = -15, 1 = 1, 2 = -13, ...
unsigned int lookupIndex(int indexSensormodel)
{
  unsigned int indexVelodyneROS = 0;
  switch(indexSensormodel)
  {
  case 0:
    indexVelodyneROS = 0;
    break;
  case 1:
    indexVelodyneROS = 2;
    break;
  case 2:
    indexVelodyneROS = 4;
    break;
  case 3:
    indexVelodyneROS = 6;
    break;
  case 4:
    indexVelodyneROS = 8;
    break;
  case 5:
    indexVelodyneROS = 10;
    break;
  case 6:
    indexVelodyneROS = 12;
    break;
  case 7:
    indexVelodyneROS = 14;
    break;
  case 8:
    indexVelodyneROS = 1;
    break;
  case 9:
    indexVelodyneROS = 3;
    break;
  case 10:
    indexVelodyneROS = 5;
    break;
  case 11:
    indexVelodyneROS = 7;
    break;
  case 12:
    indexVelodyneROS = 9;
    break;
  case 13:
    indexVelodyneROS = 11;
    break;
  case 14:
    indexVelodyneROS = 13;
    break;
  case 15:
    indexVelodyneROS = 15;
    break;
  }

  return indexVelodyneROS;
}

// test with selfmade matrix here, later with incoming matrix* M & transformation into sensor coordinate system
void backProject(obvious::Matrix* M, int* indices, obvious::Matrix* T = 0)
{
  // std::cout << __PRETTY_FUNCTION__ << " getRows = " << M->getRows() << std::endl;
  // std::cout << __PRETTY_FUNCTION__ << " getCols = " << M->getCols() << std::endl;

  obvious::Matrix copyM = *M;

  // LATER GETROWS WIE BEI MAY WEGEN TRANSPOSE??
  // for(unsigned int i = 0; i < M->getCols(); i++)
  // {
  //   double x = copyM(0, i);
  //   double y = copyM(1, i);
  //   double z = copyM(2, i);
  //   std::cout << __PRETTY_FUNCTION__ << " x = " << x << std::endl;
  //   std::cout << __PRETTY_FUNCTION__ << " y = " << y << std::endl;
  //   std::cout << __PRETTY_FUNCTION__ << " z = " << z << std::endl;
  //   double inclAngle = 0.0;
  //   double azimAngle = 0.0;
  //   returnAngles(x, y, z, &inclAngle, &azimAngle);
  //   std::cout << __PRETTY_FUNCTION__ << " inclAngle = " << inclAngle << std::endl;
  //   std::cout << __PRETTY_FUNCTION__ << " azimAngle = " << azimAngle << std::endl;

  //   // 1: calculate incoming azimuth = ROW index of indexMap
  //   // 2: calculate incoming inclination = COLUMN of indexMap
  //   unsigned int row    = 0;
  //   unsigned int column = 0;
  //   returnRayIndex(azimAngle, inclAngle, &row, &column);
  //   std::cout << __PRETTY_FUNCTION__ << " row = " << row << std::endl;
  //   std::cout << __PRETTY_FUNCTION__ << " column = " << column << std::endl;
  //   // map column from sensor model to Velodyne firing sequence (order of vertical rays differs between sensormodel and velodyne ros input)
  //   unsigned int columnMapped = lookupIndex(column);
  //   std::cout << __PRETTY_FUNCTION__ << " columnMapped = " << columnMapped << std::endl;

  //   // push value into int* that is in indexMap[row][column]
  //   std::cout << __PRETTY_FUNCTION__ << "_indexMap = " << _indexMap[row][columnMapped] << std::endl;

  //   // probe: index ausrechnen
  //   unsigned int idxCheck = columnMapped + 16 * row;
  //   std::cout << __PRETTY_FUNCTION__ << " idxCheck = " << idxCheck << std::endl;

  //   indices[i] = _indexMap[row][columnMapped];
  //   std::cout << "end of loop, content of vector indices = " << indices[i] << std::endl;
  // }

  // transform M into sensor coordinate system
  //  obvious::Matrix PoseInv = obvious::getTransformation();
  //  PoseInv.invert();
  //  if(T)
  //    PoseInv *= *T;
  //
  //  //multiply PoseInv with M where poseInv is not transposed but M is transposed (true)
  //  obvious::Matrix coords3D = obvious::Matrix::multiply(PoseInv, *M, false, true);
  //  obvious::Matrix coords3D;
  //
  //  //loop through all 3d points in matrix and calculate inclination angle phi, azimuth angle theta and length vector (length from 3D pt to 0)
  //  //WARUM DURCH ROWS? UND NICHT COLUMNS? WEIL TRANSPOSED WURDE? ne du idiot, i steht ja an der column stelle. hääää
  //  for(unsigned int i=0; i<M->getRows(); i++)
  //    {
  //      double phi = atan2(coords3D(2,i), coords3D(0,i)) - M_PI;
  //      if(phi>M_PI) phi -= M_PI;
  //      if(phi<-M_PI) phi += M_PI;
  //
  //      double r = sqrt(coords3D(0,i) * coords3D(0,i) + coords3D(1,i) * coords3D(1,i) + coords3D(2,i) * coords3D(2,i));
  //      double theta = acos(coords3D(1,i) / r);
  //      if(coords3D(2,i)>0)
  //        theta = -theta;
  //
  //      //nächste codeschritte aus may's code: rechnen für den swiipenden scanner, bei mir 360°
  //    }
}

int main(void)
{
  double              bg[3]  = {1.0, 1.0, 1.0};
  obvious::Obvious3D* viewer = new obvious::Obvious3D("VLP16 Raycast", 1024, 768, 0, 0, bg);
  double**            coordsStart;
  double**            coordsEnd;

  //  double** coordsStartColor;
  //  double** coordsEndColor;

  double rgbGreen[3] = {0.0, 180.0, 153.0};
  //  double rgbRed[3] = {255.0,0.0,51.0};

  float        angleInclThetaMin = -15.0; // smallest inclination angle - "lowest" ray
  double       thetaSphere       = 0.0;   // für Definition des Winkels für Kugelkoordinaten
  double       angleAzimuth      = 0.0;
  unsigned int azimuthRays       = static_cast<unsigned>(360 / azimuthResolution);
  unsigned int verticalRays      = 16;

  unsigned int totalRays = (azimuthRays * verticalRays) - 1;    // Anzahl der totalRays -1 weil Zählung bei 0 beginnt
  obvious::System<double>::allocate(totalRays, 3, coordsStart); //(row - for each point, column - xyz, nameArray)
  obvious::System<double>::allocate(totalRays, 3, coordsEnd);

  viewer->showAxes();

  // TEST BACKPROJECTION
  // draw sphere
  double centerSphere[3] = {0.2, -0.2, -0.007};
  double radiusSphere    = 0.05;
  double z               = centerSphere[2];
  double y               = centerSphere[1];
  double x               = centerSphere[0];
  viewer->addSphere(centerSphere, radiusSphere, rgbGreen);

  double inclAngle = 0.0;
  double azimAngle = 0.0;

  returnAngles(x, y, z, &inclAngle, &azimAngle);
  std::cout << __PRETTY_FUNCTION__ << " inclAngle = " << inclAngle << std::endl;
  std::cout << __PRETTY_FUNCTION__ << " azimAngle = " << azimAngle << std::endl;

  unsigned int azimIndex = 0;
  unsigned int inclIndex = 0;
  // returnRayIndex(azimAngle, inclAngle, &azimIndex, &inclIndex);

  // setIndexMap(azimuthRays, verticalRays);

  // test just three 3D points in matrix with method backProject
  obvious::Matrix testMatrix(3, 4);

  // p1
  testMatrix(0, 0) = 0.2;
  testMatrix(1, 0) = 0.01;
  testMatrix(2, 0) = -0.0077;
  // p2
  testMatrix(0, 1) = 0.1;
  testMatrix(1, 1) = 0.01;
  testMatrix(2, 1) = 0;
  // p3
  testMatrix(0, 2) = 0.15;
  testMatrix(1, 2) = 0.01;
  testMatrix(2, 2) = 0.0077;
  // p4
  testMatrix(0, 3) = 0.175;
  testMatrix(1, 3) = 0.01;
  testMatrix(2, 3) = 0.0077;
  obvious::Matrix* test = &testMatrix;
  /////////////////////////////////////////// DAS HIER SPÄTER RAUS UND IN AUFRUFENDER METHODE RICHTIG DIMENSIONIEREN!!!! getrows oder so
  int* indices = new int[4];
  //////////////////////////////
  // backProject(test, indices);

  obvious::Matrix* _rays;
  _rays = new obvious::Matrix(3, totalRays);

  for(unsigned int i = 0; i < azimuthRays; i++)
  {
    for(unsigned int j = 0; j < 16; j++)
    {
      if(angleInclThetaMin < 0)
      {
        thetaSphere = 90.0 + (angleInclThetaMin * (-1));
      }
      else
      {
        thetaSphere = 90.0 - angleInclThetaMin;
      }
      coordsStart[j][0] = 0.0;
      coordsStart[j][1] = 0.0;
      coordsStart[j][2] = 0.0;

      // längerer Ray, der durch Kugel geht
      if((i == azimIndex) && (j == inclIndex))

      {
        std::cout << "angleAzimuth of current ray is " << angleAzimuth << std::endl;
        std::cout << "angleInclination of current ray is " << angleInclThetaMin << std::endl;
        coordsEnd[j][0] = 0.5 * sin((thetaSphere * M_PI) / 180) * cos((angleAzimuth * M_PI) / 180);
        coordsEnd[j][1] = 0.5 * sin((thetaSphere * M_PI) / 180) * sin((angleAzimuth * M_PI) / 180);
        coordsEnd[j][2] = 0.5 * cos((thetaSphere * M_PI) / 180);
      }
      else
      {
        coordsEnd[j][0] = 0.1 * sin((thetaSphere * M_PI) / 180) * cos((angleAzimuth * M_PI) / 180);
        coordsEnd[j][1] = 0.1 * sin((thetaSphere * M_PI) / 180) * sin((angleAzimuth * M_PI) / 180);
        coordsEnd[j][2] = 0.1 * cos((thetaSphere * M_PI) / 180);
      }

      // for length calc
      obvious::Matrix ray(3, 1);
      ray(0, 0) = coordsEnd[j][0];
      ray(1, 0) = coordsEnd[j][1];
      ray(2, 0) = coordsEnd[j][2];
      // normalize length
      const double length    = sqrt(ray(0, 0) * ray(0, 0) + ray(1, 0) * ray(1, 0) + ray(2, 0) * ray(2, 0));
      const double lengthInv = 1.0 / length;
      // store in pointer _rays which is later inherited from class Sensor
      (*_rays)(0, j) = ray(0, 0) * lengthInv;
      (*_rays)(1, j) = ray(1, 0) * lengthInv;
      (*_rays)(2, j) = ray(2, 0) * lengthInv;
      // todo: indexMap

      angleInclThetaMin += inclResolution;
    }
    angleInclThetaMin = -15.0; // austauschen -15.0 gegen inclMin
    angleAzimuth += azimuthResolution;

    viewer->addLines(coordsStart, coordsEnd, totalRays, NULL);
  }

  viewer->startRendering();
  delete viewer;
} // main
