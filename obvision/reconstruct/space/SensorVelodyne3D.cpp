#include "SensorVelodyne3D.h"
#include "obcore/base/System.h"
#include "obcore/math/mathbase.h"

namespace obvious {

SensorVelodyne3D::SensorVelodyne3D(unsigned int raysIncl, double inclMin,
                                   double inclRes, double azimRes,
                                   double maxRange, double minRange,
                                   double lowReflectivityRange)
    : Sensor(3, maxRange, minRange, lowReflectivityRange) {
  _azimRes = azimRes;
  _inclRes = inclRes;
  unsigned int raysAzim = 0;
  double azimAngle = 0.0;
  double inclAngle = 0.0;
  const double resetInclMin = inclMin; // to reset variable inclMin each time
                                       // after exiting inner for-loop (=-15°
                                       // here for VLP16)

  std::cout << __PRETTY_FUNCTION__ << "_azimRes = " << _azimRes << std::endl;

  raysAzim = round(static_cast<unsigned>(360 / azimRes));
  std::cout << __PRETTY_FUNCTION__ << "raysAzim = " << raysAzim << std::endl;

  // inherited from class sensor
  // Size of first dimension, i.e., # of samples of a 2D scan or width of image
  // sensor
  _width = raysAzim;
  // Size of second dimension, i.e., 1 for a 2D scan or height of image sensor
  _height = raysIncl;
  // Number of measurement samples, i.e. _width x _height
  _size = raysAzim * raysIncl;
  std::cout << __PRETTY_FUNCTION__ << "_size = " << _size << std::endl;

  _data = new double[_size];
  _mask = new bool[_size];

  _rays = new obvious::Matrix(3, _size);
  std::cout << __PRETTY_FUNCTION__ << "raysAzim = " << raysAzim
            << " raysIncl = " << raysIncl << std::endl;

  // evtl mit _height u _width ersetzen
  obvious::System<int>::allocate(raysAzim, raysIncl, _indexMap);
  // for backproject --> SPÄTER WIRD DIE VARIABLE HIER NICHT MEHR GEBRAUCHT; DA
  // INDICES IN AUFRUFENDER METHODE (velodyne3d_ros_wrapper) DIMENSIONIERT WIRD
  // _indices = new int[4];

  // TEST
  double x = 0.2;
  double y = 0.01;
  double z = -0.007;
  returnAngles(x, y, z, &inclAngle, &azimAngle);
  std::cout << __PRETTY_FUNCTION__ << " inclAngle = " << inclAngle << std::endl;
  std::cout << __PRETTY_FUNCTION__ << " azimAngle = " << azimAngle << std::endl;

  unsigned int azimIndex = 0;
  unsigned int inclIndex = 0;
  returnRayIndex(azimAngle, inclAngle, &azimIndex, &inclIndex);

  setIndexMap(raysAzim, raysIncl);
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
  obvious::Matrix *test = &testMatrix;

  //////////////////////////////////////////////////////////////////////////////////////////
  // die methode wird nachher ja in PUSH aufgerufen, hie rnur zum test
  // backProject(test, _indices);

  unsigned int n = 0; // counts all rays of inner loop and outer loop to store
                      // in matrix _rays
  for (unsigned int i = 0; i < raysAzim; i++) {
    for (unsigned int j = 0; j < raysIncl; j++, n++) {
      // to follow definition of theta (inclination angle) in spherical
      // coordinate system
      double thetaSphere = 0.0;
      const double piHalf = deg2rad(90.0);

      if (inclMin < 0) {
        thetaSphere = piHalf + (inclMin) * (-1);
      } else {
        thetaSphere = piHalf - inclMin;
      }

      obvious::Matrix calcRays(3, _size);

      // draw long ray for the one closest to the sphere center
      float r = 0.0;
      if ((i == azimIndex) && (j == inclIndex)) {
        r = 0.5;
      } else {
        r = 0.1;
      }

      calcRays(0, 0) = r * sin(thetaSphere) * cos(azimAngle);
      calcRays(1, 0) = r * sin(thetaSphere) * sin(azimAngle);
      calcRays(2, 0) = r * cos(thetaSphere);
      // normalize rays
      const double length = sqrt(calcRays(0, 0) * calcRays(0, 0) +
                                 calcRays(1, 0) * calcRays(1, 0) +
                                 calcRays(2, 0) * calcRays(2, 0));
      const double lengthInv = 1.0 / length;

      // store end points of current ray in matrix* _rays inherited from class
      // Sensor
      (*_rays)(0, n) = calcRays(0, 0) * lengthInv;
      (*_rays)(1, n) = calcRays(1, 0) * lengthInv;
      (*_rays)(2, n) = calcRays(2, 0) * lengthInv;

      inclMin += inclRes;
    }
    inclMin = resetInclMin; // reset inclination angle to minimum
    azimAngle += azimRes;
  }
}

SensorVelodyne3D::~SensorVelodyne3D() { delete _rays; }

// return by reference - inclAngle u azimuth in DEGREEs
void SensorVelodyne3D::returnAngles(double xCoord, double yCoord, double zCoord,
                                    double *inclAngle, double *azimAngle) {
  // Inclination
  double theta = 0.0; // careful, this is angle between z-axis and x-y-plane as
                      // defined in polar coordinates --> no distinction of
                      // cases for acos necessary, bec. only values from 75° -
                      // 105° for VLP16
  double length = sqrt(xCoord * xCoord + yCoord * yCoord + zCoord * zCoord);
  double lengthInv = 1.0 / length;
  theta = obvious::rad2deg(acos(zCoord * lengthInv));

  if ((theta < 75.0) && (theta > 105.0)) {
    std::cout
        << __PRETTY_FUNCTION__
        << "3D coordinates are not within field of vision of laser scanner"
        << std::endl;
  } else {
    if (theta > 90.0) // translate theta into inclination angle "aperture angle"
                      // from -15° to +15°
    {
      *inclAngle = -(theta - 90.0); // -15° -> 0°
    } else {
      *inclAngle = 90.0 - theta; // 0° -> +15°
    }

    // Azimuth
    *azimAngle = obvious::rad2deg(atan2(yCoord, xCoord));
    if (*azimAngle < 0) {
      *azimAngle += 360.0; // express angles positively in 3rd and 4th
    }
  }
}

void SensorVelodyne3D::returnRayIndex(double azimAngle, double inclAngle,
                                      unsigned int *azimIndex,
                                      unsigned int *inclIndex) {
  *azimIndex = round(azimAngle / _azimRes);
  std::cout << __PRETTY_FUNCTION__ << " azimIndex = " << *azimIndex
            << std::endl;

  // assignment will always be the same for VLP16 - inclination resolution is
  // fixed 2° and nbr of rays is 16
  double mapInclination =
      inclAngle + 15.0; // map inclination angles (-15° -> +15°) up to positive
  // range 0° - 30° --> change so this also works for E32
  *inclIndex = round(mapInclination / _inclRes);
}

void SensorVelodyne3D::setIndexMap(unsigned int raysAzim,
                                   unsigned int raysIncl) {
  unsigned int column = 0;
  for (unsigned int row = 0; row < raysAzim; row++) {
    for (column = 0; column < raysIncl; column++) {

      _indexMap[row][column] = row * (raysIncl) + column;

      //      std::cout << "_indexMap = " << _indexMap[row][column] <<
      //      std::endl;
    }
    column = 0; // iterate over 16 vertical rays for each azimuth ray
  }
}

unsigned int SensorVelodyne3D::lookupIndex(int indexSensormodel) {
  unsigned int indexVelodyneROS = 0;
  switch (indexSensormodel) {
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

// first test with selfmade matrix before ROS pointcloud matrix comes in later
void SensorVelodyne3D::backProject(obvious::Matrix *M, int *indices, obvious::Matrix *T) 
{
  Matrix PoseInv = getTransformation();
  PoseInv.invert();
  if (T)
    PoseInv *= *T;

  // multiply  PoseInv with M where poseInv is not transposed but M is
  // transposed (true)
  // obvious::Matrix copyM = *M; das war erster test mit fake daten aus
  // constructor
  obvious::Matrix coords3D =
      obvious::Matrix::multiply(PoseInv, *M, false, true);

  std::cout << __PRETTY_FUNCTION__ << " getRows = " << M->getRows()
            << std::endl;
  std::cout << __PRETTY_FUNCTION__ << " getCols = " << M->getCols()
            << std::endl;

  // LATER GETROWS WIE BEI MAY WEGEN TRANSPOSE??
  // for(unsigned int i = 0; i < M->getCols(); i++)    //todo WARUM HIER
  // GETROWS??? muss das nicht getcolumns sein?
  for (unsigned int i = 0; i < M->getRows(); i++) 
  {
    double x = coords3D(0, i);
    double y = coords3D(1, i);
    double z = coords3D(2, i);
    std::cout << __PRETTY_FUNCTION__ << " x = " << x << std::endl;
    std::cout << __PRETTY_FUNCTION__ << " y = " << y << std::endl;
    std::cout << __PRETTY_FUNCTION__ << " z = " << z << std::endl;
    double inclAngle = 0.0;
    double azimAngle = 0.0;
    returnAngles(x, y, z, &inclAngle, &azimAngle);
    std::cout << __PRETTY_FUNCTION__ << " inclAngle = " << inclAngle
              << std::endl;
    std::cout << __PRETTY_FUNCTION__ << " azimAngle = " << azimAngle
              << std::endl;
//leave current loop if incl angle out of measurement area -15° --> +15.0°
//if( (inclAngle < -15.0) && (inclAngle > 15.0) ) 
  if(inclAngle < -15.0)
     {
      std::cout << "CONTINUEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << std::endl; 
      indices[i] = -1;
     continue;
   }
    // 1: calculate incoming azimuth = ROW index of indexMap
    // 2: calculate incoming inclination = COLUMN of indexMap
    unsigned int row = 0;
    unsigned int column = 0;
    returnRayIndex(azimAngle, inclAngle, &row, &column);
    std::cout << __PRETTY_FUNCTION__ << " row = " << row << std::endl;
    std::cout << __PRETTY_FUNCTION__ << " column = " << column << std::endl;
    // map column from sensor model to Velodyne firing sequence (order of
    // vertical rays differs between sensormodel and velodyne ros input)
    unsigned int columnMapped = lookupIndex(column);
    std::cout << __PRETTY_FUNCTION__ << " columnMapped = " << columnMapped
              << std::endl;

    //    push value into int* that is in indexMap[row][column]
    std::cout << __PRETTY_FUNCTION__
              << "_indexMap = " << _indexMap[row][columnMapped] << std::endl;

    // probe: index ausrechnen
    unsigned int idxCheck = columnMapped + 16 * row;
    std::cout << __PRETTY_FUNCTION__ << " idxCheck = " << idxCheck << std::endl;

    if (inclAngle > -15.0) 
    {
      if (inclAngle < 15.0) 
      {
        indices[i] = _indexMap[row][columnMapped];
        std::cout << "end of loop, content of vector indices = " << indices[i]
                  << std::endl;
      }
      else
      {
        indices[i] = -1;
      }
    }
      else
      {
        indices[i] = -1;
      }
      
    std::cout << "laufvariable = " << i << std::endl;
  }
  //abort();
}
} //namespace obvious